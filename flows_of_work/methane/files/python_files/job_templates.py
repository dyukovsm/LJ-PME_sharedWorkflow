import math
import mbuild as mb
from foyer import Forcefield
import foyer
import pandas as pd
import numpy as np
import random
import time
from parmed import residue
import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
import parmed
from parmed import gromacs
import signac
from flow import FlowProject, aggregator
from flow.environment import DefaultSlurmEnvironment
import flow
import subprocess
import os
from jinja2 import Environment, FileSystemLoader
import shutil
from files.python_files import names, job_tester
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from collections import defaultdict
from pathlib import Path

#GMX_PREFIX = 'gmx' # for grid
GMX_PREFIX = '/usr/local/gromacs/bin/gmx' # for potoff cluster


# Function to map continuous z to bin index
def z_to_bin_index(z_bins, z_val):
    """Map z-coordinate to nearest bin index."""
    bin_indices = np.abs(z_bins - z_val).argmin()
    return bin_indices

def simple_mdp_writer(job,mdp_name,parameters,constraints=None,templates_dir=None,template_name=None):
    loader = FileSystemLoader('.')
    env = Environment(loader=loader)
    path = os.path.relpath(f'{templates_dir}')
    MDP_NAME = template_name
    
    if constraints is None:
        update_dict = {
            'constraints_string' : ';',
            'constraints' : 'whatever',
            'constraint_algorithm_string' : ';',
            'constraint_algorithm' : 'whatever',
            'lincs_order_string' : ';',
            'lincs_order' : 'whatever'
        }
    elif 'lincs' in constraints:
        update_dict = {
            'constraints_string' : 'constraints         = ',
            'constraints' : 'all-bonds',
            'constraint_algorithm_string' : 'constraint-algorithm = ',
            'constraint_algorithm' : 'LINCS',
            'lincs_order_string' : 'lincs-order           = ',
            'lincs_order' : '6'
        }
    elif 'shake' in constraints:
        update_dict = {
            'constraints_string' : 'constraints         = ',
            'constraints' : 'all-angles',
            'constraint_algorithm_string' : 'constraint-algorithm = ',
            'constraint_algorithm' : 'SHAKE',
            'lincs_order_string' : 'shake-tol           = ',
            'lincs_order' : '0.00001'
        }
    parameters.update(update_dict)
    
    template_data = parameters
    template = env.get_template(f'{path}/{MDP_NAME}')
    
    output = template.render(template_data)
    with open(f'workspace/{job}/{mdp_name}','w') as f:
        f.write(output)
        
    
def gimme_dir(job):
    current_dir = os.getcwd()
    job_dir = f'{current_dir}/workspace/{job}' 
    return current_dir, job


def write_gmxINDEX_forRESIDUES(job, top_file = 'init.top', gro_file = 'init.gro', index_file_name = 'whacky_index_file.ndx'):

    #system_pmdTop = gromacs.GromacsTopologyFile('init.top')
    #gmx_gro = gromacs.GromacsGroFile.parse(f'init.gro')
    with(job):
        system_pmdTop = gromacs.GromacsTopologyFile(f'{top_file}')
        gmx_gro = gromacs.GromacsGroFile.parse(f'{gro_file}')
        system_pmdTop.box = gmx_gro.box
        system_pmdTop.positions = gmx_gro.positions

        angles4Gromacs = open(f'{index_file_name}','w')
        angles4Gromacs.write('[ WAT ] ;index1, atom_type\n')
        some_angles_written = False
        for i in system_pmdTop.residues:
            comments = [] 
            for j in i.atoms:
                correct_index = j.idx + 1
                comments.append(j)
                angles4Gromacs.write(f'{correct_index}\t')
                some_angles_written = True
            angles4Gromacs.write(f' \t;\t{str(comments)}\n')

        if some_angles_written:
            angles4Gromacs.write('\n;index file written correctly \n')
        angles4Gromacs.close() 
        
        
def manual_gmx_index_file_make(job,gro_file = 'init.gro', index_file_name = 'whacky_index_file.ndx',skip_residues_from_ncompounds=1000):
    with(job):
        skip_guess = skip_residues_from_ncompounds
        #skip_guess = math.ceil(math.log10(skip_guess))
        skip_guess = len(str(skip_guess))

        # Open the file and read the third line
        with open(f"{gro_file}", 'r') as f:
            for _ in range(2):
                next(f)  # skip first two lines
            line = f.readline()

        # Determine the column widths
        end_positions = [i for i, char in enumerate(line) if char != ' ' and (i == len(line) - 1 or line[i+1] == ' ')]
        column_widths = [end_positions[0] + 1] + [end_positions[i] - end_positions[i-1] for i in range(1, len(end_positions))]

        # Use numpy's genfromtxt to read the data with the determined column widths
        data = np.genfromtxt(f"{gro_file}", dtype=None, skip_header=2, delimiter=column_widths, encoding='utf-8')

        data=data[:-1]
        print(data)
        index_column = data['f2']

        result_dict = defaultdict(list)


        for record in data:
            #if record['f0'] == 
            key = record['f0'].strip()
            value = record['f2']#.strip()

            result_dict[key].append(value)

        index_file = open(f'{index_file_name}','w')
        header_preper = record['f0'].strip()#[0]
        header_preper = header_preper[skip_guess:-1]
        index_file.write(f'[ {header_preper} ]\n')

        #print('________________________________________________')
        #for i in range(10):
        #    print(header_preper)

        #sprint('________________________________________________')

        for i in result_dict.keys():
            #print(result_dict[i]) 
            dummy_list = result_dict[i]
            for j in dummy_list:
                index_file.write(f'{j}\t')
            index_file.write(f'\t ; {i} \n')

        index_file.close()


###############################################################################
######################## VORONOI TESSELLATION FUNCTIONS #######################
###############################################################################

def parse_gro_file(gro_file, residue_name):
    """
    Parse GROMACS .gro file and extract atoms for specified residue.

    Returns:
        molecules_dict: {resid: [(atom_name, x, y, z), ...], ...}
        box_vectors: (Lx, Ly, Lz) in nm
    """
    molecules_dict = {}
    box_vectors = None

    with open(gro_file, 'r') as f:
        lines = f.readlines()

    # Skip first 2 lines (title and atom count)
    data_lines = lines[2:]

    # Find box vectors: last line with exactly 3 floats containing "."
    for line in reversed(data_lines):
        parts = line.split()
        if len(parts) == 3:
            try:
                # Check if all parts are floats and contain "."
                if all('.' in p for p in parts):
                    box_vectors = (float(parts[0]), float(parts[1]), float(parts[2]))
                    break
            except ValueError:
                continue

    # Parse atom data (all lines except last box line)
    for line in data_lines[:-1]:
        if len(line) < 44:
            continue

        # Parse using fixed-width format
        resid = int(line[0:5].strip())
        resname = line[5:10].strip()
        atomname = line[10:15].strip()
        atomnr = int(line[15:20].strip())
        x = float(line[20:28])
        y = float(line[28:36])
        z = float(line[36:44])
        # Ignore velocities (columns 44+)

        # Only keep atoms matching the specified residue name
        if resname == residue_name:
            if resid not in molecules_dict:
                molecules_dict[resid] = []
            molecules_dict[resid].append((atomname, x, y, z))

    return molecules_dict, box_vectors


def calculate_center_of_mass(atoms_list, atom_mass_dict):
    """
    Calculate center of mass for a molecule.

    Args:
        atoms_list: [(atom_name, x, y, z), ...]
        atom_mass_dict: {'O': 15.999, 'H': 1.008, ...}

    Returns:
        (com_x, com_y, com_z) in nm
    """
    total_mass = 0.0
    com_x = 0.0
    com_y = 0.0
    com_z = 0.0

    for atom_name, x, y, z in atoms_list:
        mass = atom_mass_dict[atom_name]
        total_mass += mass
        com_x += mass * x
        com_y += mass * y
        com_z += mass * z

    com_x /= total_mass
    com_y /= total_mass
    com_z /= total_mass

    return (com_x, com_y, com_z)


def extract_all_coms(gro_file, residue_name, atom_mass_dict):
    """
    Extract centers of mass for all molecules in a .gro file.

    Returns:
        com_array: numpy array of shape (N, 3) with COM coordinates
        resid_list: list of residue IDs
        resname_list: list of residue names
        box_vectors: (Lx, Ly, Lz)
    """
    molecules_dict, box_vectors = parse_gro_file(gro_file, residue_name)

    com_list = []
    resid_list = []
    resname_list = []

    for resid in sorted(molecules_dict.keys()):
        atoms = molecules_dict[resid]
        com = calculate_center_of_mass(atoms, atom_mass_dict)
        com_list.append(com)
        resid_list.append(resid)
        resname_list.append(residue_name)

    com_array = np.array(com_list)

    return com_array, resid_list, resname_list, box_vectors


def calculate_voronoi_volumes_freud(gro_file, residue_name, atom_mass_dict):
    """
    Calculate Voronoi volumes using freud library with PBC.

    Args:
        gro_file: path to .gro file
        residue_name: name of residue to analyze (e.g., 'WAT')
        atom_mass_dict: dictionary mapping atom names to masses

    Returns:
        resid_list: list of residue IDs
        resname_list: list of residue names
        volumes: numpy array of Voronoi volumes (nm^3)
        coms: numpy array of centers of mass
    """
    import freud

    # Extract COMs
    coms, resid_list, resname_list, box_vectors = extract_all_coms(gro_file, residue_name, atom_mass_dict)

    # Create freud box (Lx, Ly, Lz, 0, 0, 0 for orthorhombic)
    box = freud.box.Box(Lx=box_vectors[0], Ly=box_vectors[1], Lz=box_vectors[2])

    # Compute Voronoi tessellation with PBC
    voro = freud.locality.Voronoi()
    voro.compute(system=(box, coms))

    # Get volumes
    volumes = voro.volumes

    return resid_list, resname_list, volumes, coms


def gmx_density_profile(job, trr_or_gro, index_file, tpr_file, output_xvg_name, first_frame, last_frame, slices=128):
    'return a profile of density. remember to use gpu node if tpr was also made on gpu'
    with(job):
        #comman_string = str(f'{GMX_PREFIX}') + str(' -f density ') + str(f'{trr_or_gro}') + str(' -n ') + str(f'{index_file}') + str(' -s ') + str(f'{tpr_file}') + str(' -o ') + str(f'{output_xvg_name}') + str(' -sl ')# + str(f'{slices})
        #comman_string = comman_string + str(f'{slices})
        #subprocess.run(comman=str(f'{GMX_PREFIX}') + str(' -f density ') + f'{trr_or_gro}' + ' -n ' + f'{index_file}' + ' -s ' + f'{tpr_file}', ' -o ', f'{output_xvg_name}', ' -sl ', f'{slices}'),shell=True)
        #f"{var} {var} text {var}"
        subprocess.run((f'{GMX_PREFIX}') + str(' density -f ') + str(f'{trr_or_gro}') + str(' -n ') + str(f'{index_file}') + str(' -s ') + str(f'{tpr_file}') + str(' -o ') + str(f'{output_xvg_name}') + str(' -sl ') + str(f'{slices}'),shell=True)
        #subprocess.run()#str(f'{GMX_PREFIX}') + str(' -f density ') + f'{trr_or_gro}' + ' -n ' + f'{index_file}' + ' -s ' + f'{tpr_file}', ' -o ', f'{output_xvg_name}', ' -sl ', f'{slices}'),shell=True)

        #p = subprocess.Popen([f'{GMX_PREFIX}', '-f', 'density', f'{trr_or_gro}', f'-n', f'{index_file}', f'-s', f'{tpr_file}', '-o', f'{output_xvg_name}', '-sl', f'{slices}'], stdin=subprocess.PIPE,stdout=subprocess.PIPE, cwd=os.getcwd())
        #out,err = p.communicate(input=str_sorrounding)
        #capture = p.decode() 
        

def build_slab_from_template(job, template_file_trr, template_file_tpr, template_file_mdp, 
                             output_name, pick_randomTrue_pick_allFalse=True):
        #with(job):
        
        # === FILE PATHS ===
        #base_path = Path("../../files/coordinates/slab_template/CUT")
        #base_path = Path("./")
        mdp_file = f"{template_file_mdp}" # mdp_file = base_path / "PRO_SURFTEN.mdp"
        tpr = f"{template_file_tpr}" # tpr = base_path / "PRO_SURFTEN_CHUNK_1.tpr"
        trr = f"{template_file_trr}" # trr = base_path / "PRO_SURFTEN_CHUNK_1.trr"
        output_gro = f"{output_name}" # output_gro = "random_frame.gro"

        # === PARSE MDP PARAMETERS (fast method) ===
        params = {}
        with open(mdp_file) as f:
            for line in f:
                if "=" not in line:
                    continue
                key, val = line.split("=", 1)
                key = key.strip()
                val = val.strip().split()[0]  # first token after '='
                if key in {"dt", "nsteps", "nstxout", "nstvout"}:
                    params[key] = float(val)

        dt = params["dt"]
        nsteps = int(params["nsteps"])
        nstxout = int(params["nstxout"])

        # === COMPUTE FRAME TIMES ===
        frame_interval = nstxout * dt
        max_time = nsteps * dt
        valid_times = [i * frame_interval for i in range(0, int(nsteps / nstxout) + 1)]
        valid_times[-1] = max_time  # ensure exact last time

        # === EXECUTION ===
        if pick_randomTrue_pick_allFalse:
            random_time = random.choice(valid_times)
            print(f"Extracting random frame at t = {random_time:.3f} ps")
            subprocess.run([
                "/usr/local/gromacs/bin/gmx", "trjconv",
                "-s", str(tpr),
                "-f", str(trr),
                "-o", output_gro,
                "-dump", str(random_time)
            ], input="0\n", text=True, check=True)
        else:
            print(f"Extracting all {len(valid_times)} frames individually...")
            for i, t in enumerate(valid_times):
                out_name = f"{output_gro}_frame_{i:05d}.gro"
                print(f"  → dumping t = {t:.3f} ps to {out_name}")
                subprocess.run([
                    "/usr/local/gromacs/bin/gmx", "trjconv",
                    "-s", str(tpr),
                    "-f", str(trr),
                    "-o", out_name,
                    "-dump", str(t)
                ], input="0\n", text=True, check=True)


###############################################################################
##################### DENSITY PROFILE AND SLAB ALIGNMENT #######################
###############################################################################

def align_slab_to_center(job, tpr_file, trr_file, gro_file, output_trr, output_gro, bin_size):
    """
    Align liquid slab to box center.

    This function:
    1. Computes density profile along z-axis
    2. Finds interfaces via derivative analysis
    3. Shifts slab so liquid phase is centered at Z_LEN/2

    Args:
        job: signac job context
        tpr_file: input TPR file
        trr_file: input TRR trajectory
        gro_file: input GRO file
        output_trr: output aligned TRR filename
        output_gro: output aligned GRO filename
        bin_size: bin size for density profile (nm)

    Returns:
        dict with alignment data (shift vector, interface positions, etc.)
    """
    import matplotlib.pyplot as plt

    with(job):
        # Get box dimensions
        dummy_mb = mb.load(gro_file)
        Z_LEN = dummy_mb.box.Lz
        tanh_bin_count = int(np.round(Z_LEN/bin_size, 0))

        def parse_density_xvg(xvg_file):
            """Parse XVG file and return z_bins and density arrays."""
            z_bins, density = [], []
            with open(xvg_file, 'r') as f:
                for line in f:
                    if not (line.startswith('#') or line.startswith('@')):
                        parts = line.split()
                        if len(parts) >= 2:
                            z_bins.append(float(parts[0]))
                            density.append(float(parts[1]))
            return np.array(z_bins), np.array(density)

        # Step 1: Compute initial density profile
        cmd = f"echo 0 | {names.GMX_PREFIX} density -s {tpr_file} -f {trr_file} -o {names.NAME_DENSITY_XVG} -d Z -sl {tanh_bin_count}"
        subprocess.run(cmd, shell=True, capture_output=True, text=True)

        z_current, density_current = parse_density_xvg(names.NAME_DENSITY_XVG)

        # Step 2: Find interfaces via derivative
        ddens_dz = np.gradient(density_current, z_current)
        max_idx = np.argmax(ddens_dz)
        min_idx = np.argmin(ddens_dz)

        z_left_interface = z_current[max_idx]
        z_right_interface = z_current[min_idx]
        z_midpoint = Z_LEN / 2.0
        z_mean_interface = (z_left_interface + z_right_interface) / 2.0
        middle2zOver2vector = z_midpoint - z_mean_interface

        # Store pre-shift data
        z_pre = z_current.copy()
        density_pre = density_current.copy()
        ddens_dz_pre = ddens_dz.copy()

        # Step 3: Shift slab to center
        cmd_shift = f"echo 0 | {names.GMX_PREFIX} trjconv -f {trr_file} -s {tpr_file} -o {output_trr} -trans 0.0 0.0 {middle2zOver2vector} -pbc mol"
        subprocess.run(cmd_shift, shell=True, capture_output=True, text=True)

        cmd_shift_gro = f"echo 0 | {names.GMX_PREFIX} trjconv -f {gro_file} -s {tpr_file} -o {output_gro} -trans 0.0 0.0 {middle2zOver2vector} -pbc mol"
        subprocess.run(cmd_shift_gro, shell=True, capture_output=True, text=True)

        # Step 4: Compute post-shift density profile
        cmd_post = f"echo 0 | {names.GMX_PREFIX} density -s {tpr_file} -f {output_trr} -o {names.NAME_DENSITY_XVG_POST} -d Z -sl {tanh_bin_count}"
        subprocess.run(cmd_post, shell=True, capture_output=True, text=True)

        z_post, density_post = parse_density_xvg(names.NAME_DENSITY_XVG_POST)
        ddens_dz_post = np.gradient(density_post, z_post)

        max_idx_post = np.argmax(ddens_dz_post)
        min_idx_post = np.argmin(ddens_dz_post)
        z_left_interface_post = z_post[max_idx_post]
        z_right_interface_post = z_post[min_idx_post]
        z_mean_interface_post = (z_left_interface_post + z_right_interface_post) / 2.0

        # Generate alignment comparison plot
        fig, ((ax1, ax3), (ax2, ax4)) = plt.subplots(2, 2, figsize=(20, 10), sharex=True, sharey='row')

        # PRE-SHIFT plots
        ax1.plot(z_pre, density_pre, 'b-', linewidth=2, label='Density')
        ax1.axvline(z_left_interface, color='red', linestyle='--', linewidth=2,
                   label=f'Left interface: {z_left_interface:.3f} nm')
        ax1.axvline(z_right_interface, color='green', linestyle='--', linewidth=2,
                   label=f'Right interface: {z_right_interface:.3f} nm')
        ax1.axvline(z_midpoint, color='black', linestyle='-', linewidth=2,
                   label=f'Z_LEN/2: {z_midpoint:.3f} nm')
        ax1.axvline(z_mean_interface, color='purple', linestyle='--', linewidth=2,
                   label=f'Mean interface: {z_mean_interface:.3f} nm')
        ax1.set_ylabel('Density (kg/m³)', fontsize=12, fontweight='bold')
        ax1.set_title('PRE-SHIFT Density Profile', fontsize=14, fontweight='bold')
        ax1.legend(fontsize=8)
        ax1.grid(True, alpha=0.3)

        ax2.plot(z_pre, ddens_dz_pre, 'k-', linewidth=2, label='dρ/dz')
        ax2.axhline(0, color='gray', linestyle=':', linewidth=1)
        ax2.set_xlabel('z (nm)', fontsize=12, fontweight='bold')
        ax2.set_ylabel('dρ/dz (kg/m⁴)', fontsize=12, fontweight='bold')
        ax2.set_title('PRE-SHIFT Derivative', fontsize=14, fontweight='bold')
        ax2.legend(fontsize=8)
        ax2.grid(True, alpha=0.3)

        # POST-SHIFT plots
        ax3.plot(z_post, density_post, 'b-', linewidth=2, label='Density')
        ax3.axvline(z_left_interface_post, color='red', linestyle='--', linewidth=2,
                   label=f'Left interface: {z_left_interface_post:.3f} nm')
        ax3.axvline(z_right_interface_post, color='green', linestyle='--', linewidth=2,
                   label=f'Right interface: {z_right_interface_post:.3f} nm')
        ax3.axvline(z_midpoint, color='black', linestyle='-', linewidth=2,
                   label=f'Z_LEN/2: {z_midpoint:.3f} nm')
        ax3.axvline(z_mean_interface_post, color='purple', linestyle='--', linewidth=2,
                   label=f'Mean interface: {z_mean_interface_post:.3f} nm')
        ax3.set_ylabel('Density (kg/m³)', fontsize=12, fontweight='bold')
        ax3.set_title('POST-SHIFT Density Profile', fontsize=14, fontweight='bold')
        ax3.legend(fontsize=8)
        ax3.grid(True, alpha=0.3)

        ax4.plot(z_post, ddens_dz_post, 'k-', linewidth=2, label='dρ/dz')
        ax4.axhline(0, color='gray', linestyle=':', linewidth=1)
        ax4.set_xlabel('z (nm)', fontsize=12, fontweight='bold')
        ax4.set_ylabel('dρ/dz (kg/m⁴)', fontsize=12, fontweight='bold')
        ax4.set_title('POST-SHIFT Derivative', fontsize=14, fontweight='bold')
        ax4.legend(fontsize=8)
        ax4.grid(True, alpha=0.3)

        plt.tight_layout()
        plt.savefig(names.NAME_SHIFT_COMPARISON_PNG, dpi=300)
        plt.close()

        # Return alignment data
        alignment_data = {
            'shift_vector': middle2zOver2vector,
            'z_len': Z_LEN,
            'z_left_interface_post': z_left_interface_post,
            'z_right_interface_post': z_right_interface_post,
            'z_mean_interface_post': z_mean_interface_post
        }

        print(f"Slab alignment complete. Shift vector: {middle2zOver2vector:.3f} nm")

        return alignment_data


def classify_phases_from_aligned_slab(job, tpr_file, aligned_trr, aligned_gro, bin_size, d_multiplier, noise_multiplier):
    """
    Robustly classify liquid, gas, and transition phases from aligned slab.

    This function:
    1. Computes density profile from aligned trajectory
    2. Fits tanh function to get interface parameters
    3. Robustly classifies phases using derivative-based noise analysis
    4. Saves comprehensive phase data to file

    Args:
        job: signac job context
        tpr_file: input TPR file
        aligned_trr: aligned TRR trajectory
        aligned_gro: aligned GRO file
        bin_size: bin size for density profile (nm)
        d_multiplier: multiplier for interface thickness
        noise_multiplier: multiplier for derivative noise threshold

    Returns:
        dict with phase classification data
    """
    from scipy.optimize import curve_fit
    import matplotlib.pyplot as plt

    with(job):
        # Get box dimensions
        dummy_mb = mb.load(aligned_gro)
        Z_LEN = dummy_mb.box.Lz
        tanh_bin_count = int(np.round(Z_LEN/bin_size, 0))

        # Tanh fit function
        def tanh_fit_function(z, rho_liq, rho_gas, z0, d):
            return 0.5 * (rho_liq + rho_gas) - 0.5 * (rho_liq - rho_gas) * np.tanh(2.0 * (z - z0) / d)

        def parse_density_xvg(xvg_file):
            """Parse XVG file and return z_bins and density arrays."""
            z_bins, density = [], []
            with open(xvg_file, 'r') as f:
                for line in f:
                    if not (line.startswith('#') or line.startswith('@')):
                        parts = line.split()
                        if len(parts) >= 2:
                            z_bins.append(float(parts[0]))
                            density.append(float(parts[1]))
            return np.array(z_bins), np.array(density)

        # Compute density profile from aligned trajectory
        cmd = f"echo 0 | {names.GMX_PREFIX} density -s {tpr_file} -f {aligned_trr} -o {names.NAME_DENSITY_XVG_POST} -d Z -sl {tanh_bin_count}"
        subprocess.run(cmd, shell=True, capture_output=True, text=True)

        z_post, density_post = parse_density_xvg(names.NAME_DENSITY_XVG_POST)
        ddens_dz_post = np.gradient(density_post, z_post)

        max_idx_post = np.argmax(ddens_dz_post)
        min_idx_post = np.argmin(ddens_dz_post)
        z_left_interface_post = z_post[max_idx_post]
        z_right_interface_post = z_post[min_idx_post]
        z_mean_interface_post = (z_left_interface_post + z_right_interface_post) / 2.0

        # Step 5: Fit tanh (rotate data for better fit)
        z_rotated = z_post.copy()
        density_rotated = density_post.copy()

        # Rotate: flip left half to right
        for i in range(len(z_post)):
            if z_post[i] < Z_LEN / 2.0:
                z_rotated[i] = Z_LEN - z_post[i]

        # Sort by z after rotation
        sort_idx = np.argsort(z_rotated)
        z_rotated = z_rotated[sort_idx]
        density_rotated = density_rotated[sort_idx]

        # Fit tanh
        popt, pcov = curve_fit(tanh_fit_function, z_rotated, density_rotated,
                              p0=[np.max(density_rotated), np.min(density_rotated), z_right_interface_post, 1.0],
                              maxfev=10000)
        rho_liq_fit, rho_gas_fit, z0_fit, d_fit = popt

        # Unrotate to get both interfaces
        left_z0_fit = Z_LEN - z0_fit
        right_z0_fit = z0_fit

        left_lower_limit = left_z0_fit - d_multiplier * d_fit
        left_upper_limit = left_z0_fit + d_multiplier * d_fit
        right_lower_limit = right_z0_fit - d_multiplier * d_fit
        right_upper_limit = right_z0_fit + d_multiplier * d_fit

        # Step 6: ROBUST PHASE CLASSIFICATION USING DERIVATIVE ANALYSIS

        # Compute absolute derivative
        abs_ddens_dz = np.abs(ddens_dz_post)

        # Find baseline noise level: use median of smallest 50% of |dρ/dz|
        sorted_abs_deriv = np.sort(abs_ddens_dz)
        baseline_noise = np.median(sorted_abs_deriv[:len(sorted_abs_deriv)//2])

        # Define transition threshold
        transition_threshold = baseline_noise * noise_multiplier

        # Classify each z-point
        is_transition = abs_ddens_dz > transition_threshold
        is_stable = ~is_transition

        # Among stable regions, classify as liquid or gas based on proximity to max/min
        if np.sum(is_stable) > 0:
            stable_densities = density_post[is_stable]
            max_stable_density = np.max(stable_densities)
            min_stable_density = np.min(stable_densities)

            # Classify based on distance from max vs min
            # Liquid: closer to max than to min
            # Gas: closer to min than to max
            dist_from_max = max_stable_density - density_post
            dist_from_min = density_post - min_stable_density

            is_liquid = is_stable & (dist_from_max < dist_from_min)
            is_gas = is_stable & (dist_from_min < dist_from_max)
        else:
            # Fallback: no stable regions found
            max_density = np.max(density_post)
            min_density = np.min(density_post)
            dist_from_max = max_density - density_post
            dist_from_min = density_post - min_density
            is_liquid = dist_from_max < dist_from_min
            is_gas = dist_from_min < dist_from_max

        # Calculate phase densities using derivative-based classification
        # These are "robust" because they use stable regions only
        if np.sum(is_liquid) > 0:
            rho_liquid_robust = np.mean(density_post[is_liquid])
            z_liquid_coords = z_post[is_liquid]
        else:
            rho_liquid_robust = rho_liq_fit
            z_liquid_coords = np.array([])

        if np.sum(is_gas) > 0:
            rho_gas_robust = np.mean(density_post[is_gas])
            z_gas_coords = z_post[is_gas]
        else:
            rho_gas_robust = rho_gas_fit
            z_gas_coords = np.array([])

        if np.sum(is_transition) > 0:
            rho_transition_robust = np.mean(density_post[is_transition])
            z_transition_coords = z_post[is_transition]
        else:
            rho_transition_robust = 0.5 * (rho_liquid_robust + rho_gas_robust)
            z_transition_coords = np.array([])

        # Step 7: Save comprehensive phase data
        with open(names.NAME_PHASE_DATA, 'w') as f:
            f.write("=" * 80 + "\n")
            f.write("SLAB ALIGNMENT AND PHASE CLASSIFICATION DATA\n")
            f.write("=" * 80 + "\n\n")

            f.write("--- MIN-MAX ---\n")
            max_density = np.max(density_post)
            min_density = np.min(density_post)
            f.write(f"minDens_raw (kg/m³): {max_density}\n")
            f.write(f"maxDens_raw (kg/m³): {min_density}\n\n")

            f.write("--- TANH FIT PARAMETERS ---\n")
            f.write(f"rho_liq_fit (kg/m³):     {rho_liq_fit:.10f}\n")
            f.write(f"rho_gas_fit (kg/m³):     {rho_gas_fit:.10f}\n")
            f.write(f"z0_fit (nm):             {z0_fit:.10f}\n")
            f.write(f"d_fit (nm):              {d_fit:.10f}\n\n")

            f.write("--- TRANSITION ZONE BOUNDARIES ---\n")
            f.write(f"Left transition start (nm):   {left_lower_limit:.10f}\n")
            f.write(f"Left transition stop (nm):    {left_upper_limit:.10f}\n")
            f.write(f"Right transition start (nm):  {right_lower_limit:.10f}\n")
            f.write(f"Right transition stop (nm):   {right_upper_limit:.10f}\n\n")

            f.write("--- ROBUST PHASE DENSITIES (Derivative-Based Classification) ---\n")
            f.write(f"Baseline noise (|dρ/dz| median):  {baseline_noise:.10f} kg/m⁴\n")
            f.write(f"Transition threshold:              {transition_threshold:.10f} kg/m⁴\n")
            f.write(f"Liquid density (robust mean):      {rho_liquid_robust:.10f} kg/m³\n")
            f.write(f"Gas density (robust mean):         {rho_gas_robust:.10f} kg/m³\n")
            f.write(f"Transition density (robust mean):  {rho_transition_robust:.10f} kg/m³\n")
            f.write(f"Num liquid points:                 {np.sum(is_liquid)}\n")
            f.write(f"Num gas points:                    {np.sum(is_gas)}\n")
            f.write(f"Num transition points:             {np.sum(is_transition)}\n\n")

            f.write("--- Z-COORDINATES BY PHASE ---\n")
            f.write(f"Liquid z-coords (nm): {z_liquid_coords.tolist()}\n\n")
            f.write(f"Gas z-coords (nm): {z_gas_coords.tolist()}\n\n")
            f.write(f"Transition z-coords (nm): {z_transition_coords.tolist()}\n\n")

        # Convert z-coordinates from nm to Å for consistency with MDAnalysis (which uses Å)
        z_post_angstrom = z_post * 10.0
        bin_size_angstrom = bin_size * 10.0

        # Save phase classification data as NPZ for volume classification
        # Store: phase masks (boolean arrays), z-bin positions (in Å), and bin edges
        np.savez_compressed(
            names.NAME_PHASE_Z_COORDS,
            z_bins=z_post_angstrom,           # Z-bin center positions in Å
            is_liquid=is_liquid,              # Boolean mask for liquid bins
            is_gas=is_gas,                    # Boolean mask for gas bins
            is_transition=is_transition,      # Boolean mask for transition bins
            bin_size=bin_size_angstrom        # Bin size in Å for mapping continuous z to bins
        )

        print(f"Saved phase z-coordinates (converted to Å):")
        print(f"  z-range: {z_post_angstrom.min():.1f} - {z_post_angstrom.max():.1f} Å")
        print(f"  bin_size: {bin_size_angstrom:.2f} Å")

        # Step 8: Generate phase classification plots
        fig, (ax3, ax4) = plt.subplots(2, 1, figsize=(20, 10), sharex=True)

        # Get density range for vertical shading
        y_min_ax3, y_max_ax3 = density_post.min(), density_post.max()
        y_range_ax3 = y_max_ax3 - y_min_ax3
        y_min_ax3 -= 0.05 * y_range_ax3
        y_max_ax3 += 0.05 * y_range_ax3

        # POST-SHIFT plots with phase classification
        # Add vertical grey bands for TRANSITION bins only (shade each bin individually)
        if np.sum(is_transition) > 0:
            tra_z = z_post[is_transition]
            half_bin = bin_size / 2.0
            # Shade each transition bin individually
            for i, z_center in enumerate(tra_z):
                label = 'Transition region' if i == 0 else None
                ax3.fill_betweenx([y_min_ax3, y_max_ax3],
                                 z_center - half_bin, z_center + half_bin,
                                 color='lightgrey', alpha=0.3, label=label)

        ax3.plot(z_post, density_post, 'b-', linewidth=2, label='Density')
        ax3.scatter(z_post[is_liquid], density_post[is_liquid], c='blue', s=40, alpha=0.6, label='Liquid')
        ax3.scatter(z_post[is_gas], density_post[is_gas], c='orange', s=40, alpha=0.6, label='Gas')
        ax3.scatter(z_post[is_transition], density_post[is_transition], c='red', s=40, alpha=0.6, label='Transition')
        ax3.axhline(rho_liquid_robust, color='darkblue', linestyle=':', linewidth=2,
                   label=f'ρ_liq (robust): {rho_liquid_robust:.1f} kg/m³')
        ax3.axhline(rho_gas_robust, color='darkorange', linestyle=':', linewidth=2,
                   label=f'ρ_gas (robust): {rho_gas_robust:.1f} kg/m³')
        ax3.set_ylabel('Density (kg/m³)', fontsize=12, fontweight='bold')
        ax3.set_title('POST-SHIFT Density with Phase Classification', fontsize=14, fontweight='bold')
        ax3.set_ylim(y_min_ax3, y_max_ax3)
        ax3.legend(fontsize=8)
        ax3.grid(True, alpha=0.3)

        # Get derivative range for vertical shading
        y_min_ax4, y_max_ax4 = ddens_dz_post.min(), ddens_dz_post.max()
        y_range_ax4 = y_max_ax4 - y_min_ax4
        y_min_ax4 -= 0.05 * y_range_ax4
        y_max_ax4 += 0.05 * y_range_ax4

        # Add vertical grey bands for TRANSITION bins only (shade each bin individually)
        if np.sum(is_transition) > 0:
            for z_center in tra_z:
                ax4.fill_betweenx([y_min_ax4, y_max_ax4],
                                 z_center - half_bin, z_center + half_bin,
                                 color='lightgrey', alpha=0.3)

        ax4.plot(z_post, ddens_dz_post, 'k-', linewidth=2, label='dρ/dz')
        ax4.axhline(transition_threshold, color='red', linestyle='--', linewidth=1,
                   label=f'Transition threshold: {transition_threshold:.2f}')
        ax4.axhline(-transition_threshold, color='red', linestyle='--', linewidth=1)
        ax4.axhline(0, color='gray', linestyle=':', linewidth=1)
        ax4.fill_between(z_post, -transition_threshold, transition_threshold, color='lightgreen', alpha=0.2,
                        label='Stable region (horizontal)')
        ax4.set_xlabel('z (nm)', fontsize=12, fontweight='bold')
        ax4.set_ylabel('dρ/dz (kg/m⁴)', fontsize=12, fontweight='bold')
        ax4.set_title('POST-SHIFT Derivative with Noise Threshold', fontsize=14, fontweight='bold')
        ax4.set_ylim(y_min_ax4, y_max_ax4)
        ax4.legend(fontsize=8)
        ax4.grid(True, alpha=0.3)

        plt.tight_layout()
        plt.savefig(names.NAME_PHASE_CLASSIFICATION_PNG, dpi=300)
        plt.close()

        # Return phase data dictionary
        phase_data = {
            'rho_liq_fit': rho_liq_fit,
            'rho_gas_fit': rho_gas_fit,
            'rho_liquid_robust': rho_liquid_robust,
            'rho_gas_robust': rho_gas_robust,
            'rho_transition_robust': rho_transition_robust,
            'z0_fit': z0_fit,
            'd_fit': d_fit,
            'left_lower_limit': left_lower_limit,
            'left_upper_limit': left_upper_limit,
            'right_lower_limit': right_lower_limit,
            'right_upper_limit': right_upper_limit,
            'baseline_noise': baseline_noise,
            'transition_threshold': transition_threshold
        }

        print(f"Phase classification complete:")
        print(f"  ρ_liquid (robust) = {rho_liquid_robust:.2f} kg/m³")
        print(f"  ρ_gas (robust) = {rho_gas_robust:.2f} kg/m³")
        print(f"  ρ_transition (robust) = {rho_transition_robust:.2f} kg/m³")

        return phase_data


###############################################################################
######################## MDANALYSIS VORONOI FUNCTIONS #########################
###############################################################################

def compute_voronoi_volumes_with_mda(job, gro_file, trr_file, output_npz, max_frames=None):
    """
    Compute Voronoi volumes on ALL residues (full system tessellation only).

    CRITICAL: Voronoi tessellation MUST be computed on the full system to get correct
    volumes (each cell needs all its neighbors). Phase classification is done separately.

    Args:
        job: signac job context
        gro_file: path to .gro file
        trr_file: path to .trr trajectory file
        output_npz: output filename for ALL residues (volumes + COMs)
        max_frames: maximum number of frames to process (None for all)

    Returns:
        None (saves data to output_npz with com_trajectory and volumes)
    """
    import MDAnalysis as mda
    import freud

    with(job):
        print(f"Loading: {gro_file} and {trr_file}")
        u = mda.Universe(gro_file, trr_file)

        n_frames_total = len(u.trajectory)
        n_residues = len(u.residues)

        # Limit frames based on max_frames parameter
        if max_frames is None:
            n_frames = n_frames_total
            print(f"Found {n_frames_total} total frames. Will process all frames.")
        else:
            n_frames = min(n_frames_total, max_frames)
            print(f"Found {n_frames_total} total frames. Will process first {n_frames} frames.")

        print(f"Found {n_residues} residues.")

        # Initialize storage for ALL residues
        all_coms = np.zeros((n_frames, n_residues, 3), dtype=np.float32)
        all_volumes = np.zeros((n_frames, n_residues), dtype=np.float32)

        # Initialize freud
        voronoi = freud.locality.Voronoi()

        print("Computing Voronoi tessellation on FULL system (all residues)...")

        # Process trajectory
        for i, ts in enumerate(u.trajectory[:n_frames]):
            if i % 100 == 0:
                print(f"  Processing frame {i}/{n_frames}...")

            # Get box and COMs for ALL residues
            box = freud.Box.from_box(ts.dimensions)
            coms_this_frame = np.array([residue.atoms.center_of_mass() for residue in u.residues])

            # Compute Voronoi on ALL residues (critical for correct volumes)
            voronoi.compute(system=(box, coms_this_frame))

            # Store ALL COMs and volumes
            all_coms[i] = coms_this_frame
            all_volumes[i] = voronoi.volumes

        print(f"Tessellation complete. Saving full results to {output_npz}...")
        np.savez_compressed(
            output_npz,
            com_trajectory=all_coms,
            volumes=all_volumes
        )

        print(f"Voronoi tessellation complete: {n_frames} frames, {n_residues} residues")


def classify_voronoi_by_phase(job, input_npz_all, output_npz_liq, output_npz_gas, output_npz_tra, output_z_hist_png):
    """
    Classify pre-computed Voronoi volumes into phases based on z-coordinate binning.

    Uses phase classification data (z-bins and phase masks) to assign each residue's
    Voronoi volume to liquid, gas, or transition phase based on its COM z-coordinate.

    Args:
        job: signac job context
        input_npz_all: NPZ file with all Voronoi volumes and COMs
        output_npz_liq: output NPZ for liquid phase
        output_npz_gas: output NPZ for gas phase
        output_npz_tra: output NPZ for transition phase
        output_z_hist_png: output PNG showing z-position distribution by phase

    Returns:
        None (saves 3 phase-specific NPZ files + histogram)
    """
    import matplotlib.pyplot as plt

    with(job):
        # Load Voronoi data
        print(f"Loading Voronoi data from {input_npz_all}...")
        voronoi_data = np.load(input_npz_all)
        all_coms = voronoi_data['com_trajectory']
        all_volumes = voronoi_data['volumes']
        n_frames, n_residues, _ = all_coms.shape

        # Load phase classification data
        print("Loading phase classification data...")
        phase_data = np.load(names.NAME_PHASE_Z_COORDS)

        # Check if we have new format (z_bins + masks) or old format (z_liquid, z_gas, z_transition)
        if 'z_bins' in phase_data:
            # New format: use phase masks
            z_bins = phase_data['z_bins']
            is_liquid = phase_data['is_liquid']
            is_gas = phase_data['is_gas']
            is_transition = phase_data['is_transition']
            bin_size = float(phase_data['bin_size'])
            print(f"Using NEW phase classification format (z-bins + masks)")
        else:
            # Old format: need to regenerate - throw helpful error
            raise ValueError(
                f"ERROR: {names.NAME_PHASE_Z_COORDS} has OLD format!\n"
                f"Please rerun CLASSIFY_PHASES operation to regenerate with new format:\n"
                f"  python project.py run -o CLASSIFY_PHASES\n"
                f"Old format had: {list(phase_data.keys())}\n"
                f"New format needs: ['z_bins', 'is_liquid', 'is_gas', 'is_transition', 'bin_size']"
            )

        print(f"Phase bin statistics:")
        print(f"  Total bins: {len(z_bins)}")
        print(f"  Liquid bins: {np.sum(is_liquid)}")
        print(f"  Gas bins: {np.sum(is_gas)}")
        print(f"  Transition bins: {np.sum(is_transition)}")
        print(f"  Bin size: {bin_size} nm")



        # Classify residues by phase for each frame
        print("Classifying volumes by phase...")
        liquid_coms_list = []
        liquid_volumes_list = []
        gas_coms_list = []
        gas_volumes_list = []
        transition_coms_list = []
        transition_volumes_list = []

        # For z-position histogram
        all_z_positions = {'gas': [], 'transition': [], 'liquid': []}

        for frame_idx in range(n_frames):
            if frame_idx % 100 == 0:
                print(f"  Classifying frame {frame_idx}/{n_frames}...")

            # Get z-coordinates for all residues in this frame
            z_coords = all_coms[frame_idx, :, 2]

            # Map each residue to a bin and classify
            liquid_mask = np.zeros(n_residues, dtype=bool)
            gas_mask = np.zeros(n_residues, dtype=bool)
            transition_mask = np.zeros(n_residues, dtype=bool)

            for res_idx in range(n_residues):
                z_val = z_coords[res_idx]
                bin_idx = z_to_bin_index(z_bins, z_val)

                if is_liquid[bin_idx]:
                    liquid_mask[res_idx] = True
                    all_z_positions['liquid'].append(z_val)
                elif is_gas[bin_idx]:
                    gas_mask[res_idx] = True
                    all_z_positions['gas'].append(z_val)
                elif is_transition[bin_idx]:
                    transition_mask[res_idx] = True
                    all_z_positions['transition'].append(z_val)

            # Extract phase-specific data
            liquid_coms_list.append(all_coms[frame_idx][liquid_mask])
            liquid_volumes_list.append(all_volumes[frame_idx][liquid_mask])
            gas_coms_list.append(all_coms[frame_idx][gas_mask])
            gas_volumes_list.append(all_volumes[frame_idx][gas_mask])
            transition_coms_list.append(all_coms[frame_idx][transition_mask])
            transition_volumes_list.append(all_volumes[frame_idx][transition_mask])

        # Save phase-specific NPZ files
        for phase_name, coms_list, vols_list, output_file in [
            ('GAS', gas_coms_list, gas_volumes_list, output_npz_gas),
            ('TRANSITION', transition_coms_list, transition_volumes_list, output_npz_tra),
            ('LIQUID', liquid_coms_list, liquid_volumes_list, output_npz_liq)
        ]:
            max_phase_residues = max(len(vols) for vols in vols_list) if vols_list else 0

            if max_phase_residues > 0:
                phase_coms_padded = np.full((n_frames, max_phase_residues, 3), np.nan, dtype=np.float32)
                phase_volumes_padded = np.full((n_frames, max_phase_residues), np.nan, dtype=np.float32)

                for frame_idx, (coms, vols) in enumerate(zip(coms_list, vols_list)):
                    n_res = len(vols)
                    if n_res > 0:
                        phase_coms_padded[frame_idx, :n_res, :] = coms
                        phase_volumes_padded[frame_idx, :n_res] = vols

                np.savez_compressed(
                    output_file,
                    com_trajectory=phase_coms_padded,
                    volumes=phase_volumes_padded
                )

                total_residues = sum(len(vols) for vols in vols_list)
                print(f"  {phase_name}: {total_residues} total residue-frames → {output_file}")
            else:
                print(f"  {phase_name}: No residues found!")

        # Generate z-position histogram to verify classification
        print(f"Generating z-position histogram: {output_z_hist_png}")
        fig, ax = plt.subplots(1, 1, figsize=(12, 8))

        # Plot in order: GAS → TRANSITION → LIQUID
        colors = {'gas': 'orange', 'transition': 'red', 'liquid': 'blue'}
        for phase in ['gas', 'transition', 'liquid']:
            if len(all_z_positions[phase]) > 0:
                ax.hist(all_z_positions[phase], bins=100, alpha=0.6,
                       label=f'{phase.upper()} (n={len(all_z_positions[phase])})',
                       color=colors[phase])

        ax.set_xlabel('Z-coordinate (Å)', fontsize=12, fontweight='bold')
        ax.set_ylabel('Count', fontsize=12, fontweight='bold')
        ax.set_title('Z-Position Distribution by Phase Classification', fontsize=14, fontweight='bold')
        ax.legend(fontsize=10)
        ax.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig(output_z_hist_png, dpi=300)
        plt.close()

        # Calculate tessellation-based densities from Voronoi volumes
        # Collect all volumes per phase (flatten across frames)
        liq_volumes_all = np.concatenate([vols for vols in liquid_volumes_list if len(vols) > 0]) if any(len(v) > 0 for v in liquid_volumes_list) else np.array([])
        gas_volumes_all = np.concatenate([vols for vols in gas_volumes_list if len(vols) > 0]) if any(len(v) > 0 for v in gas_volumes_list) else np.array([])
        tra_volumes_all = np.concatenate([vols for vols in transition_volumes_list if len(vols) > 0]) if any(len(v) > 0 for v in transition_volumes_list) else np.array([])

        # Calculate densities: ρ = m / V
        # Volumes are in Å³, need to convert to m³ for kg/m³
        # molecule_mass_kg = molar_mass / Avogadro * 1e-3
        avogadro = 6.022e23
        molar_mass = names.CURRENT_MOLAR_MASS
        molecule_mass_kg = (molar_mass / avogadro) * 1e-3

        # Average volume per phase (in Å³)
        avg_vol_liq = np.mean(liq_volumes_all) if len(liq_volumes_all) > 0 else 0
        avg_vol_gas = np.mean(gas_volumes_all) if len(gas_volumes_all) > 0 else 0
        avg_vol_tra = np.mean(tra_volumes_all) if len(tra_volumes_all) > 0 else 0

        # Convert to density: ρ = m / V, where V is in Å³ = 1e-30 m³
        dens_Tess_Liq = molecule_mass_kg / (avg_vol_liq * 1e-30) if avg_vol_liq > 0 else 0
        dens_Tess_Gas = molecule_mass_kg / (avg_vol_gas * 1e-30) if avg_vol_gas > 0 else 0
        dens_Tess_Tra = molecule_mass_kg / (avg_vol_tra * 1e-30) if avg_vol_tra > 0 else 0

        print(f"Tessellation-based densities:")
        print(f"  Liquid: {dens_Tess_Liq:.2f} kg/m³ (avg vol: {avg_vol_liq:.2f} Å³)")
        print(f"  Gas: {dens_Tess_Gas:.2f} kg/m³ (avg vol: {avg_vol_gas:.2f} Å³)")
        print(f"  Transition: {dens_Tess_Tra:.2f} kg/m³ (avg vol: {avg_vol_tra:.2f} Å³)")

        # Append tessellation-based densities to the phase classification data file
        with open(names.NAME_PHASE_DATA, 'a') as f:
            f.write("--- TESSELLATION-BASED DENSITIES (Voronoi Volume Classification) ---\n")
            f.write(f"dens_Tess_Liq (kg/m³): {dens_Tess_Liq:.10f}\n")
            f.write(f"dens_Tess_Gas (kg/m³): {dens_Tess_Gas:.10f}\n")
            f.write(f"dens_Tess_Tra (kg/m³): {dens_Tess_Tra:.10f}\n")
            f.write(f"avg_vol_Liq (Å³):      {avg_vol_liq:.10f}\n")
            f.write(f"avg_vol_Gas (Å³):      {avg_vol_gas:.10f}\n")
            f.write(f"avg_vol_Tra (Å³):      {avg_vol_tra:.10f}\n")
            f.write(f"N_Liq_residues:        {len(liq_volumes_all)}\n")
            f.write(f"N_Gas_residues:        {len(gas_volumes_all)}\n")
            f.write(f"N_Tra_residues:        {len(tra_volumes_all)}\n\n")

        print("Volume classification by phase complete!")


def generate_phase_histograms_individual(job, gro_file, npz_liq, npz_gas, npz_tra, output_png, molar_mass=18.015):
    """
    Generate individual phase histograms from NPZ files.

    Args:
        job: signac job context
        gro_file: .gro file for reading box dimensions
        npz_liq, npz_gas, npz_tra: paths to phase-specific NPZ files
        output_png: output histogram filename
        molar_mass: molar mass in g/mol (default: water)
    """
    import matplotlib.pyplot as plt

    with(job):
        # Read box volume from .gro file
        box_volume = None
        box_Lx = box_Ly = box_Lz = None

        try:
            with open(gro_file, 'r') as f:
                lines = f.readlines()
                box_line = lines[-1].split()
                Lx, Ly, Lz = float(box_line[0]), float(box_line[1]), float(box_line[2])
                box_volume = Lx * Ly * Lz * 1000.0  # Convert to Å³
                box_Lx, box_Ly, box_Lz = Lx, Ly, Lz
                print(f"Box dimensions from {gro_file}:")
                print(f"  Lx = {box_Lx:.3f} nm, Ly = {box_Ly:.3f} nm, Lz = {box_Lz:.3f} nm")
                print(f"  Box volume = {box_volume:.2f} Å³")
        except Exception as e:
            print(f"Warning: Could not read box from {gro_file}: {e}")

        # Load all three NPZ files
        liq_volumes_flat = np.array([])
        gas_volumes_flat = np.array([])
        tra_volumes_flat = np.array([])

        # Load LIQ
        #if os.path.isfile(npz_liq):
        data_liq = np.load(npz_liq)
        liq_volumes = data_liq['volumes']
        liq_volumes_flat = liq_volumes.flatten() if liq_volumes.ndim > 1 else liq_volumes
        liq_volumes_flat = liq_volumes_flat[np.isfinite(liq_volumes_flat)]
        print(f"Loaded {len(liq_volumes_flat)} LIQ volumes")

        # Load GAS
        #if os.path.isfile(npz_gas):
        data_gas = np.load(npz_gas)
        gas_volumes = data_gas['volumes']
        gas_volumes_flat = gas_volumes.flatten() if gas_volumes.ndim > 1 else gas_volumes
        gas_volumes_flat = gas_volumes_flat[np.isfinite(gas_volumes_flat)]
        print(f"Loaded {len(gas_volumes_flat)} GAS volumes")

        # Load TRA
        #if os.path.isfile(npz_tra):
        data_tra = np.load(npz_tra)
        tra_volumes = data_tra['volumes']
        tra_volumes_flat = tra_volumes.flatten() if tra_volumes.ndim > 1 else tra_volumes
        tra_volumes_flat = tra_volumes_flat[np.isfinite(tra_volumes_flat)]
        print(f"Loaded {len(tra_volumes_flat)} TRA volumes")

        # Calculate densities
        avogadro = 6.022e23  # molecules/mol
        molecule_mass_kg = (molar_mass / avogadro) * 1e-3

        avg_volume_liq = np.mean(liq_volumes_flat) if len(liq_volumes_flat) > 0 else 0
        avg_volume_gas = np.mean(gas_volumes_flat) if len(gas_volumes_flat) > 0 else 0
        avg_volume_tra = np.mean(tra_volumes_flat) if len(tra_volumes_flat) > 0 else 0

        density_liq_kg_m3 = molecule_mass_kg / (avg_volume_liq * 1e-30) if avg_volume_liq > 0 else 0
        density_gas_kg_m3 = molecule_mass_kg / (avg_volume_gas * 1e-30) if avg_volume_gas > 0 else 0
        density_tra_kg_m3 = molecule_mass_kg / (avg_volume_tra * 1e-30) if avg_volume_tra > 0 else 0

        # Calculate total statistics
        total_N = len(liq_volumes_flat) + len(gas_volumes_flat) + len(tra_volumes_flat)
        total_volume_liq = np.sum(liq_volumes_flat) if len(liq_volumes_flat) > 0 else 0
        total_volume_gas = np.sum(gas_volumes_flat) if len(gas_volumes_flat) > 0 else 0
        total_volume_tra = np.sum(tra_volumes_flat) if len(tra_volumes_flat) > 0 else 0
        total_volume_sum = total_volume_liq + total_volume_gas + total_volume_tra

        # Create plots - Order: GAS → TRANSITION → LIQUID
        fig, axes = plt.subplots(1, 3, figsize=(18, 8), sharey=False)

        # GAS subplot (first)
        if len(gas_volumes_flat) > 0:
            axes[0].hist(gas_volumes_flat, bins=50, edgecolor='black', alpha=0.7, color='orange')
            axes[0].axvline(avg_volume_gas, color='darkorange', linestyle='--', linewidth=2,
                           label=f'Mean: {avg_volume_gas:.2f} Å³')
            axes[0].set_xlabel('Volume (Å³)', fontsize=11)
            axes[0].set_ylabel('Frequency', fontsize=11)
            axes[0].set_title(f'GAS Phase\\n{len(gas_volumes_flat)} points, ρ={density_gas_kg_m3:.1f} kg/m³',
                             fontsize=12, fontweight='bold')
            axes[0].legend(fontsize=9)
            axes[0].grid(True, alpha=0.3)
        else:
            axes[0].text(0.5, 0.5, 'No GAS data', ha='center', va='center',
                        transform=axes[0].transAxes, fontsize=14)
            axes[0].set_title('GAS Phase', fontsize=12, fontweight='bold')

        # TRANSITION subplot (middle)
        if len(tra_volumes_flat) > 0:
            axes[1].hist(tra_volumes_flat, bins=50, edgecolor='black', alpha=0.7, color='red')
            axes[1].axvline(avg_volume_tra, color='darkred', linestyle='--', linewidth=2,
                           label=f'Mean: {avg_volume_tra:.2f} Å³')
            axes[1].set_xlabel('Volume (Å³)', fontsize=11)
            axes[1].set_title(f'TRANSITION Phase\\n{len(tra_volumes_flat)} points, ρ={density_tra_kg_m3:.1f} kg/m³',
                             fontsize=12, fontweight='bold')
            axes[1].legend(fontsize=9)
            axes[1].grid(True, alpha=0.3)
        else:
            axes[1].text(0.5, 0.5, 'No TRANSITION data', ha='center', va='center',
                        transform=axes[1].transAxes, fontsize=14)
            axes[1].set_title('TRANSITION Phase', fontsize=12, fontweight='bold')

        # LIQUID subplot (last)
        if len(liq_volumes_flat) > 0:
            axes[2].hist(liq_volumes_flat, bins=50, edgecolor='black', alpha=0.7, color='blue')
            axes[2].axvline(avg_volume_liq, color='darkblue', linestyle='--', linewidth=2,
                           label=f'Mean: {avg_volume_liq:.2f} Å³')
            axes[2].set_xlabel('Volume (Å³)', fontsize=11)
            axes[2].set_title(f'LIQUID Phase\\n{len(liq_volumes_flat)} points, ρ={density_liq_kg_m3:.1f} kg/m³',
                             fontsize=12, fontweight='bold')
            axes[2].legend(fontsize=9)
            axes[2].grid(True, alpha=0.3)
        else:
            axes[2].text(0.5, 0.5, 'No LIQUID data', ha='center', va='center',
                        transform=axes[2].transAxes, fontsize=14)
            axes[2].set_title('LIQUID Phase', fontsize=12, fontweight='bold')

        # Overall title
        overall_title = f'Voronoi Volume Distribution by Phase\\n'
        overall_title += f'Total N={total_N} | '
        if total_volume_sum > 0:
            overall_title += f'Vol: LIQ={total_volume_liq/total_volume_sum*100:.1f}%, '
            overall_title += f'GAS={total_volume_gas/total_volume_sum*100:.1f}%, '
            overall_title += f'TRA={total_volume_tra/total_volume_sum*100:.1f}% | '
        overall_title += f'ΣVol={total_volume_sum:.1f} Å³'
        if box_volume is not None:
            overall_title += f' | BoxVol={box_volume:.1f} Å³'
        fig.suptitle(overall_title, fontsize=13, fontweight='bold', y=0.98)

        plt.tight_layout(rect=[0, 0, 1, 0.96])
        plt.savefig(output_png, dpi=300)
        print(f"Individual histograms saved to: {output_png}")
        plt.close()


###############################################################################
################## MULTI-METHOD DENSITY RESULTS AGGREGATION ###################
###############################################################################

def many_method_extractor(job, list_of_files2read, look_up_strs, destination_file):
    """
    Extract density data from multiple files and aggregate into a single results file.

    Args:
        job: signac job context
        list_of_files2read: list of filenames to read (relative to job directory)
        look_up_strs: list of strings to search for (extract float after each)
        destination_file: full path to output file (in project directory)

    Returns:
        dict: {filename: {lookup_str: float_value}}
    """
    import time
    import os
    import fcntl

    with(job):
        results = {}

        # Loop through files to read
        for filename in list_of_files2read:
            if not os.path.isfile(filename):
                print(f"Warning: {filename} not found in job {job.id}")
                continue

            results[filename] = {}

            with open(filename, 'r') as f:
                content = f.read()

            # Loop through lookup strings
            for lookup_str in look_up_strs:
                if lookup_str in content:
                    # Find the line containing the lookup string
                    for line in content.split('\n'):
                        if lookup_str in line:
                            # Extract the float value after the lookup string
                            try:
                                # Split by the lookup string and get the part after
                                after_str = line.split(lookup_str)[-1]
                                # Extract first number from the remaining string
                                import re
                                numbers = re.findall(r'[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?', after_str)
                                if numbers:
                                    results[filename][lookup_str] = float(numbers[0])
                                    break
                            except (ValueError, IndexError) as e:
                                print(f"Warning: Could not extract float for '{lookup_str}' in {filename}: {e}")

        # Build job info string
        job_info = f"{job.id:<42} {job.sp.r_cut:<8} {job.sp.cut_type:<8} {job.sp.replicas:<8} TEMP{job.sp.temperature:<8}"

        # Build data string from results
        data_values = []
        FIELD_WIDTH = 8
        SEP = " \t\t "

        #for filename in list_of_files2read:
        #    if filename in results:
        #        for lookup_str in look_up_strs:
        #            if lookup_str in results[filename]:
        #                data_values.append(f"{results[filename][lookup_str]:.6f}")
        #            else:
        #                data_values.append("N/A")

        for filename in list_of_files2read:
            if filename in results:
                for lookup_str in look_up_strs:
                    if lookup_str in results[filename]:
                        value = results[filename][lookup_str]
                        # right-align float with fixed width and 6 decimals
                        data_values.append(f"{value:>{FIELD_WIDTH}.6f}{SEP}")
                    else:
                        # right-align N/A in the same field width
                        data_values.append(f"{'N/A':>{FIELD_WIDTH}}{SEP}")

        data_str = " ".join(data_values)
        output_line = f"{job_info} {data_str}\n"

        # Try to write to destination file with file locking
        max_retries = 3600  # 1 hour of retries
        retry_delay = 1.0   # 1 second between retries

        for attempt in range(max_retries):
            try:
                # Open file in append mode
                with open(destination_file, 'a') as f:
                    # Try to acquire exclusive lock
                    fcntl.flock(f.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
                    try:
                        f.write(output_line)
                        print(f"Successfully wrote to {destination_file}")
                        return results
                    finally:
                        # Release lock
                        fcntl.flock(f.fileno(), fcntl.LOCK_UN)
            except (IOError, BlockingIOError):
                if attempt < max_retries - 1:
                    print(f"File {destination_file} is busy, retrying in {retry_delay}s... (attempt {attempt+1}/{max_retries})")
                    time.sleep(retry_delay)
                else:
                    print(f"ERROR: Could not write to {destination_file} after {max_retries} attempts")
                    raise

        return results
