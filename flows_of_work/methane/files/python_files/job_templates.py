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
from files.python_files import names
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import MDAnalysis as mda
import matplotlib
matplotlib.use("Agg")
from collections import defaultdict
import fcntl

#GMX_PREFIX = 'gmx' # for grid
GMX_PREFIX = '/usr/local/gromacs/bin/gmx' # for potoff cluster

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


###############################################################################
##################### HELPER FUNCTIONS #######################################
###############################################################################

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


###############################################################################
##################### SLAB ALIGNMENT #########################################
###############################################################################

def align_slab_to_center(job, tpr_file, trr_file, gro_file, output_trr, output_gro, bin_size=0.45):
    """
    Align liquid slab to box center using PBC-aware math.

    This function:
    1. Computes density profile along z-axis
    2. Finds interfaces via derivative analysis (max/min of dρ/dz)
    3. Uses modulo arithmetic to find slab center (handles wrapped slabs)
    4. Shifts slab so liquid phase is centered at Z_LEN/2

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
    with(job):
        # Get box dimensions
        dummy_mb = mb.load(gro_file)
        Z_LEN = dummy_mb.box.Lz
        bin_count = int(np.round(Z_LEN / bin_size, 0))
        
        # Step 1: Compute initial density profile
        density_xvg_pre = 'density_profile_pre_shift.xvg'
        cmd = f"echo 0 | {names.GMX_PREFIX} density -s {tpr_file} -f {trr_file} -o {density_xvg_pre} -d Z -sl {bin_count}"
        subprocess.run(cmd, shell=True, capture_output=True, text=True)
       
        z_current, density_current = parse_density_xvg(density_xvg_pre)
       
        # Step 2: Find interfaces via derivative
        # max derivative = rising edge (gas→liquid)
        # min derivative = falling edge (liquid→gas)
        ddens_dz = np.gradient(density_current, z_current)
        max_idx = np.argmax(ddens_dz)
        min_idx = np.argmin(ddens_dz)
       
        z_rising = z_current[max_idx]   # gas→liquid transition
        z_falling = z_current[min_idx]  # liquid→gas transition
       
        # Step 3: PBC-aware slab center calculation (no branching)
        # Liquid spans from z_rising → z_falling going forward (with PBC wrap)
        # The modulo handles both normal and wrapped cases automatically
        slab_width = (z_falling - z_rising) % Z_LEN
        slab_center = (z_rising + slab_width / 2.0) % Z_LEN
       
        # Calculate shift to put slab center at Z_LEN/2
        z_midpoint = Z_LEN / 2.0
        shift_vector = z_midpoint - slab_center
       
        # Store pre-shift data for plotting
        z_pre = z_current.copy()
        density_pre = density_current.copy()
       
        # Step 4: Shift slab to center
        cmd_shift_trr = f"echo 0 | {names.GMX_PREFIX} trjconv -f {trr_file} -s {tpr_file} -o {output_trr} -trans 0.0 0.0 {shift_vector} -pbc mol"
        subprocess.run(cmd_shift_trr, shell=True, capture_output=True, text=True)
       
        cmd_shift_gro = f"echo 0 | {names.GMX_PREFIX} trjconv -f {gro_file} -s {tpr_file} -o {output_gro} -trans 0.0 0.0 {shift_vector} -pbc mol"
        subprocess.run(cmd_shift_gro, shell=True, capture_output=True, text=True)
       
        # Step 5: Compute post-shift density profile for verification
        density_xvg_post = 'density_profile_post_shift.xvg'
        cmd_post = f"echo 0 | {names.GMX_PREFIX} density -s {tpr_file} -f {output_trr} -o {density_xvg_post} -d Z -sl {bin_count}"
        subprocess.run(cmd_post, shell=True, capture_output=True, text=True)
       
        z_post, density_post = parse_density_xvg(density_xvg_post)
       
        # Generate comparison plot
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
       
        ax1.plot(z_pre, density_pre, 'b-', lw=2)
        ax1.axvline(z_rising, color='red', ls='--', label=f'Rising (gas→liq): {z_rising:.2f} nm')
        ax1.axvline(z_falling, color='green', ls='--', label=f'Falling (liq→gas): {z_falling:.2f} nm')
        ax1.axvline(z_midpoint, color='black', ls='-', label=f'Z/2: {z_midpoint:.2f} nm')
        ax1.axvline(slab_center, color='purple', ls=':', label=f'Slab center: {slab_center:.2f} nm')
        ax1.set_xlabel('z (nm)')
        ax1.set_ylabel('ρ (kg/m³)')
        ax1.set_title(f'PRE-SHIFT (width={slab_width:.2f} nm)')
        ax1.legend(fontsize=8)
        ax1.grid(True, alpha=0.3)
       
        ax2.plot(z_post, density_post, 'b-', lw=2)
        ax2.axvline(z_midpoint, color='black', ls='-', label=f'Z/2: {z_midpoint:.2f} nm')
        ax2.set_xlabel('z (nm)')
        ax2.set_ylabel('ρ (kg/m³)')
        ax2.set_title(f'POST-SHIFT (shifted by {shift_vector:.2f} nm)')
        ax2.legend(fontsize=8)
        ax2.grid(True, alpha=0.3)
       
        plt.tight_layout()
        plt.savefig(names.NAME_SHIFT_COMPARISON_PNG, dpi=150)
        plt.close()
       
        print(f"Slab alignment complete:")
        print(f"  z_rising={z_rising:.3f}, z_falling={z_falling:.3f}")
        print(f"  slab_width={slab_width:.3f}, slab_center={slab_center:.3f}")
        print(f"  shift_vector={shift_vector:.3f} nm")
        print(f"  Output: {output_trr}, {output_gro}")

        ##
        ##
        ##
        ##

        # ==========================================================
        # POST-SHIFT ANALYSIS: COM DRIFT + TIME-RESOLVED DENSITIES
        # ==========================================================

        # ---- Load aligned trajectory ----
        u = mda.Universe(tpr_file, output_trr)
        com_list = []
        time_list = []
        for ts in u.trajectory:
            com_list.append(u.atoms.center_of_mass())  # mass-weighted
            time_list.append(ts.time)                  # ps

        com = np.array(com_list)
        times_ps = np.array(time_list)
        nframes_actual = len(times_ps)
        A_to_nm = 10.0
        com = com / A_to_nm

        # ---- Mean drift relative to first frame ----
        drift_x = com[:, 0] - com[0, 0]
        drift_y = com[:, 1] - com[0, 1]
        drift_z = com[:, 2] - com[0, 2]

        mean_drift_x = np.mean(drift_x)
        mean_drift_y = np.mean(drift_y)
        mean_drift_z = np.mean(drift_z)

        ### ---- Center of Mass over time ----
        ##masses = np.array([atom.element.mass for atom in traj.topology.atoms])
        ##masses = masses.reshape(1, -1, 1)  # reshape for broadcasting
  
        ##coords = traj.xyz  # (frames, atoms, 3)
        ##com = np.sum(coords * masses, axis=1) / np.sum(masses)
  
        ##mean_drift_x = np.mean(com[:, 0] - com[0, 0])
        ##mean_drift_y = np.mean(com[:, 1] - com[0, 1])
        ##mean_drift_z = np.mean(com[:, 2] - com[0, 2])

        # ---- Save COM subplot figure ----
        fig_com, axs = plt.subplots(3, 1, sharex=True, figsize=(6, 8))

        #axs[0].plot(times_ps, com[:, 0])
        axs[0].plot(times_ps, drift_x)
        axs[0].set_ylabel("drift x from t=0 (nm)")
        #axs[0].set_ylabel("COM X (nm)")

        #axs[1].plot(times_ps, com[:, 1])
        axs[1].plot(times_ps, drift_y)
        axs[1].set_ylabel("drift y from t=0 (nm)")

        #axs[2].plot(times_ps, com[:, 2])
        axs[2].plot(times_ps, drift_z)
        axs[2].set_ylabel("drift z from t=0 (nm)")
        axs[2].set_xlabel("Time (ps)")

        plt.tight_layout()
        plt.savefig(names.NAME_COM_DRIFT_PNG, dpi=150)
        plt.close(fig_com)

        # ---- Time interval splitting (10 segments) ----
        total_time_ps = times_ps[-1]
        interval_edges = np.linspace(0.0, total_time_ps, 11)

        density_profiles = []
        interval_fit_center = []
        interval_fit_thicc = []

        for i in range(10):
            t_begin = interval_edges[i]
            t_end = interval_edges[i + 1]

            xvg_name = f"density_interval_post_shift_{i}.xvg"

            cmd_density = (
                f"echo 0 | {names.GMX_PREFIX} density "
                f"-f {output_trr} -s {tpr_file} "
                f"-o {xvg_name} -d Z -sl {bin_count} "
                f"-b {t_begin} -e {t_end}"
            )
            subprocess.run(cmd_density, shell=True, capture_output=True, text=True)
            print(f'analysing profile drift for job {job.id}')

            z_tmp, rho_tmp = parse_density_xvg(xvg_name)
            popt, pcov = curve_fit(_dual_tanh_profile, z_tmp, rho_tmp, p0=[np.min(rho_tmp),np.max(rho_tmp),np.mean(z_tmp),10,1.0,1.0])
            dummy1, dummy2, center_z, thicc, dummy4, dummy5 = popt
            interval_fit_center.append(center_z)
            interval_fit_thicc.append(thicc)
            density_profiles.append((z_tmp, rho_tmp))

        # ---- Plot all density profiles (red → violet) ----
        fig_den = plt.figure(figsize=(6, 5))
        cmap = plt.get_cmap("rainbow")
        colors = cmap(np.linspace(0, 1, 10))

        for i, (z_tmp, rho_tmp) in enumerate(density_profiles):
            plt.plot(
                z_tmp,
                rho_tmp,
                color=colors[i],
                label=f"{interval_edges[i]:.0f}-{interval_edges[i+1]:.0f} ps",
            )

        plt.xlabel("z (nm)")
        plt.ylabel("Density (kg/m³)")
        plt.title("Post-Shift Density Profiles (10 Time Segments)")
        plt.legend(fontsize=7)
        plt.tight_layout()
        plt.savefig(names.NAME_PROFILE_DRIFT_PNG, dpi=150)
        plt.close(fig_den)

        print(f"  Mean COM drift (nm):")
        print(f"    X: {mean_drift_x:.6f}")
        print(f"    Y: {mean_drift_y:.6f}")
        print(f"    Z: {mean_drift_z:.6f}")

        output_line = (f"{job.id}\t{job.sp.r_cut}\t{job.sp.cut_type}\t{job.sp.temperature} "
                       f"\t\tmeanCOM={np.mean(com[:, 2]):.6f} driftCOM={drift_z:.6f}"
                       f"\t\t10IntervalTanhSTDcenter={np.std(interval_fit_center):.6f}"
                       f"\t\t10IntervalTanhMEANthicc={np.mean(interval_fit_thicc):.6f} \t10IntervalTanhSTDthicc={np.std(interval_fit_thicc):.6f} \n")

        output_txt_project_dir = f'../../{names.NAME_ALIGNMENT_ANALYSIS}'
        max_retries = 32
        for attempt in range(max_retries):
            try:
                with open(output_txt_project_dir, 'a') as f:
                    fcntl.flock(f.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
                    try:
                        f.write(output_line)
                        print(f"Results written to {output_txt_project_dir}")
                    finally:
                        fcntl.flock(f.fileno(), fcntl.LOCK_UN)
                break
            except (IOError, BlockingIOError):
                if attempt < max_retries - 1:
                    time.sleep(1)
                else:
                    print(f"ERROR: Could not write to {output_txt_project_dir}")
        #log_file_of_alignment.write(f'{job.id}\t{job.sp.r_cut}\t{job.sp.cut_type}\t{job.sp.temperature}'
        #                            f'\t\tmeanCOM={np.mean(com[:, 2]):.6f} driftCOM={mean_drift_z:.6f} '
        #                            f'\t\t10IntervalTanhSTDcenter={np.std(interval_fit_center):.6f} '
        #                            f'\t\t10IntervalTanhMEANthicc={np.mean(interval_fit_thicc):.6f} \t10IntervalTanhSTDthicc={np.std(interval_fit_thicc):.6f}')
        #log_file_of_alignment.close()

        ##
        ##
        ##
        ##

        return {
            #'shift_vector': shift_vector,
            #'z_len': Z_LEN,
            #'slab_width': slab_width,
            #'slab_center': slab_center,
            #'z_rising': z_rising,
            #'z_falling': z_falling,
            'mean_com_drift_x': mean_drift_x,
            'mean_com_drift_y': mean_drift_y,
            'mean_com_drift_z': mean_drift_z
        }


###############################################################################
##################### DUAL-TANH DENSITY PROFILE FITTING #######################
###############################################################################

def _dual_tanh_profile(z, rho_G, rho_L, z_center, w, delta_L, delta_R):
    """
    Dual-tanh density profile for liquid slab.

    Parameters:
        z: z-coordinate array (nm)
        rho_G: gas phase density (kg/m³)
        rho_L: liquid phase density (kg/m³)
        z_center: center of slab (nm)
        w: slab width (nm)
        delta_L: left interface thickness (nm)
        delta_R: right interface thickness (nm)

    Returns:
        density profile (kg/m³)
    """
    return (
        rho_G
        + (rho_L-rho_G) / 2.0
        * (
            np.tanh((z_center + w/2.0 - z) / delta_R)
            + np.tanh((z - z_center + w/2.0) / delta_L)
        )
    )


def dual_tanh_fit(job, tpr_file, trr_file, gro_file, output_png, output_txt_project_dir, bin_size=0.45):
    """
    Fit dual-tanh to ALIGNED density profile (slab centered at Z/2).

    Matches practice_dualfit.py exactly - simple curve_fit, no bounds.

    Args:
        job: signac job context
        tpr_file: input TPR file
        trr_file: input TRR trajectory (should be ALIGNED)
        gro_file: input GRO file (should be ALIGNED)
        output_png: output plot filename (in job dir)
        output_txt_project_dir: full path to output txt file (in project dir)
        bin_size: bin size for density profile (nm)

    Returns:
        dict with fit results: rho_G, rho_L, z_center, w, dL, dR
    """
    import fcntl

    with(job):
        # Get box dimensions
        dummy_mb = mb.load(gro_file)
        Z_LEN = dummy_mb.box.Lz
        bin_count = int(np.round(Z_LEN / bin_size, 0))

        # Run gmx density on aligned trajectory
        density_xvg = 'density_profile_for_fit.xvg'
        cmd = f"echo 0 | {names.GMX_PREFIX} density -s {tpr_file} -f {trr_file} -o {density_xvg} -d Z -sl {bin_count}"
        subprocess.run(cmd, shell=True, capture_output=True, text=True)

        # Parse XVG file
        z_bins, density = [], []
        with open(density_xvg, 'r') as f:
            for line in f:
                if not (line.startswith('#') or line.startswith('@')):
                    parts = line.split()
                    if len(parts) >= 2:
                        z_bins.append(float(parts[0]))
                        density.append(float(parts[1]))
        z = np.array(z_bins)
        rho = np.array(density)

        # Simple initial guesses (matching practice_dualfit.py)
        # Since slab is centered: gas at edges, liquid in center
        z_mid = Z_LEN / 2.0
        edge_fraction = 0.2  # use outer 20% for gas guess
        center_fraction = 0.3  # use central 30% for liquid guess

        gas_mask = (z < z.min() + edge_fraction * Z_LEN) | (z > z.max() - edge_fraction * Z_LEN)
        liq_mask = (z > z_mid - center_fraction * Z_LEN / 2) & (z < z_mid + center_fraction * Z_LEN / 2)

        p0 = [
            np.mean(rho[gas_mask]),   # rho_G - gas at edges
            np.mean(rho[liq_mask]),   # rho_L - liquid in center
            z_mid,                     # z_center - slab is centered
            10.0,                      # w - initial width guess
            1.0,                       # delta_L
            1.0,                       # delta_R
        ]

        # Simple fit - no bounds, just like practice_dualfit.py
        try:
            popt, pcov = curve_fit(_dual_tanh_profile, z, rho, p0=p0)
            rho_G, rho_L, z_center, w, dL, dR = popt
            fit_success = True

            # Calculate R² (coefficient of determination)
            rho_predicted = _dual_tanh_profile(z, *popt)
            ss_res = np.sum((rho - rho_predicted) ** 2)
            ss_tot = np.sum((rho - np.mean(rho)) ** 2)
            r_squared = 1.0 - (ss_res / ss_tot)
        except Exception as e:
            print(f"Fit failed: {e}")
            rho_G, rho_L, z_center, w, dL, dR = p0
            fit_success = False
            r_squared = 0.0

        # Generate plot (matching practice_dualfit.py style)
        z_fine = np.linspace(z.min(), z.max(), 1500)
        rho_fit = _dual_tanh_profile(z_fine, rho_G, rho_L, z_center, w, dL, dR)

        plt.figure(figsize=(10, 6))
        plt.scatter(z, rho, s=25, label="MD data")
        plt.plot(z_fine, rho_fit, 'k-', lw=2, label=(
            f"Fit\n"
            f"ρ_G={rho_G:.3f}, ρ_L={rho_L:.3f}\n"
            f"z_c={z_center:.3f}, w={w:.3f}\n"
            f"δ_L={dL:.3f}, δ_R={dR:.3f}"
        ))

        # Horizontal densities
        plt.axhline(rho_G, ls='--', color='C1', label=f"ρ_G={rho_G:.3f}")
        plt.axhline(rho_L, ls='--', color='C2', label=f"ρ_L={rho_L:.3f}")

        # Vertical geometry
        z_left = z_center - w / 2.0
        z_right = z_center + w / 2.0

        plt.axvspan(z_left - dR, z_left + dR, color="pink", alpha=0.3,
                   label=f"Left interface (δ_R={dR:.3f})")
        plt.axvspan(z_right - dL, z_right + dL, color="lightblue", alpha=0.3,
                   label=f"Right interface (δ_L={dL:.3f})")

        plt.axvline(z_center, ls=':', color='k', label=f"z_c={z_center:.3f}")
        plt.axvline(z_left, ls='--', color='gray')
        plt.axvline(z_right, ls='--', color='gray')

        plt.xlabel("z (nm)")
        plt.ylabel("ρ(z) (kg/m³)")
        plt.title(f"Dual-Tanh Fit - {job.id[:8]}")
        plt.legend(fontsize=9)
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig(output_png, dpi=300)
        plt.close()

        print(f"Dual-tanh fit complete:")
        print(f"  ρ_G = {rho_G:.4f} kg/m³")
        print(f"  ρ_L = {rho_L:.4f} kg/m³")
        print(f"  w = {w:.4f} nm")
        print(f"  δ_L = {dL:.4f} nm, δ_R = {dR:.4f} nm")
        print(f"  R² = {r_squared:.6f}")

        # Write to project-level txt file with locking
        output_line = f"{job.id}\t{job.sp.r_cut}\t{job.sp.cut_type}\t{job.sp.temperature}\t{rho_G:.6f}\t{rho_L:.6f}\t{w:.6f}\t{dL:.6f}\t{dR:.6f}\t{r_squared:.6f}\n"

        max_retries = 60
        for attempt in range(max_retries):
            try:
                with open(output_txt_project_dir, 'a') as f:
                    fcntl.flock(f.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
                    try:
                        f.write(output_line)
                        print(f"Results written to {output_txt_project_dir}")
                    finally:
                        fcntl.flock(f.fileno(), fcntl.LOCK_UN)
                break
            except (IOError, BlockingIOError):
                if attempt < max_retries - 1:
                    time.sleep(1)
                else:
                    print(f"ERROR: Could not write to {output_txt_project_dir}")

        return {
            'rho_G': rho_G,
            'rho_L': rho_L,
            'z_center': z_center,
            'w': w,
            'dL': dL,
            'dR': dR,
            'r_squared': r_squared,
            'fit_success': fit_success
        }


###############################################################################
##################### GENERALIZED GAUSSIAN DENSITY PROFILE FIT ################
###############################################################################

from scipy.special import gamma as gamma_func

def _gen_gaussian_profile(z, rho_G, rho_L, z_center, sigma, beta):
    """
    Generalized Gaussian (generalized normal) density profile for a liquid slab.

    The normalized generalized Gaussian is:
    f(x) = (β / (2σΓ(1/β))) * exp(-|(x - μ)/σ|^β)

    For density fitting with amplitude A (peak height above baseline):
    ρ(z) = ρ_G + A * (β / (2σΓ(1/β))) * exp(-|(z - z_center)/σ|^β)

    At z = z_center, ρ = ρ_L, so:
    A = (ρ_L - ρ_G) * (2σΓ(1/β)) / β

    Parameters:
        z: z-coordinate array (nm)
        rho_G: gas phase density (kg/m³)
        rho_L: liquid phase density (kg/m³)
        z_center: center of slab (nm)
        sigma: width parameter (nm)
        beta: shape parameter
              β = 2: standard Gaussian
              β < 2: heavier tails (more peaked)
              β > 2: lighter tails (flat-topped, more slab-like)

    Returns:
        density profile (kg/m³)
    """
    # Normalization factor for generalized Gaussian
    norm = beta / (2.0 * sigma * gamma_func(1.0 / beta))
    # Amplitude to achieve ρ_L at center
    A = (rho_L - rho_G) / norm
    return rho_G + A * norm * np.exp(-np.abs((z - z_center) / sigma)**beta)

import numpy as np

def _rational_profile(z, rho_G, rho_L, z_center, sigma, n):
    """
    Rational (Butterworth) density profile.
    n controls the sharpness of the shoulder. 
    """
    # Using 2*n and abs to ensure a symmetric, smooth, positive denominator
    denominator = 1 + np.abs((z - z_center) / sigma)**(2 * n)
    return rho_G + (rho_L - rho_G) / denominator

def gen_gaussian_fit(job, tpr_file, trr_file, gro_file, output_png, output_txt_project_dir, bin_size=0.45):
    """
    Fit generalized Gaussian profile to ALIGNED density data.

    ρ(z) = ρ_G + (ρ_L - ρ_G) * exp(-|(z - z_center) / σ|^β)

    The shape parameter β allows the fit to adapt:
    - β = 2: standard Gaussian
    - β > 2: flat-topped slab (more realistic for liquid)

    Args:
        job: signac job context
        tpr_file: input TPR file
        trr_file: input TRR trajectory (should be ALIGNED)
        gro_file: input GRO file (should be ALIGNED)
        output_png: output plot filename (in job dir)
        output_txt_project_dir: full path to output txt file (in project dir)
        bin_size: bin size for density profile (nm)

    Returns:
        dict with fit results: rho_G, rho_L, z_center, sigma, beta, r_squared
    """
    import fcntl

    with(job):
        # Get box dimensions
        dummy_mb = mb.load(gro_file)
        Z_LEN = dummy_mb.box.Lz
        bin_count = int(np.round(Z_LEN / bin_size, 0))

        # Run gmx density on aligned trajectory
        density_xvg = 'density_profile_for_gengauss_fit.xvg'
        cmd = f"echo 0 | {names.GMX_PREFIX} density -s {tpr_file} -f {trr_file} -o {density_xvg} -d Z -sl {bin_count}"
        subprocess.run(cmd, shell=True, capture_output=True, text=True)

        # Parse XVG file
        z, rho = parse_density_xvg(density_xvg)

        # Initial guesses
        z_mid = Z_LEN / 2.0
        edge_fraction = 0.2
        center_fraction = 0.3

        gas_mask = (z < z.min() + edge_fraction * Z_LEN) | (z > z.max() - edge_fraction * Z_LEN)
        liq_mask = (z > z_mid - center_fraction * Z_LEN / 2) & (z < z_mid + center_fraction * Z_LEN / 2)

        p0 = [
            np.mean(rho[gas_mask]),   # rho_G
            np.mean(rho[liq_mask]),   # rho_L
            z_mid,                    # z_center
            5.0,                      # sigma
            2.0,                      # beta (start with Gaussian)
        ]

        # Fit generalized Gaussian to density profile
        try:
            popt, pcov = curve_fit(_gen_gaussian_profile, z, rho, p0=p0)
            rho_G, rho_L, z_center, sigma, beta = popt
            fit_success = True

            # Calculate R²
            rho_predicted = _gen_gaussian_profile(z, *popt)
            ss_res = np.sum((rho - rho_predicted) ** 2)
            ss_tot = np.sum((rho - np.mean(rho)) ** 2)
            r_squared = 1.0 - (ss_res / ss_tot)
        except Exception as e:
            print(f"Generalized Gaussian fit failed: {e}")
            rho_G, rho_L, z_center, sigma, beta = p0
            fit_success = False
            r_squared = 0.0

        # Generate plot
        z_fine = np.linspace(z.min(), z.max(), 1500)
        rho_fit = _gen_gaussian_profile(z_fine, rho_G, rho_L, z_center, sigma, beta)

        plt.figure(figsize=(10, 6))
        plt.scatter(z, rho, s=25, label="MD data")
        plt.plot(z_fine, rho_fit, 'k-', lw=2, label=(
            f"Gen. Gaussian fit (R²={r_squared:.4f})\n"
            f"ρ_G={rho_G:.3f}, ρ_L={rho_L:.3f}\n"
            f"z_c={z_center:.3f}, σ={sigma:.3f}, β={beta:.3f}"
        ))

        # Horizontal densities
        plt.axhline(rho_G, ls='--', color='C1', label=f"ρ_G={rho_G:.3f}")
        plt.axhline(rho_L, ls='--', color='C2', label=f"ρ_L={rho_L:.3f}")

        # Vertical center and width region
        plt.axvline(z_center, ls=':', color='k', label=f"z_c={z_center:.3f}")
        plt.axvspan(z_center - sigma, z_center + sigma, color="pink", alpha=0.3,
                   label=f"σ={sigma:.3f}")

        plt.xlabel("z (nm)")
        plt.ylabel("ρ(z) (kg/m³)")
        plt.title(f"Generalized Gaussian Fit - {job.id[:8]}")
        plt.legend(fontsize=9)
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig(output_png, dpi=300)
        plt.close()

        print(f"Generalized Gaussian fit complete:")
        print(f"  ρ_G = {rho_G:.4f} kg/m³")
        print(f"  ρ_L = {rho_L:.4f} kg/m³")
        print(f"  z_center = {z_center:.4f} nm")
        print(f"  σ = {sigma:.4f} nm")
        print(f"  β = {beta:.4f}")
        print(f"  R² = {r_squared:.6f}")

        # Write to project-level txt file with locking
        output_line = f"{job.id}\t{job.sp.r_cut}\t{job.sp.cut_type}\t{job.sp.temperature}\t{rho_G:.6f}\t{rho_L:.6f}\t{z_center:.6f}\t{sigma:.6f}\t{beta:.6f}\t{r_squared:.6f}\n"

        max_retries = 6
        for attempt in range(max_retries):
            try:
                with open(output_txt_project_dir, 'a') as f:
                    fcntl.flock(f.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
                    try:
                        f.write(output_line)
                        print(f"Results written to {output_txt_project_dir}")
                    finally:
                        fcntl.flock(f.fileno(), fcntl.LOCK_UN)
                break
            except (IOError, BlockingIOError):
                if attempt < max_retries - 1:
                    time.sleep(1)
                else:
                    print(f"ERROR: Could not write to {output_txt_project_dir}")

        return {
            'rho_G': rho_G,
            'rho_L': rho_L,
            'z_center': z_center,
            'sigma': sigma,
            'beta': beta,
            'r_squared': r_squared,
            'fit_success': fit_success
        }
    
def gen_rational_fit(job, tpr_file, trr_file, gro_file, output_png, output_txt_project_dir, bin_size=0.45):
    """
    Fit Rational (Butterworth-like) profile to ALIGNED density data.

    ρ(z) = ρ_G + (ρ_L - ρ_G) / (1 + |(z - z_c) / σ|^(2n))

    The parameter n replaces beta:
    - Low n (1-2): Very smooth, heavy tails, rounded "cone"
    - High n (>4): Sharp shoulders, very flat top

    Returns:
        dict with fit results: rho_G, rho_L, z_center, sigma, n_order, r_squared
    """
    with(job):
        import fcntl
        # Get box dimensions (assuming mb is mbuild, as in your original)
        import mbuild as mb 
        dummy_mb = mb.load(gro_file)
        Z_LEN = dummy_mb.box.Lz
        bin_count = int(np.round(Z_LEN / bin_size, 0))

        # Run gmx density on aligned trajectory
        density_xvg = 'density_profile_for_rational_fit.xvg'
        # Assuming names.GMX_PREFIX is defined globally in your script
        cmd = f"echo 0 | gmx density -s {tpr_file} -f {trr_file} -o {density_xvg} -d Z -sl {bin_count}"
        subprocess.run(cmd, shell=True, capture_output=True, text=True)

        # Parse XVG file (assuming parse_density_xvg is defined in your script)
        z, rho = parse_density_xvg(density_xvg)

        # Initial guesses
        z_mid = Z_LEN / 2.0
        edge_fraction = 0.2
        center_fraction = 0.3

        gas_mask = (z < z.min() + edge_fraction * Z_LEN) | (z > z.max() - edge_fraction * Z_LEN)
        liq_mask = (z > z_mid - center_fraction * Z_LEN / 2) & (z < z_mid + center_fraction * Z_LEN / 2)

        # p0 = [rho_G, rho_L, z_center, sigma, n]
        p0 = [
            np.mean(rho[gas_mask]),   # rho_G
            np.mean(rho[liq_mask]),   # rho_L
            z_mid,                    # z_center
            Z_LEN / 4.0,              # sigma (starting width)
            2.0,                      # n (start with smooth n=2)
        ]

        # Fit Rational profile to density profile
        try:
            from scipy.optimize import curve_fit
            # Bounds: [rho_G, rho_L, z_c, sigma, n]
            # sigma and n must be > 0. n min=0.5 prevents singularity.
            lower_bounds = [-np.inf, -np.inf, 0, 0.01, 0.5]
            upper_bounds = [np.inf, np.inf, Z_LEN, Z_LEN, 20.0]
            
            popt, pcov = curve_fit(_rational_profile, z, rho, p0=p0, bounds=(lower_bounds, upper_bounds))
            rho_G, rho_L, z_center, sigma, n_order = popt
            fit_success = True

            # Calculate R²
            rho_predicted = _rational_profile(z, *popt)
            ss_res = np.sum((rho - rho_predicted) ** 2)
            ss_tot = np.sum((rho - np.mean(rho)) ** 2)
            r_squared = 1.0 - (ss_res / ss_tot)
        except Exception as e:
            print(f"Rational fit failed: {e}")
            rho_G, rho_L, z_center, sigma, n_order = p0
            fit_success = False
            r_squared = 0.0

        # Generate plot
        z_fine = np.linspace(z.min(), z.max(), 1500)
        rho_fit = _rational_profile(z_fine, rho_G, rho_L, z_center, sigma, n_order)

        plt.figure(figsize=(10, 6))
        plt.scatter(z, rho, s=25, label="MD data", alpha=0.5)
        plt.plot(z_fine, rho_fit, 'r-', lw=2, label=(
            f"Rational Fit (R²={r_squared:.4f})\n"
            f"ρ_G={rho_G:.3f}, ρ_L={rho_L:.3f}\n"
            f"z_c={z_center:.3f}, σ={sigma:.3f}, n={n_order:.3f}"
        ))

        # Horizontal densities
        plt.axhline(rho_G, ls='--', color='C1', label=f"ρ_G={rho_G:.3f}")
        plt.axhline(rho_L, ls='--', color='C2', label=f"ρ_L={rho_L:.3f}")

        # Vertical center and width region
        plt.axvline(z_center, ls=':', color='k', label=f"z_c={z_center:.3f}")
        plt.axvspan(z_center - sigma, z_center + sigma, color="cyan", alpha=0.1,
                   label=f"σ range")

        plt.xlabel("z (nm)")
        plt.ylabel("ρ(z) (kg/m³)")
        plt.title(f"Rational (Butterworth) Fit - {job.id[:8]}")
        plt.legend(fontsize=9)
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig(output_png, dpi=300)
        plt.close()

        print(f"Rational fit complete: ρ_G={rho_G:.3f}, ρ_L={rho_L:.3f}, R²={r_squared:.6f}")

        # Write to project-level txt file with locking
        output_line = f"{job.id}\t{job.sp.r_cut}\t{job.sp.cut_type}\t{job.sp.temperature}\t{rho_G:.6f}\t{rho_L:.6f}\t{z_center:.6f}\t{sigma:.6f}\t{n_order:.6f}\t{r_squared:.6f}\n"

        max_retries = 2
        for attempt in range(max_retries):
            try:
                with open(output_txt_project_dir, 'a') as f:
                    fcntl.flock(f.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
                    try:
                        f.write(output_line)
                    finally:
                        fcntl.flock(f.fileno(), fcntl.LOCK_UN)
                break
            except (IOError, BlockingIOError):
                time.sleep(0.5)

        return {
            'rho_G': rho_G,
            'rho_L': rho_L,
            'z_center': z_center,
            'sigma': sigma,
            'n_order': n_order,
            'r_squared': r_squared,
            'fit_success': fit_success
        }
    

