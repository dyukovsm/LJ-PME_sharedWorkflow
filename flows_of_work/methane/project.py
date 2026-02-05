import math
import mbuild as mb
import numpy as np
import signac
from flow import FlowProject
from flow.environment import DefaultSlurmEnvironment
import os
import shutil
from files.python_files import names, job_tester, job_templates
import forcefield_utilities
import mbuild as mb
import gmso
from gmso.external.convert_mbuild import from_mbuild
from gmso.parameterization import apply
from gmso.formats.top import write_top
from gmso.formats.gro import write_gro
import matplotlib.pyplot as plt
import pandas as pd
import re
import subprocess
import io
import freud

PROJECT_FILES_DIR = os.path.abspath('files')
PROJECT_DIR = os.path.abspath('.')
MDP_DIR = 'mdp'; XYZ_DIR = 'coordinates'; XML_DIR = 'xml/trappe'

MIN_CORES = 1; BUILD_CORES = 2; MAX_CORES = 4; SIMULATION_GPU = 1; PADDING_CORES_4_CLASSIFYING_PHASES = 16
TINNY_MEM = 0.512; LOW_MEM = 1.024; HIGH_MEM = 2.048; FOUR_GIGS = 4.096; SIXTEEN_GIGS = 16.384; VORONOI_MEMORY = 16.0
SHORT_WAIT = 2.0; HALF_DAY = 8.0; DAY_WAIT = 24.0; MED_WAIT = 96.0; ONE_WORKWEEK = 111.0; TWO_WEEKS = 222.0

PRINT_MY_NODE = 'echo -e "Hello World\nHello World upcomming hostname"; hostname'

# for testing
#MID_EQ_STEPS        = int(1000)     # int(2000000) 
#LONG_EQ_STEPS       = int(1000)     # int(20000000)
#SLOW_OUTPUT         = int(100)      # int(10000)     
#SLOW_CALC           = int(100)      # int(100)      
#
#PRO_STEPS           = int(1000)     # int(10000000) 
#FAST_OUTPUT         = int(100)      # int(100)     
#FAST_CALC           = int(100)      # int(100)

# for running
MID_EQ_STEPS        = int(2000000) 
LONG_EQ_STEPS       = int(50000000)
SLOW_OUTPUT         = int(100000)
SLOW_CALC           = int(100)   

PRO_STEPS           = int(10000000) 
FAST_OUTPUT         = int(5000)     
FAST_CALC           = int(100)      

PLANNED_Z_ELONGATION = 74.2444
N_MOLECULES = int(40496)
INIT_CUBELENGTH = 14.2

GMX_PREFIX = names.GMX_PREFIX

PLANNED_Z_ELONGATION_PME = 14.2
N_MOLECULES_PME = int(3346)
INIT_CUBELENGTH_PME = 6.2 

current_directory = os.getcwd()
current_directory_name = os.path.basename(current_directory)
project = signac.get_project()

class Custom_environment(DefaultSlurmEnvironment):  

    hostname_pattern = r".*\.grid\.wayne\.edu"
    template = "gmx_grid_fall2025.sh"

    
###################################################################################################

@FlowProject.pre(job_tester.build_input_starter)
@FlowProject.post(job_tester.inits_written)
@FlowProject.post(job_tester.mdps_written)
@FlowProject.operation(directives={ "np": BUILD_CORES,  "ngpu": 0, "memory": HIGH_MEM, "walltime": SHORT_WAIT})
def build_input(job):
        
    with(job):
        
        if 'Cut-off' in job.sp.cut_type:
            working_n = N_MOLECULES
            working_l = INIT_CUBELENGTH
        elif 'PME' in job.sp.cut_type:
            working_n = N_MOLECULES_PME
            working_l = INIT_CUBELENGTH_PME
        
        methane = mb.load(f'{PROJECT_FILES_DIR}/{XYZ_DIR}/Methane.mol2')
        methane.name = 'MET'
        
        starting_box = mb.fill_box(compound = methane, n_compounds = working_n, box=[working_l, working_l, working_l])
        trappe_ff_xml = forcefield_utilities.GMSOFFs().load_xml(f'{PROJECT_FILES_DIR}/{XML_DIR}/gmx-units-trappe-mie.xml').to_gmso_ff()
        
        gmso_starting_box = from_mbuild(starting_box)
        
        #### force_field_dict = {
        ####     'MET' : trappe_ff_xml
        #### }

        for dummy in gmso_starting_box.sites:
            if 'MET' in dummy.molecule.name:
                dummy.label = '_CH4'
                #dummy.molecule.isrigid = True
                dummy.molecule.name = 'MET'
        
        print(f'before apply')
        
        apply(top=gmso_starting_box,
                forcefields= trappe_ff_xml, # force_field_dict,
                ) #identify_connections=True)
        
        #apply(top=gmso_starting_box,
        #        forcefields=force_field_dict,
        #        identify_connections=True)

        write_gro(gmso_starting_box, filename='init.gro')
        write_top(gmso_starting_box, filename='init.top')
    
    #eqNVT
    parameters = {
        'integrator' : 'md',
        'nsteps' : MID_EQ_STEPS,
        'output_control' : SLOW_OUTPUT,
        'nstcalcenergy' : SLOW_CALC,
        'nstlist' : 10,
        'rcoulomb' : job.sp.r_cut,
        'coulombtype' : 'Cut-off', #METHANE METHANE METHANE METHANE
        'coulomb_modifier' : 'None',
        'rcoulomb_switch' : 0.0,
        'vdw_type' : job.sp.cut_type,
        'vdw_modifier' : 'None',
        'rvdw' : job.sp.r_cut,
        'rvdw_switch' : 0.0,
        'DispCorr' : 'No',
        'tcouple' : 'nose-hoover',
        'ref_t' : job.sp.temperature 
        }
    
    job_templates.simple_mdp_writer(job,mdp_name=f'{names.NAME_EQ_NVT}.mdp',parameters=parameters,constraints=None,templates_dir=f'{PROJECT_FILES_DIR}/mdp/',template_name='NVT_template_generic.mdp')

    #eqSURFTEN
    parameters.update({
        'nsteps' : LONG_EQ_STEPS,
        'output_control' : SLOW_OUTPUT,
        'nstcalcenergy' : FAST_CALC,
        })
    
    job_templates.simple_mdp_writer(job,mdp_name=f'{names.NAME_EQ_SURFTEN}.mdp',parameters=parameters,constraints=None,templates_dir=f'{PROJECT_FILES_DIR}/mdp/',template_name='NVT_template_generic.mdp')

    #proSURFTEN
    parameters.update({
        'nsteps' : PRO_STEPS,
        'output_control' : FAST_OUTPUT,
        'nstcalcenergy' : FAST_CALC,
        })
    
    job_templates.simple_mdp_writer(job,mdp_name=f'{names.NAME_PRO_SURFTEN}.mdp',parameters=parameters,constraints=None,templates_dir=f'{PROJECT_FILES_DIR}/mdp/',template_name='NVT_template_generic.mdp')


###################################################################################################


@FlowProject.pre(job_tester.inits_written)
@FlowProject.pre(job_tester.mdps_written)
@FlowProject.post(job_tester.eq_nvt_post_em_done)
@FlowProject.operation(directives={ "np": MAX_CORES,  "ngpu": SIMULATION_GPU, "memory": HIGH_MEM, "walltime": MED_WAIT},with_job=True,cmd=True)
def EQ_NVT(job):
    
    build_mdp = str(GMX_PREFIX + ' grompp -f ' + f'{names.NAME_EQ_NVT}.mdp -c ' + f'init.gro -p ' + 'init.top -o ' + f'{names.NAME_EQ_NVT}.tpr -maxwarn 100')
    run_mdp = str(GMX_PREFIX + f' mdrun -nt ' + f'{MAX_CORES}' + ' -deffnm' + f' {names.NAME_EQ_NVT}')
    run_command = str(PRINT_MY_NODE + '; ' + 'sleep 6' + '; ' + build_mdp + '; ' + 'sleep 16' + '; ' + run_mdp)
    
    return run_command


###################################################################################################


@FlowProject.pre(job_tester.inits_written) 
@FlowProject.pre(job_tester.mdps_written) 
@FlowProject.pre(job_tester.eq_nvt_post_em_done) 
@FlowProject.post(job_tester.build_surfTen_nvt_done) 
@FlowProject.operation(directives={ "np": MIN_CORES,  "ngpu": 0, "memory": LOW_MEM, "walltime": SHORT_WAIT})
def ELONGATE_FOR_SURFTEN(job):
    with(job):   
        initialBox =  mb.load(f'{names.NAME_EQ_NVT}.gro')
        boxLength = initialBox.box.lengths
        
        outputFile = open(f'{names.NAME_ELONGATED}.gro','w')
        
        dummyFile = open(f'{names.NAME_EQ_NVT}.gro','r')
        lines = dummyFile.readlines(); dummyFile.close()
        
        if 'Cut-off' in job.sp.cut_type:
            zLength = PLANNED_Z_ELONGATION # INIT_CUBELENGTH
        elif 'PME' in job.sp.cut_type:
            zLength = PLANNED_Z_ELONGATION_PME # INIT_CUBELENGTH_PME
        
        for i in range(len(lines)):
            if i == len(lines)-1:
                #zLength=PLANNED_Z_ELONGATION #boxLength[2]*3 # otherwise get a tuple can't be changed error
                outputFile.write(f'   {boxLength[0]}   {boxLength[1]}   {zLength}\n')
            else:
                outputFile.write(f'{lines[i]}')
                
        outputFile.close()


###################################################################################################


@FlowProject.pre(job_tester.inits_written) 
@FlowProject.pre(job_tester.mdps_written) 
@FlowProject.pre(job_tester.eq_nvt_post_em_done)
@FlowProject.pre(job_tester.build_surfTen_nvt_done)
@FlowProject.post(job_tester.eq_nvt_surften_done)
@FlowProject.operation(directives={ "np": MAX_CORES,  "ngpu": SIMULATION_GPU, "memory": HIGH_MEM, "walltime": TWO_WEEKS},with_job=True,cmd=True)
def EQ_SURFTEN(job):
    
    build_mdp = str(GMX_PREFIX + ' grompp -f ' + f'{names.NAME_EQ_SURFTEN}.mdp -c ' + f'{names.NAME_ELONGATED}.gro -p ' + 'init.top -o ' + f'{names.NAME_EQ_SURFTEN}.tpr -maxwarn 100')    
    # build_mdp = str(GMX_PREFIX + ' grompp -f ' + f'{names.NAME_EQ_SURFTEN}.mdp -c ' + f'{names.NAME_ELONGATED}.gro -p ' + 'init.top -o ' + f'{names.NAME_EQ_SURFTEN}.tpr -maxwarn 100')
    run_mdp = str(GMX_PREFIX + f' mdrun -nt ' + f'{MAX_CORES}' + ' -deffnm' + f' {names.NAME_EQ_SURFTEN}')
    run_command = str(PRINT_MY_NODE + '; ' + 'sleep 6' + '; ' + build_mdp + '; ' + 'sleep 16' + '; ' + run_mdp)
    
    return run_command


###################################################################################################


@FlowProject.pre(job_tester.inits_written) 
@FlowProject.pre(job_tester.mdps_written) 
@FlowProject.pre(job_tester.eq_nvt_post_em_done)
@FlowProject.pre(job_tester.build_surfTen_nvt_done)
@FlowProject.pre(job_tester.eq_nvt_surften_done)
@FlowProject.post(job_tester.pro_nvt_surften_done)
@FlowProject.operation(directives={ "np": MAX_CORES,  "ngpu": SIMULATION_GPU, "memory": HIGH_MEM, "walltime": ONE_WORKWEEK},with_job=True,cmd=True)
def PRO_SURFTEN(job):
    
    build_mdp = str(GMX_PREFIX + ' grompp -f ' + f'{names.NAME_PRO_SURFTEN}.mdp -c ' + f'{names.NAME_EQ_SURFTEN}.gro -p ' + 'init.top -o ' + f'{names.NAME_PRO_SURFTEN}.tpr -maxwarn 100')
    run_mdp = str(GMX_PREFIX + f' mdrun -nt ' + f'{MAX_CORES}' + ' -deffnm' + f' {names.NAME_PRO_SURFTEN}')
    run_command = str(PRINT_MY_NODE + '; ' + 'sleep 6' + '; ' + build_mdp + '; ' + 'sleep 16' + '; ' + run_mdp)
    
    return run_command


###################################################################################################

@FlowProject.pre(job_tester.pro_nvt_surften_done)
@FlowProject.post(job_tester.data_collected)
@FlowProject.operation(directives={ "np": BUILD_CORES,  "ngpu": SIMULATION_GPU, "memory": LOW_MEM, "walltime": SHORT_WAIT})
def GRAPH_AND_COLLECT_PROPERTIES(job):
    with(job):
        
        output_file = names.NAME_PRO_SURFTEN
        #output_file = f'{names.NAME_TEMP_RAMP_STOP}' # names.PRO_SURFTEN_CHUNK_TO_STARTING_GRO_FILE[last_completed_chunk+1]

        # try to match both this and the strings below with gromacs promp
        properties_of_interest = ["Potential", "LJ-(SR)", "Coulomb-(SR)", "Coul.-recip.", "Total-Energy", "Vir-ZZ", "Pres-ZZ", "#Surf*SurfTen"]
        
        # properties_of_interest_to_search_string_dict[properties_of_interest] has the exact strings to search for
        properties_of_interest_to_search_string_dict = {
            # ex:
            # Disper.-corr.           : ['Disper. corr.','(kJ/mol)'] 
            # Pres-ZZ                 : ['Pres-ZZ','(bar)']
            # Total-Energy            : ['Total Energy','(kJ/mol)']
            properties_of_interest[0] : ['Potential','(kJ/mol)'],
            properties_of_interest[1] : ['LJ (SR)','(kJ/mol)'],
            properties_of_interest[2] : ['Coulomb (SR)','(kJ/mol'],
            properties_of_interest[3] : ['Coul. recip.','(kJ/mol'],
            properties_of_interest[4] : ['Total Energy','(kJ/mol)'],
            properties_of_interest[5] : ['Vir-ZZ','(kJ/mol)'], # gotto love how gmx swaps them ;p
            properties_of_interest[6] : ['Pres-ZZ','(bar)'],
            properties_of_interest[7] : ['#Surf*SurfTen','(bar nm)']
        }
        
        properties_of_interest_storage_dict = {
            properties_of_interest[0] : 0.0 ,
            properties_of_interest[1] : 0.0 ,
            properties_of_interest[2] : 0.0 ,
            properties_of_interest[3] : 0.0 ,
            properties_of_interest[4] : 0.0 ,
            properties_of_interest[5] : 0.0 ,
            properties_of_interest[6] : 0.0,
            properties_of_interest[7] : 0.0
        }

        gromacs_input = b'1\n0\n'
        result = subprocess.run(
        [f"{names.GMX_PREFIX}", "energy", "-f", f"{output_file}.edr", "-o", "dummy_data.xvg"],
        input=gromacs_input,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT 
        )
        
        with open("gmx_energy_index_reader.txt", "wb") as f:
            f.write(result.stdout)
            
        with open("gmx_energy_index_reader.txt", "r") as f:
            text = f.read()
            
        pattern = r'(\d+)\s+(\S+)'

        matches = re.findall(pattern, text, re.MULTILINE)
        index_map = {}
        for index, name in matches:
            clean_name = name.strip()
            index_map[clean_name] = int(index)

        results = {}
        for prop in properties_of_interest:
            if prop in index_map:
                results[prop] = index_map[prop]
                print(f'index_map[prop] {index_map[prop]}')
            elif "Total" in index_map or "Total" in prop:
                print(f'index_map: {index_map} did not match with prop {prop}')
                
        
        newline_string = "\n".join(str(results[prop]) for prop in properties_of_interest if prop in results)
        
        # run gmx density just for lolz
        p = subprocess.Popen([f'{names.GMX_PREFIX}', 'density', '-f', f'{output_file}.trr', '-s', f'{output_file}.tpr', '-o', 'dens_profile.xvg','-sl','88'], stdin=subprocess.PIPE,stdout=subprocess.PIPE)
        out,err = p.communicate(f'0\n0\n'.encode('utf-8'))#(b'13\n14\n15\n16\n17\n18\n19\n20\n0\n')
        capture = out.decode()

        ###### --------- Ensure the indexes correspond to properteis we need  --------- ######
        
        p = subprocess.Popen([f'{names.GMX_PREFIX}', '-quiet', 'energy', '-f', f'{output_file}.edr', '-o', f'{names.GENERAL_LOCAL_DATA}_{output_file}.xvg'], stdin=subprocess.PIPE,stdout=subprocess.PIPE)
        out,err = p.communicate(f'{newline_string}'.encode('utf-8'))#(b'13\n14\n15\n16\n17\n18\n19\n20\n0\n')
        capture = out.decode()
        
        Dummy_GMX_output = open(f'{names.GENERAL_LOCAL_DATA}_{output_file}.txt','w')
        Dummy_GMX_output.write(capture)
        Dummy_GMX_output.close()
        
        Dummy_GMX_output = open(f'{names.GENERAL_LOCAL_DATA}_{output_file}.txt','r')
        aggregate_surTenFile = open(f"../../{names.GENERAL_GLOBAL_DATA}.txt",'a')
        
                
        for a_single_line in Dummy_GMX_output:
            for property_str in properties_of_interest:
                search_property_str_dict = properties_of_interest_to_search_string_dict[property_str]
                search_str_start = search_property_str_dict[0]
                search_str_end = search_property_str_dict[1]
                
                if (search_str_start in a_single_line) and (search_str_end in a_single_line):
                    numpyCatcher=np.fromstring(a_single_line.strip(f'{search_str_start}{search_str_end}'),dtype=float,sep=' ')[0]
                    properties_of_interest_storage_dict[property_str] = numpyCatcher
                                    
            # read Dummy_GMX_output write to aggregate_surTenFile
            
        aggregate_surTenFile.write(f"{job.id:<42} {job.sp.r_cut:<8} {job.sp.cut_type:<8} {job.sp.replicas:<8} TEMP{job.sp.temperature:<8}"
                                   f" {properties_of_interest_storage_dict[properties_of_interest[0]]:<42} " 
                                   f" {properties_of_interest_storage_dict[properties_of_interest[1]]:<42} " 
                                   f" {properties_of_interest_storage_dict[properties_of_interest[2]]:<42} " 
                                   f" {properties_of_interest_storage_dict[properties_of_interest[3]]:<42} "
                                   f" {properties_of_interest_storage_dict[properties_of_interest[4]]:<42} "
                                   f" {properties_of_interest_storage_dict[properties_of_interest[5]]:<42} "
                                   f" {properties_of_interest_storage_dict[properties_of_interest[6]]:<42} "
                                   f" {properties_of_interest_storage_dict[properties_of_interest[7]]:<42} "
                                   "\n")
        
        
        
        ### graph the data of simulation to see if it's converged

        Dummy_GMX_output.close(); aggregate_surTenFile.close()
        
        xvg_png_datasource = open(f'{names.GENERAL_LOCAL_DATA}_{output_file}.xvg','r')
        #lines = xvg_png_datasource.strip().split('\n')
        lines = xvg_png_datasource.readlines()
        
        header_lines = []
        data_lines = []
        in_data_section = False
        
        
        for line in lines:
            if line.startswith('@') or line.startswith('#'):
                header_lines.append(line)
                if "@TYPE xy" in line:
                    in_data_section = True
            else : # in_data_section # #and line.strip(): #and not line.startswith('#'):
                data_lines.append(line.strip())

        # Extract column names from header
        column_names = {}
        xaxis_label = "Time (ps)" # Default, will be updated if found
        yaxis_label = ""
        title = ""

        for line in header_lines:
            if line.startswith('@ s'):
                match = re.search(r'@ s(\d+) legend "(.+)"', line)
                if match:
                    col_index = int(match.group(1))
                    legend_name = match.group(2)
                    column_names[col_index] = legend_name
            elif line.startswith('@ xaxis'):
                match = re.search(r'@ xaxis\s+label "(.+)"', line)
                if match:
                    xaxis_label = match.group(1)
            elif line.startswith('@ yaxis'):
                match = re.search(r'@ yaxis\s+label "(.+)"', line)
                if match:
                    yaxis_label = match.group(1)
            elif line.startswith('@ title'):
                match = re.search(r'@\s+title "(.+)"', line)
                if match:
                    title = match.group(1)

        # Create a list of column names in order
        # The first column is always the x-axis label.
        # Then append the other column names based on their index.
        ordered_column_names = [xaxis_label]
        max_col_index = max(column_names.keys()) if column_names else -1
        for i in range(max_col_index + 1):
            if i in column_names:
                ordered_column_names.append(column_names[i])

        # Read data into a pandas DataFrame
        # Using io.StringIO to treat the list of data lines as a file
        #print(f'data_lines: {data_lines}')
        df = pd.read_csv(io.StringIO("\n".join(data_lines)), sep=r'\s+', header=None)

        # Assign column names to the DataFrame
        df.columns = ordered_column_names[:len(df.columns)]

        # Plotting the data
        num_cols = len(df.columns) - 1  # Exclude the 'Time (ps)' column
        fig, axes = plt.subplots(num_cols, 1, figsize=(10, 5 * num_cols), sharex=True)

        if num_cols == 1:
            axes = [axes] # Ensure axes is iterable even for a single subplot
        
        for i, col_name in enumerate(df.columns[1:]): # Iterate over data columns, skipping time
            axes[i].plot(df[xaxis_label], df[col_name])
            axes[i].set_ylabel(f'{col_name} {yaxis_label}') # Add yaxis_label to each subplot's y-label
            axes[i].grid(True)
            #axes[i].set_title(f'{title}: {col_name}') # <- from chat
            key_to_mean_data = ''
            for key, value_list in properties_of_interest_to_search_string_dict.items(): 
                if col_name in value_list[0]:
                    key_to_mean_data = key
            if key_to_mean_data != '':
                axes[i].set_title(f'{col_name}; mean {properties_of_interest_storage_dict[key_to_mean_data]}') # <- hooman edited.

        axes[-1].set_xlabel(xaxis_label) # Set x-label only on the last subplot

        plt.tight_layout()
        plt.savefig(f'{names.GENERAL_LOCAL_DATA}_{output_file}.png')
        plt.close()


###################################################################################################
###################################################################################################
#### DENSITY PROFILE ALIGNMENT AND PHASE CLASSIFICATION
###################################################################################################
###################################################################################################

@FlowProject.pre(job_tester.pro_nvt_surften_done)
@FlowProject.post(job_tester.slab_aligned)
@FlowProject.operation(directives={ "np": BUILD_CORES,  "ngpu": 0, "memory": HIGH_MEM, "walltime": SHORT_WAIT})
def ALIGN_SLAB(job):
    """
    Align liquid slab to box center.

    This operation:
    1. Computes density profile and finds liquid-gas interfaces
    2. Shifts slab so liquid phase is centered at Z_LEN/2 (handles PBC wrapping)
    3. Generates aligned trajectory files for phase classification and Voronoi analysis
    """
    with(job):
        output_file = names.NAME_PRO_SURFTEN

        # Input files
        tpr_file = f'{output_file}.tpr'
        trr_file = f'{output_file}.trr'
        gro_file = f'{output_file}.gro'

        # Output files (aligned)
        output_trr = names.NAME_ALIGNED_TRR
        output_gro = names.NAME_ALIGNED_GRO

        # Run alignment
        alignment_data = job_templates.align_slab_to_center(
            job,
            tpr_file=tpr_file,
            trr_file=trr_file,
            gro_file=gro_file,
            output_trr=output_trr,
            output_gro=output_gro,
            bin_size=names.BIN_SIZE
        )

        print(f"Slab alignment complete. Aligned files created:")
        print(f"  {output_trr}")
        print(f"  {output_gro}")
        print(f"  {names.NAME_SHIFT_COMPARISON_PNG}")


@FlowProject.pre(job_tester.slab_aligned)
@FlowProject.post(job_tester.phases_classified)
@FlowProject.operation(directives={ "np": BUILD_CORES,  "ngpu": 0, "memory": HIGH_MEM, "walltime": SHORT_WAIT})
def CLASSIFY_PHASES(job):
    """
    Robustly classify liquid, gas, and transition phases from aligned slab.

    This operation:
    1. Computes density profile from aligned trajectory
    2. Fits tanh function to extract interface parameters
    3. Uses derivative-based noise analysis to robustly classify phases
    4. Saves comprehensive phase data with liquid/gas/transition densities
    5. Saves z-coordinates for each phase (for Voronoi filtering)
    """
    with(job):
        output_file = names.NAME_PRO_SURFTEN

        # Input files (aligned)
        tpr_file = f'{output_file}.tpr'
        aligned_trr = names.NAME_ALIGNED_TRR
        aligned_gro = names.NAME_ALIGNED_GRO

        # Run phase classification
        phase_data = job_templates.classify_phases_from_aligned_slab(
            job,
            tpr_file=tpr_file,
            aligned_trr=aligned_trr,
            aligned_gro=aligned_gro,
            bin_size=names.BIN_SIZE,
            d_multiplier=names.D_TANH_THICKNESS_MULTIPLIER,
            noise_multiplier=names.DERIVATIVE_NOISE_THRESHOLD_MULTIPLIER
        )

        file1 = names.NAME_PHASE_DATA
        file2 = names.NAME_PHASE_Z_COORDS
        file3 = names.NAME_PHASE_CLASSIFICATION_PNG

        if job.isfile(file1) and job.isfile(file2) and job.isfile(file3):
            print(f"Phase classification complete. Output files created:")
            print(f"  {file1}")
            print(f"  {file2}")
            print(f"  {file3}")
        else :
            for i in range(6):
                print(f'totally not there -->>{file1} {file2} {file3}')


###################################################################################################
###################################################################################################
#### MDANALYSIS VORONOI WORKFLOW
###################################################################################################
###################################################################################################

@FlowProject.pre(job_tester.phases_classified)
@FlowProject.post(job_tester.mda_voronoi_tessellated)
@FlowProject.operation(directives={ "np": BUILD_CORES,  "ngpu": SIMULATION_GPU, "memory": VORONOI_MEMORY, "walltime": ONE_WORKWEEK})
def TESSELLATE_VORONOI(job):
    """
    Compute Voronoi tessellation on the FULL system (all residues).

    CRITICAL: Voronoi must be computed on ALL residues to get correct volumes.
    Phase classification happens in a separate operation.
    """
    with(job):
        # Use aligned files from ALIGN_SLAB operation
        gro_file = names.NAME_ALIGNED_GRO
        trr_file = names.NAME_ALIGNED_TRR

        # Compute Voronoi on full system only
        job_templates.compute_voronoi_volumes_with_mda(
            job,
            gro_file=gro_file,
            trr_file=trr_file,
            output_npz=names.NAME_RESIDUE_DATA_ALL,
            max_frames=names.MAX_FRAMES_TO_READ
        )

        print(f"Voronoi tessellation complete: {names.NAME_RESIDUE_DATA_ALL}")


# this is not optimized (eats 16Gb RAM),  and not needed for tanh.
@FlowProject.pre(job_tester.phases_classified)
@FlowProject.pre(job_tester.mda_voronoi_tessellated)
@FlowProject.post(job_tester.mda_voronoi_classified)
@FlowProject.operation(directives={ "np": PADDING_CORES_4_CLASSIFYING_PHASES,  "ngpu": 0, "memory": VORONOI_MEMORY, "walltime": SHORT_WAIT})
def CLASSIFY_VORONOI_BY_PHASE(job):
    """
    Classify pre-computed Voronoi volumes into liquid, gas, and transition phases.

    Uses z-coordinate binning to assign each residue's volume to the appropriate phase
    based on the phase classification from the density profile analysis.
    """
    with(job):
        # Classify volumes by phase using z-binning
        job_templates.classify_voronoi_by_phase(
            job,
            input_npz_all=names.NAME_RESIDUE_DATA_ALL,
            output_npz_liq=names.NAME_RESIDUE_DATA_LIQ,
            output_npz_gas=names.NAME_RESIDUE_DATA_GAS,
            output_npz_tra=names.NAME_RESIDUE_DATA_TRA,
            output_z_hist_png=names.NAME_Z_POSITION_HISTOGRAM
        )

        print(f"Volume classification complete:")
        print(f"  LIQ: {names.NAME_RESIDUE_DATA_LIQ}")
        print(f"  GAS: {names.NAME_RESIDUE_DATA_GAS}")
        print(f"  TRA: {names.NAME_RESIDUE_DATA_TRA}")
        print(f"  Z-histogram: {names.NAME_Z_POSITION_HISTOGRAM}")

# this is not needed for tanh.
@FlowProject.pre(job_tester.phases_classified)
@FlowProject.pre(job_tester.mda_voronoi_classified)
@FlowProject.post(job_tester.mda_individual_histograms_complete)
@FlowProject.operation(directives={ "np": MAX_CORES,  "ngpu": 0, "memory": VORONOI_MEMORY, "walltime": SHORT_WAIT})
def GENERATE_PHASE_HISTOGRAMS(job):
    """
    Generate individual phase histograms from NPZ files.
    """
    with(job):
        gro_file = names.NAME_ALIGNED_GRO

        job_templates.generate_phase_histograms_individual(
            job,
            gro_file=gro_file,
            npz_liq=names.NAME_RESIDUE_DATA_LIQ,
            npz_gas=names.NAME_RESIDUE_DATA_GAS,
            npz_tra=names.NAME_RESIDUE_DATA_TRA,
            output_png=names.NAME_PHASE_HISTOGRAM_INDIVIDUAL,
            molar_mass=names.CURRENT_MOLAR_MASS
        )


###################################################################################################
###################################################################################################
#### MULTI-METHOD DENSITY RESULTS AGGREGATION
###################################################################################################
###################################################################################################

@FlowProject.pre(job_tester.mda_individual_histograms_complete)
@FlowProject.post(job_tester.density_results_exported)
@FlowProject.operation(directives={ "np": BUILD_CORES,  "ngpu": 0, "memory": LOW_MEM, "walltime": SHORT_WAIT})
def EXPORT_DENSITY_RESULTS(job):
    """
    Export density results from multiple methods to a single aggregation file.

    Extracts data from slab_phase_classification_data.txt and writes to density_results.txt
    in the project directory with job info and all density values.
    """
    with(job):
        # Define files to read and lookup strings
        files_to_read = [names.NAME_PHASE_DATA]

        # Lookup strings for extracting values from slab_phase_classification_data.txt
        # Order matches the header in density_results.txt
        lookup_strings = [
            "minDens_raw (kg/m³):",       # minDens_raw (for minMaxDensGas)
            "maxDens_raw (kg/m³):",       # maxDens_raw (for minMaxDensLiq)
            "rho_gas_fit (kg/m³):",       # tanhGas
            "rho_liq_fit (kg/m³):",       # tanhLiq
            "d_fit (nm):",                # tanhd_fit
            "Gas density (robust mean):", # dens_stableDeriv_Tess_Gas
            "Transition density (robust mean):",  # dens_stableDeriv_Tess_Tra
            "Liquid density (robust mean):"      # dens_stableDeriv_Tess_Liq
        ]

        # Destination file in project directory
        destination_file = os.path.join(PROJECT_DIR, names.MANY_METHODS_JOBS_DATA)

        # Extract and write results
        job_templates.many_method_extractor(
            job,
            list_of_files2read=files_to_read,
            look_up_strs=lookup_strings,
            destination_file=destination_file
        )


if __name__ == '__main__':
    FlowProject().main()
    
    

