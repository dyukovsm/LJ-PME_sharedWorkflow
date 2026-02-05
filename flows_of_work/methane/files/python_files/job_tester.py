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


extension_list_of_common_files = [".gro", ".trr", ".log", ".edr", ".tpr"]   # former extension_list_list
#extension_list_inits = [".gro",".top"]
extension_list_inits = [".top"]


def return_file_with_extensions(file_names,extension_list):
    file_names_with_extensions = [name + ext for name in file_names for ext in extension_list]

    return file_names_with_extensions


def test_existence_simple(job,file_lsit):
    with(job):
        test_passed = False
        for i in file_lsit:
            if job.isfile(i):
                test_passed = True
            elif not(job.isfile(i)):
                test_passed = False
                break
        return test_passed
    
    
## def look_in_file(job,file_list,look_string,debug=False,check_for_not=False,check_for_not_str=['']):
##     with(job):
##         test_passed = False
##         #if debug:
##         #    missing_file = open(f"debug_look_IN_file_{file_list[0]}.txt",'w')
##         #    missing_file.write('test')
##         #    missing_file.close()
##         #    #close the file before
##         #    missing_file = open(f"debug_look_IN_file_{file_list[0]}.txt",'a')
##         #    #for i in file_list:
##         #    #    missing_file.write(f'{i}\n')
##         for i in file_list:
##             if job.isfile(i):
##                 file_with_lines = open(f'{i}','r')
##                 lines = file_with_lines.readlines()
##                 for j in lines:
##                     #f debug:
##                     #   if look_string not in j:
##                     #       missing_file.write(f'{look_string} not found in \t\t\t {j}\n')
##                     #   else:
##                     #       missing_file.write(f'{look_string} WAS FOUND in \t\t\t {j}\n')
##                     if check_for_not:
##                         for single_string in check_for_not_str:
##                             if single_string in j:
##                                 #print(f'single_string : {single_string} for job {job.id}')
##                                 return test_passed
##                         
##                     if look_string in j:
##                         test_passed = True
##                         break
##                 file_with_lines.close()
##             #elif debug:
##             #    missing_file.write(f'ERROR {i} not found.\n')
##     return test_passed

def look_in_file(job,file_list,look_string,debug=False,check_for_not=False,check_for_not_str=['']):
    with(job):
        test_passed = False

        for i in file_list:
            if job.isfile(i):
                file_with_lines = open(f'{i}','r')
                lines = file_with_lines.readlines()
                for j in lines:

                    if check_for_not:
                        for single_string in check_for_not_str:
                            if single_string in j:

                                return test_passed
                        
                    if look_string in j:
                        test_passed = True
                        break
                file_with_lines.close()

    return test_passed

def run_only_one(job):
    test_passed = False
    if job.sp.replica < 1:
        test_passed = True
    return test_passed


def important_jobs(job):
    test_passed = False
    if job.sp.r_cut < 6.0:
        #if 'PME' in job.sp.cut_type:
            test_passed = True
    return test_passed

############################__JOB_SPECIFIC_FUNCTIONS__############################


@FlowProject.label
def inits_written(job):
    with(job):
        check_these_files = return_file_with_extensions(file_names=['init'],extension_list = extension_list_inits)
        return test_existence_simple(job,check_these_files)
    

@FlowProject.label
def mdps_written(job):
    with(job):
        check_these_files = return_file_with_extensions(file_names=[names.NAME_EQ_NVT,
                                                                    names.NAME_EQ_SURFTEN,
                                                                    names.NAME_PRO_SURFTEN],extension_list = ['.mdp'])
        return test_existence_simple(job,check_these_files)
    
    
@FlowProject.label
def build_input_starter(job):
    with(job):
        starter_bool = not(inits_written(job)) and not(mdps_written(job))
        return starter_bool
    

# equilibrated slabs copies.


##################################################################################

# EQ_NVT -> ELONGATE -> EQ_SURFTEN -> PRO_SURFTEN -> ANALYSIS

##################################################################################

@FlowProject.label
def eq_nvt_post_em_done(job):
    with(job):
        bool_job = look_in_file(job,[f"{names.NAME_EQ_NVT}.log"],"Finished",check_for_not=True,check_for_not_str=['Received the TERM', 'Received the INT'])
        return bool_job


##################################################################################

@FlowProject.label
def build_surfTen_nvt_done(job):
    with(job):
        bool_job = False
        files_to_check = return_file_with_extensions(
            file_names=[f'{names.NAME_ELONGATED}'],
            extension_list=['.gro']
        )
        bool_job = test_existence_simple(job,files_to_check)
        return bool_job

##################################################################################


@FlowProject.label
def eq_nvt_surften_done(job):
    with(job):
        bool_job = look_in_file(job,[f"{names.NAME_EQ_SURFTEN}.log"],"Finished",check_for_not=True,check_for_not_str=['Received the TERM', 'Received the INT'])
        return bool_job


##################################################################################


@FlowProject.label
def pro_nvt_surften_done(job):
    with(job):
        bool_job = look_in_file(job,[f"{names.NAME_PRO_SURFTEN}.log"],"Finished",check_for_not=True,check_for_not_str=['Received the TERM', 'Received the INT'])
        return bool_job
    
    
##################################################################################


@FlowProject.label
def data_collected(job):
    test_passed = False
    local_name_of_file = f'{names.GENERAL_GLOBAL_DATA}.txt'
    if os.path.exists(local_name_of_file):
        with open(local_name_of_file, "r") as f:
            contents = f.read()
            if job.id in contents:
                test_passed = True

    return test_passed


@FlowProject.label
def frames_extracted_for_tessellation(job):
    with(job):
        test_passed = False
        pseudolog_file = 'EXTRACT_FRAMES_FOR_TESSELLATION_pseudolog.txt'

        # Check if pseudolog exists
        if not os.path.isfile(pseudolog_file):
            return False

        # Read pseudolog to get nsteps and nstxout
        nsteps = None
        nstxout = None

        with open(pseudolog_file, 'r') as f:
            for line in f:
                if 'nsteps=' in line:
                    parts = line.split(',')
                    for part in parts:
                        if 'nsteps=' in part:
                            nsteps = int(part.split('=')[1].strip())
                        elif 'nstxout=' in part:
                            nstxout = int(part.split('=')[1].strip())

        if nsteps is None or nstxout is None:
            return False

        # Calculate number of frames
        num_frames = int(nsteps / nstxout)

        # Check all directories and their first/last .gro files
        all_dirs_valid = True

        for start_frame in range(0, num_frames, 100):
            end_frame = min(start_frame + 99, num_frames - 1)
            dir_name = f'step_{start_frame}to{end_frame}'

            # Check if directory exists
            if not os.path.isdir(dir_name):
                all_dirs_valid = False
                break

            # Check first .gro file in this directory
            first_gro = f'{dir_name}/production_frame{start_frame}.gro'
            if not os.path.isfile(first_gro):
                all_dirs_valid = False
                break

            # Check last .gro file in this directory
            last_gro = f'{dir_name}/production_frame{end_frame}.gro'
            if not os.path.isfile(last_gro):
                all_dirs_valid = False
                break

        test_passed = all_dirs_valid

        return test_passed

@FlowProject.label
def voronoi_com_built(job):
    with(job):
        # Check if COM extraction is complete
        return os.path.isfile(names.NAME_COM_DATA) and os.path.isfile(names.NAME_COM_LOG)


@FlowProject.label
def voronoi_tess_complete(job):
    with(job):
        # Check if Voronoi tessellation is complete
        return os.path.isfile(names.NAME_VORO_DATA) and os.path.isfile(names.NAME_VORO_LOG)


@FlowProject.label
def voronoi_analysis_complete(job):
    with(job):
        # Check if analysis outputs exist
        required_files = [
            names.NAME_VORO_FINAL_LOG, # 'tessellation_summary.txt',
            names.NAME_VORO_FINAL_RAW, # 'tessellation_volumes_raw.txt',
            names.NAME_VORO_VOL_HIST_PIC, # 'tessellation_volumes_histogram.png',
            names.NAME_VORO_VOL_TIME_PIC, # 'tessellation_volumes_timeseries.png',
            names.NAME_VORO_VOL_2D_PIC, # 'tessellation_volumes_2d.png'
            names.NAME_VORO_STD_PIC
        ]
        return all(os.path.isfile(f) for f in required_files)


###############################################################################
################## MDANALYSIS VORONOI WORKFLOW LABELS #########################
###############################################################################

@FlowProject.label
def mda_voronoi_tessellated(job):
    """Check if Voronoi tessellation is complete (full system)."""
    with(job):
        return os.path.isfile(names.NAME_RESIDUE_DATA_ALL)


@FlowProject.label
def mda_voronoi_classified(job):
    """Check if Voronoi volumes are classified by phase."""
    with(job):
        required_files = [
            names.NAME_RESIDUE_DATA_LIQ,
            names.NAME_RESIDUE_DATA_GAS,
            names.NAME_RESIDUE_DATA_TRA,
            names.NAME_Z_POSITION_HISTOGRAM
        ]
        return all(os.path.isfile(f) for f in required_files)


@FlowProject.label
def mda_individual_histograms_complete(job):
    """Check if individual phase histogram exists."""
    with(job):
        return os.path.isfile(names.NAME_PHASE_HISTOGRAM_INDIVIDUAL)


@FlowProject.label
def slab_aligned(job):
    """Check if slab alignment is complete."""
    with(job):
        required_files = [
            names.NAME_ALIGNED_GRO,
            names.NAME_ALIGNED_TRR,
            # names.NAME_SHIFT_COMPARISON_PNG # uncomment later
        ]
        return all(os.path.isfile(f) for f in required_files)


@FlowProject.label
def phases_classified(job):
    """Check if phase classification is complete."""
    with(job):
        required_files = [
            names.NAME_PHASE_DATA,
            names.NAME_PHASE_Z_COORDS,
            names.NAME_PHASE_CLASSIFICATION_PNG
        ]
        return all(os.path.isfile(f) for f in required_files)


## @FlowProject.label
## def density_results_exported(job):
##     """Check if density results have been exported to the aggregation file."""
##     # This checks if the job's data has been written to the aggregation file
##     # We use a marker file in the job directory to track this
##     with(job):
##         return os.path.isfile('.density_results_exported')


@FlowProject.label
def density_results_exported(job):
    test_passed = False
    local_name_of_file = f'{names.MANY_METHODS_JOBS_DATA}'
    
    if os.path.exists(local_name_of_file):
        #print(f'{names.MANY_METHODS_JOBS_DATA}opened')
        with open(local_name_of_file, "r") as f:
            #print(f'{names.MANY_METHODS_JOBS_DATA} read')
            contents = f.read()
            if job.id in contents:
                test_passed = True
                #print(f'job.id {job.id} test_passed {test_passed}')
        #else:
        #    print(f'{names.MANY_METHODS_JOBS_DATA} not found')

    return test_passed