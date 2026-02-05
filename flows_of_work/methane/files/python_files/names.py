# GMX_PREFIX = '/usr/local/gromacs/bin/gmx' # potoff
GMX_PREFIX = 'gmx' # grid

NAME_EQ_NVT = 'EQ_NVT'
NAME_EQ_SURFTEN = 'EQ_SURFTEN'; NAME_PRO_SURFTEN = 'PRO_SURFTEN'

NAME_ELONGATED = 'ELONGATED_BOX_PLACEHOLDER'

GENERAL_LOCAL_DATA = 'raw_general_data_for'
GENERAL_GLOBAL_DATA = 'aggregate_general_Data'

NAME_PSEUDOLOG_TESS_EXTRACT = 'EXTRACT_FRAMES_FOR_TESSELLATION_pseudolog.txt'
NAME_COM_DATA = 'voronoi_centers_of_mass.npz'; NAME_COM_LOG = 'voronoi_com_metadata.txt'
NAME_VORO_DATA = 'voronoi_volumes.npz'; NAME_VORO_LOG = 'voronoi_diagnostics.txt'

NAME_VORO_FINAL_LOG = 'tessellation_summary.txt'; NAME_VORO_FINAL_RAW = 'tessellation_volumes_raw.txt'
NAME_VORO_VOL_HIST_PIC = 'tessellation_volumes_histogram.png'; NAME_VORO_VOL_TIME_PIC = 'tessellation_volumes_timeseries.png'
NAME_VORO_VOL_2D_PIC = 'tessellation_volumes_2d.png'; NAME_VORO_STD_PIC = 'tessellation_standard_deviation.png'

# MDAnalysis-based Voronoi workflow files
NAME_ALIGNED_GRO = 'slab_alignedWith_halfZlen_PRO_SURFTEN.gro'
NAME_ALIGNED_TRR = 'slab_alignedWith_halfZlen_PRO_SURFTEN.trr'
NAME_RESIDUE_DATA_ALL = 'residue_analysis_data_ALL.npz'
NAME_RESIDUE_DATA_LIQ = 'residue_analysis_data_LIQ.npz'
NAME_RESIDUE_DATA_GAS = 'residue_analysis_data_GAS.npz'
NAME_RESIDUE_DATA_TRA = 'residue_analysis_data_TRA.npz'
NAME_Z_POSITION_HISTOGRAM = 'z_position_distribution_by_phase.png'
NAME_PHASE_HISTOGRAM_INDIVIDUAL = 'volume_histogram_individual_phases.png'
NAME_PHASE_HISTOGRAM_LIQUID = 'volume_histogram_liquid_cutoff_{cutoff}nm.png'
NAME_PHASE_HISTOGRAM_GAS = 'volume_histogram_gas_cutoff_{cutoff}nm.png'
NAME_PHASE_HISTOGRAM_GAS_LOG = 'log_scale_volume_histogram_gas_cutoff_{cutoff}nm.png'
NAME_PHASE_HISTOGRAM_TRANSITION = 'volume_histogram_transition_cutoff_{cutoff}nm.png'

# Constants for MDAnalysis workflow
MAX_FRAMES_TO_READ = 2000
MOLAR_MASS_WATER = 18.015; MOLAR_MASS_METHANE = 16.04 # g/mol
CURRENT_MOLAR_MASS = MOLAR_MASS_METHANE
HISTOGRAM_BINS = 200
GAS_HIST_CAP = 10000
LIQ_HIST_CAP = 200
TRA_HIST_CAP = 200

# Density profile and slab alignment constants
BIN_SIZE = 0.45  # nm - bin size for density profile
D_TANH_THICKNESS_MULTIPLIER = 1.0  # Multiplier for interface thickness
DERIVATIVE_NOISE_THRESHOLD_MULTIPLIER = 2.0  # For classifying stable vs transition regions
NAME_DENSITY_XVG = 'density_profile.xvg'
NAME_DENSITY_XVG_POST = 'density_profile_post_shift.xvg'
NAME_SHIFT_COMPARISON_PNG = 'PRE_POST_SHIFT_comparison.png'
NAME_TANH_FIT_PNG = 'ROTATED_UNROTATED_tanh_fit_comparison.png'
NAME_PHASE_CLASSIFICATION_PNG = 'phase_classification_with_derivative.png'
NAME_PHASE_DATA = 'slab_phase_classification_data.txt'
NAME_PHASE_Z_COORDS = 'slab_phase_z_coordinates.npz'

# Multi-method density results aggregation
MANY_METHODS_JOBS_DATA = 'density_results.txt'

# Chunk-based workflow is no longer used for methane
# Direct file names (NAME_EQ_SURFTEN, NAME_PRO_SURFTEN) are used instead

