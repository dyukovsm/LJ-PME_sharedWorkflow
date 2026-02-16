# GMX_PREFIX = '/usr/local/gromacs/bin/gmx' # potoff
GMX_PREFIX = 'gmx' # grid

NAME_EQ_NVT = 'EQ_NVT'
NAME_EQ_SURFTEN = 'EQ_SURFTEN'; NAME_PRO_SURFTEN = 'PRO_SURFTEN'

NAME_ELONGATED = 'ELONGATED_BOX_PLACEHOLDER'

GENERAL_LOCAL_DATA = 'raw_general_data_for'
GENERAL_GLOBAL_DATA = 'aggregate_general_Data'

DENS_LOCAL_DATA = 'raw_dens_Data'
DENS_GLOBAL_DATA = 'aggregate_dens_Data_fromMIN-MAX'

# Slab alignment
NAME_ALIGNED_GRO = 'aligned_PRO_SURFTEN.gro'
NAME_ALIGNED_TRR = 'aligned_PRO_SURFTEN.trr'
NAME_SHIFT_COMPARISON_PNG = 'slab_shift_comparison.png'
NAME_COM_DRIFT_PNG = 'com_drift_post_shift.png'
NAME_PROFILE_DRIFT_PNG = 'profile_drift_post_shift.png'
NAME_ALIGNMENT_ANALYSIS = 'drift_log.txt'

# Dual-tanh fit results
DUAL_TANH_FIT_PNG = 'dual_tanh_density_fit.png'
DUAL_TANH_RESULTS_TXT = 'dual_tanh_results.txt'

# Generalized Gaussian fit results
GEN_GAUSS_FIT_PNG = 'gen_gaussian_density_fit.png'
GEN_GAUSS_RESULTS_TXT = 'gen_gaussian_results.txt'