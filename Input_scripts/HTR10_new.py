##### Input parameters #####

#%% Calculation options
# Transport/depletion
transport = True
resolve_first = False
correct = False
domain_decomposition = False

# Motion
discrete_motion = False
looping = True

# Restart
restart_calculation = True
read_first_compositions = False

# Output
plotting = True
saving = True
write_restart = True
write_restart_discharged = True
write_restart_discarded = False

#%% Case
path_to_case = './Models' # Path to input folder parent
case_name = 'HTR-10_long_det' # Name of the input folder
main_input_file = 'input.inp' # Name of the main input file

#%% Serpent data (from Serpent input)

# Pebbles
pbed_file = 'HTR-10_0.inp'
r_pebbles = 3 # cm
pebbles_dict = {'u_fuel_pebble':{'mat_name':'fuel', 'pebbles_frac':0.57, 'r_fuel_kernel':0.025, 'Ntrisos':8335, 'threshold_type':'burnup', 'threshold_dir':+1, 'threshold':72},
                'u_graph_pebble':{'pebbles_frac':0.43}}

# Others
pbed_universe_name = 'u_pb'
detectors = {} #'flux_pebbles_thermal':{'E':[1e-11, 1.86e-6], 'normalized':True},
            #  'flux_pebbles_fast':{'E':[0.1, 20], 'normalized':True},
            #  'power_pebbles':{'extra_cards':['dr', -8, 'void']}}

#%% Depletion steps
power_normalization_field = 'power'
power_normalization_value = 10e6 # W
Nsteps = 1000
neutrons_per_cycle = 20000
decay_step = 3 # days

#%% Burnup cycle

#%% Motion
motion_direction = -1

# Discrete motion
if discrete_motion:
    Nrows_to_move = 6
    time_per_pass = 522/8 # days

# DEM
else:
    recirc_threshold = 200 # cm, absolute value
    positions_folder = '/global/scratch/users/co_nuclear/pebble_positions_larger/'
    DEM_step_increment = 1
    circulation_rate = 125 # pebbles/days
    DEM_circulation_step = 1564 # pebbles
    positions_scale = 100
    positions_translation = [0,0, -610]
    if looping:
    # Looper
        DEM_start = 40
        DEM_end = 129
        looper_Nr = 5
        looper_Nz = 10
        looper_method = 'rz'
import numpy as np
base_reordering = np.loadtxt('/global/scratch/users/yvesrobert/HxF_dev/Models/HTR-10_long_det/base_reordering.txt', dtype=int)

#%% Outputing
output_folder_name = 'HTR10' # Name of the output folder
verbosity_mode = 0
inventory_names = []

#%% Domain decomposition
if domain_decomposition:
    allowed_decomposition_types = 'rs'
    max_domains = [6, 8]

#%% Restart write
if write_restart:
    restart_write_every = 1
if write_restart_discharged:
    restart_discharged_write_every = 1
if write_restart_discarded:
    restart_discarded_write_every = 5

#%% Restart read
if restart_calculation:
    restart_step = 0
    restart_data = '/global/scratch/users/yvesrobert/HxF/Cases/HTR-10_large_3/Data/core_750.csv'
    restart_binary = '/global/scratch/users/yvesrobert/HxF/Cases/HTR-10_large_3/wrk_Serpent/input.inp.wrk_750'
    different_positions = True
# elif read_first_compositions:
#     restart_binary = '/global/scratch/users/yvesrobert/HxF_dev/Cases/gFHR_dm_restart_FTC/wrk_Serpent/input.wrk_250' #'/global/scratch/users/yvesrobert/HxF_dev/Cases/gFHR_dm_restart/wrk_Serpent/spectrum_restart/final_gFHR_dm_equilibrium.wrk'
#     restart_data =   '/global/scratch/users/yvesrobert/HxF_dev/Cases/gFHR_dm_restart_FTC/Data/core_250.csv' #'/global/scratch/users/yvesrobert/HxF_dev/Cases/gFHR_dm_restart/Data/first_data.csv'

if saving:
    write_global = True
    write_incore = True
    write_reinserted = False
    write_discarded = True
    write_discharged = True

if plotting:
    plot_base = False
    plot_detectors = False
    plot_inventory = False
    plot_keff = True
    plot_cumulative = False
    plot_pass = True
    plot_geom = False

extra_fields = ['burnup']