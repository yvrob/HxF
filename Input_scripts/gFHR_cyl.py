##### Input parameters #####

#%% Calculation options
# Transport/depletion
transport = True
resolve_first = False
correct = False
domain_decomposition = True

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

#%% Case
path_to_case = './Models' # Path to input folder parent
case_name = 'gFHR_nodet' # Name of the input folder
main_input_file = 'input' # Name of the main input file

#%% Serpent data (from Serpent input)

# Pebbles
pbed_file = 'fpb_pos'
r_pebbles = 2 # cm

# Others
pbed_universe_name = 'u_pb'
detector_names = [] #'flux_pebbles_thermal', 'flux_pebbles_fast', 'power_pebbles']

#%% Depletion steps
power_normalization_field = 'power'
power_normalization_value = 280e6 # W
Nsteps = 2000
neutrons_per_cycle = [5000]*134 + [20000]*134 + [200000]*(Nsteps-2*134)
decay_step = 0 # days

#%% Burnup cycle
pebbles_dict = {'u_fuel_pebble':{'mat_name':'fuel', 'pebbles_frac':1, 'r_fuel_kernel':0.02125, 'Ntrisos':9022, 'threshold_type':'passes', 'threshold_dir':+1, 'threshold':8}}

#%% Motion
motion_direction = +1
recirc_threshold = 50 # cm, absolute value

# Discrete motion
if discrete_motion:
    Nrows_to_move = 6
    time_per_pass = 522/8 # days

# DEM
else:
    positions_folder = '/global/scratch/users/yvesrobert/chrono_pbed_case/gFHR_cyl/filtered_settled_gFHR_cyl/'
    DEM_step_increment = 1
    circulation_rate = 250190/(522/8) # pebbles/days
    DEM_circulation_step = 15011 # pebbles
    positions_scale = 100 #  m->cm
    positions_translation = [0,0, 60 - 0.01 - (370-369.47)] # reflector of 60 cm - 1 mm due to DEM - error due to the fact that processing took 310 and not 309.47
    DEM_start = 68
    DEM_end = 201
    DEM_files_prefix = 'filtered_step'
    import numpy as np
    base_reordering = np.loadtxt('/global/scratch/users/yvesrobert/HxF_dev/Models/gFHR_nodet/base_reordering.txt', dtype=int)
    if looping:
    # Looper
        looper_Nr = 7
        looper_Nz = 15
        looper_method = 'rz'

#%% Outputing
output_folder_name = 'gFHR_dem_cyl' # Name of the output folder
verbosity_mode = 3
inventory_names = []

#%% Domain decomposition
if domain_decomposition:
    allowed_decomposition_types = 'rs'
    max_domains = [6, 8]

#%% Restart write
if write_restart:
    restart_write_every = 10

#%% Restart read
if restart_calculation:
    restart_step = 0
    restart_binary = '/global/scratch/users/yvesrobert/HxF_dev/Cases/gFHR_dm_restart_FTC/wrk_Serpent/input.wrk_250' #'/global/scratch/users/yvesrobert/HxF_dev/Cases/gFHR_dm_restart/wrk_Serpent/spectrum_restart/final_gFHR_dm_equilibrium.wrk'
    restart_data =   '/global/scratch/users/yvesrobert/HxF_dev/Cases/gFHR_dm_restart_FTC/Data/core_250.csv' #'/global/scratch/users/yvesrobert/HxF_dev/Cases/gFHR_dm_restart/Data/first_data.csv'
    different_positions = True
elif read_first_compositions:
    restart_binary = '/global/scratch/users/yvesrobert/HxF_dev/Cases/gFHR_dm_restart_FTC/wrk_Serpent/input.wrk_250' #'/global/scratch/users/yvesrobert/HxF_dev/Cases/gFHR_dm_restart/wrk_Serpent/spectrum_restart/final_gFHR_dm_equilibrium.wrk'
    restart_data =   '/global/scratch/users/yvesrobert/HxF_dev/Cases/gFHR_dm_restart_FTC/Data/core_250.csv' #'/global/scratch/users/yvesrobert/HxF_dev/Cases/gFHR_dm_restart/Data/first_data.csv'

if saving:
    write_global = True
    write_incore = True
    write_reinserted = False
    write_discarded = True
    write_discharged = True

if plotting:
    plot_base = True
    plot_detectors = False
    plot_inventory = False
    plot_keff = True
    plot_cumulative = True
    plot_pass = True
    plot_geom = True

extra_fields = ['burnup']

