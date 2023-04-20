##### Input parameters #####

#%% Calculation options
# Transport/depletion
transport = True
resolve_first = True
correct = False
domain_decomposition = True

# Motion
discrete_motion = True
looping = False

# Restart
restart_calculation = False
read_firt_compositions = True

# Output
plotting = True
saving = True
write_restart = True

#%% Case
path_to_case = './Models' # Path to input folder parent
case_name = 'gFHR' # Name of the input folder
main_input_file = 'input' # Name of the main input file

#%% Serpent data (from Serpent input)

# Pebbles
pbed_file = 'fpb_pos'
r_pebbles = 2 # cm
fuel_mat_name = 'fuel'
fuel_frac = 1
fuel_pebble_universe_name = 'u_fuel_pebble'
graph_pebble_universe_name = 'u_graph_pebble'

# Others
pbed_universe_name = 'u_pb'
detector_names = ['flux_pebbles_thermal', 'flux_pebbles_fast', 'power_pebbles']

#%% Depletion steps
power_normalization_field = 'power'
power_normalization_value = 280e6 # W
Nsteps = 1000
neutrons_per_cycle = [20000]*68*2+[500000]*(Nsteps-68*3)
decay_step = 0 # days

#%% Burnup cycle
threshold_type = 'passes'
threshold = 8

#%% Motion
motion_direction = +1
recirc_threshold = 1 # cm, absolute value

# Discrete motion
if discrete_motion:
    Nrows_to_move = [6]*68*2+[3]*(Nsteps-68*3)
    time_per_pass = 522/8 # days

# DEM
else:
    positions_folder = '/global/scratch/users/yvesrobert/HTR-10/pebble_positions/'
    DEM_step_increment = 1
    circulation_rate = 125 # pebbles/days
    DEM_circulation_step = 1350 # pebbles
    positions_scale = 100
    positions_translation = [0,0, -610]
    if looping:
    # Looper
        DEM_start = 60
        DEM_end = 244
        looper_Nr = 5
        looper_Nz = 10
        looper_method = 'xyz'

#%% Outputing
output_folder_name = 'gFHR' # Name of the output folder
verbosity_mode = 1
inventory_names = ['922350', '922360', '922380', '942390', '942380', '942400', '942410', '942420', '952410', '952421', '551370', '541350', '531310', '380900']

#%% Domain decomposition
if domain_decomposition:
    allowed_decomposition_types = 'rs'
    max_domains = [6, 8]

#%% Restart write
if write_restart:
    restart_write_every = 5

#%% Restart read
if restart_calculation:
    restart_step = 650
    restart_data = '/global/scratch/users/yvesrobert/HTR-10_latest/Cases/HTR-10_P1_T1/Data/data_650.csv'
    restart_binary = '/global/scratch/users/yvesrobert/HTR-10_latest/Cases/HTR-10_P1_T1/wrk_Serpent/input.inp.wrk_650'
elif read_firt_compositions:
    restart_binary = '/global/scratch/users/yvesrobert/HxF/Tools/Create_initial_compositions/restart/first_compos.wrk'
    restart_data = '/global/scratch/users/yvesrobert/HxF/Tools/Create_initial_compositions/data_first_compos.csv'

if saving:
    write_global = True
    write_incore = True
    write_reinserted = True
    write_discarded = True
    write_discharged = True

