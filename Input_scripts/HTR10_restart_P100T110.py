##### Input parameters #####

#%% Calculation options
# Transport/depletion
transport = True
resolve_first = True
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

#%% Case
path_to_case = './Models' # Path to input folder parent
case_name = 'HTR-10' # Name of the input folder
main_input_file = 'input.inp' # Name of the main input file

#%% Serpent data (from Serpent input)

# Pebbles
pbed_file = 'HTR-10_0.inp'
r_pebbles = 3 # cm
fuel_mat_name = 'fuel'
fuel_frac = 0.57
fuel_pebble_universe_name = 'u_fuel_pebble'
graph_pebble_universe_name = 'u_graph_pebble'

# Others
pbed_universe_name = 'u_pb'
detector_names = ['flux_pebbles_thermal', 'flux_pebbles_fast', 'power_pebbles']

#%% Depletion steps
power_normalization_field = 'power'
power_normalization_value = 10e6 # W
Nsteps = 1000
neutrons_per_cycle = 200000  #*20 + [10000]*20 + [10000]*(Nsteps-20-20-50) + [50000]*50
decay_step = 0 # days

#%% Burnup cycle
threshold_type = 'burnup'
threshold = 72*1.1

#%% Motion
motion_direction = -1
recirc_threshold = 200 # cm, absolute value

# Discrete motion
if discrete_motion:
    Nrows_to_move = [20]*20+[10]*20+[5]*(Nsteps-20-20)
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
output_folder_name = 'HTR-10_restart_P1T1.1' # Name of the output folder
verbosity_mode = 3


#%% Domain decomposition
if domain_decomposition:
    allowed_decomposition_types = 'ars'
    max_domains = [6, 4, 10]

#%% Restart write
if write_restart:
    restart_write_every = 10

#%% Restart read
if restart_calculation:
    restart_step = 650
    restart_data = '/global/scratch/users/yvesrobert/HTR-10_latest/Cases/HTR-10_P1_T1/Data/data_650.csv'
    restart_binary = '/global/scratch/users/yvesrobert/HTR-10_latest/Cases/HTR-10_P1_T1/wrk_Serpent/input.inp.wrk_650'
elif read_first_compositions:
    restart_binary = './Tools/Create_initial_compositions/restart/first_compos.wrk'
    restart_data = './Tools/Create_initial_compositions/data_first_compos.csv'

