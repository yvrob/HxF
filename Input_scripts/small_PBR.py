##### Input parameters #####

#%% Calculation options
# Transport/depletion
transport = True
resolve_first = False

# Motion
discrete_motion = True

# Output
plotting = True
saving = True
write_restart = True
write_restart_discharged = False
write_restart_discarded = False

#%% Case
path_to_case = './Models' # Path to input folder parent
case_name = 'small_PBR' # Name of the input folder
main_input_file = 'input' # Name of the main input file

#%% Serpent data (from Serpent input)

# Pebbles
pbed_file = 'fpb_pos'
r_pebbles = 2 # cm

# Others
pbed_universe_name = 'u_pb'
detectors = {'flux_pebbles_thermal':{'E':[1e-11, 1.86e-6], 'normalized':True},
             'flux_pebbles_fast':{'E':[0.1, 20], 'normalized':True},
             'power_pebbles':{'extra_cards':['dr', -8, 'void']}}

#%% Depletion steps
power_normalization_field = 'power'
power_normalization_value = 3e6 # W
Nsteps = 1000
neutrons_per_cycle = 100 #00
decay_step = 0 # days

#%% Burnup cycle
pebbles_dict = {'u_fuel_pebble':{'mat_name':'fuel', 'pebbles_frac':1, 'r_fuel_kernel':0.02125, 'Ntrisos':9022, 
                'threshold_type':'burnup', 'threshold_dir':+1, 'threshold':'adjustable', 'threshold_ini':120,
                'target':'keff', 'start_search':0.04, 'target_value':1.00,
                'max_threshold_delta':5, 'learning_rate':200, 'decay_rate':0.99, 'update_every':1}}

#%% Motion
motion_direction = +1
recirc_threshold = 1 # cm, absolute value

# Discrete motion
if discrete_motion:
    Nrows_to_move = 5
    time_per_pass = 100 # days

#%% Outputing
output_folder_name = 'small_PBR_test_det' # Name of the output folder
verbosity_mode = 3
inventory_names = [] #'922350', '922360', '922380', '942390', '942380', '942400', '942410', '942420', '952410', '952421', '551370', '541350', '531310', '380900']

read_first_compositions=False
restart_binary = '/global/scratch/users/yvesrobert/HxF_dev/Cases/small_PBR/first_compos.wrk'
restart_data = '/global/scratch/users/yvesrobert/HxF_dev/Cases/small_PBR/first_data.csv'

#%% Restart write
if write_restart:
    restart_write_every = 10
if write_restart_discharged:
    restart_discharged_write_every = 1
if write_restart_discarded:
    restart_discarded_write_every = 5

if saving:
    write_global = True
    write_incore = True
    write_reinserted = False
    write_discarded = True
    write_discharged = True

if plotting:
    plot_base = True
    plot_detectors = True
    plot_inventory = False
    plot_keff = True
    plot_cumulative = True
    plot_pass = True

extra_fields = ['burnup']