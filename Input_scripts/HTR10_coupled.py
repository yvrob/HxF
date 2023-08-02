##### Input parameters #####

#%% Calculation options
# Transport/depletion
transport = True
resolve_first = False
correct = False
domain_decomposition = False
use_decnfy_lib = True

# Thermal_hydraulics
thermal_coupling = True

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
write_restart_discarded = True

#%% Case
path_to_case = './Models' # Path to input folder parent
case_name = 'HTR-10_long_coupled' # Name of the input folder  A CHANGER
main_input_file = 'input.inp' # Name of the main input file

#%% Serpent data (from Serpent input)

# Pebbles
pbed_file = 'HTR-10_0.inp'
r_pebbles = 3 # cm
pebbles_dict = {'u_fuel_pebble':{'mat_name':'fuel', 'pebbles_frac':0.57, 'r_fuel_kernel':0.025, 'Ntrisos':8335, 'threshold_type':'burnup', 'threshold_dir':+1, 'threshold':72},
                'u_graph_pebble':{'pebbles_frac':0.43}}

# Others
pbed_universe_name = 'u_pb'
detectors = {'flux_pebbles_thermal':{'E':[1e-11, 1.86e-6], 'normalized':True},
             'flux_pebbles_fast':{'E':[0.1, 20], 'normalized':True},
             'power_pebbles':{'extra_cards':['dr', -8, 'void']}}

#%% Depletion steps
power_normalization_field = 'power'
power_normalization_value = 10e6 # W
Nsteps = 10000
neutrons_per_cycle = 50000  #*20 + [10000]*20 + [10000]*(Nsteps-20-20-50) + [50000]*50
decay_step = 3 # days

#%% Burnup cycle
inventory_names = ['551370'] #'10010', '10030', '20040', '30070', '40100', '50100', '50110', '60120', '60140', '80160', '80170', '290670', '300660', '300670', '300680', '300700', '300720', '310690', '310710', '310720', '320720', '320730', '320740', '320760', '320770', '320780', '330750', '330760', '330770', '330780', '330790', '330810', '340760', '340770', '340780', '340790', '340791', '340800', '340810', '340811', '340820', '340830', '340840', '340850', '340860', '350790', '350810', '350820', '350830', '350840', '350841', '350850', '350860', '350870', '350880', '360800', '360810', '360820', '360830', '360831', '360840', '360850', '360851', '360860', '360870', '360880', '360890', '360900', '360910', '370830', '370840', '370850', '370860', '370870', '370880', '370890', '370900', '370901', '370910', '370920', '370930', '380860', '380870', '380880', '380890', '380900', '380910', '380920', '380930', '380940', '380950', '390880', '390890', '390891', '390900', '390910', '390911', '390920', '390930', '390940', '390950', '390960', '390961', '390970', '400900', '400910', '400920', '400930', '400940', '400950', '400960', '400970', '400980', '400990', '401000', '401010', '401020', '410930', '410931', '410940', '410950', '410951', '410960', '410970', '410980', '410981', '410990', '410991', '411000', '411010', '420940', '420950', '420960', '420970', '420980', '420990', '421000', '421010', '421020', '421030', '421040', '421050', '421060', '430980', '430990', '430991', '431000', '431010', '431020', '431030', '431040', '431050', '431060', '431070', '440990', '441000', '441010', '441020', '441030', '441040', '441050', '441060', '441070', '441080', '441090', '451020', '451021', '451030', '451031', '451040', '451041', '451050', '451051', '451060', '451061', '451070', '451080', '451090', '461040', '461050', '461060', '461070', '461080', '461090', '461100', '461110', '461120', '471090', '471091', '471101', '471110', '471111', '471120', '471130', '471150', '481100', '481110', '481120', '481130', '481131', '481140', '481150', '481151', '481160', '481170', '481171', '481180', '491130', '491150', '491151', '491170', '491171', '491191', '501150', '501160', '501170', '501171', '501180', '501190', '501191', '501200', '501210', '501211', '501220', '501230', '501231', '501240', '501250', '501251', '501260', '501270', '501271', '501280', '501290', '501291', '501300', '501301', '501310', '501311', '501320', '511210', '511220', '511230', '511240', '511250', '511260', '511261', '511270', '511280', '511281', '511290', '511300', '511301', '511310', '511320', '511321', '511330', '521220', '521230', '521231', '521240', '521250', '521251', '521260', '521270', '521271', '521280', '521290', '521291', '521300', '521310', '521311', '521320', '521330', '521331', '521340', '521350', '521360', '531260', '531270', '531280', '531290', '531300', '531301', '531310', '531320', '531321', '531330', '531340', '531341', '531350', '531360', '531361', '531370', '531380', '541280', '541290', '541300', '541310', '541311', '541320', '541330', '541331', '541340', '541350', '541351', '541360', '541370', '541380', '541390', '541400', '551320', '551330', '551340', '551341', '551350', '551351', '551360', '551370', '551380', '551381', '551390', '551400', '551410', '561320', '561340', '561350', '561360', '561370', '561371', '561380', '561390', '561400', '561410', '561420', '561430', '561440', '561450', '571370', '571380', '571390', '571400', '571410', '571420', '571430', '571440', '571450', '571460', '571461', '581380', '581390', '581400', '581410', '581420', '581430', '581440', '581450', '581460', '581470', '581480', '591410', '591420', '591421', '591430', '591440', '591441', '591450', '591460', '591470', '591480', '591481', '591490', '591510', '601420', '601430', '601440', '601450', '601460', '601470', '601480', '601490', '601500', '601510', '601520', '601530', '611460', '611470', '611480', '611481', '611490', '611510', '611520', '611530', '611540', '621470', '621480', '621490', '621500', '621510', '621520', '621530', '621540', '621550', '621560', '621570', '621580', '631510', '631520', '631530', '631540', '631541', '631550', '631560', '631570', '631580', '631590', '641520', '641540', '641550', '641560', '641570', '641580', '641590', '641600', '651580', '651590', '651600', '651610', '661600', '661610', '661620', '661630', '661640', '661660', '671650', '671660', '671661', '681660', '681670', '681680', '681690', '681700', '691690', '691710', '701720', '902310', '902320', '902340', '912310', '912340', '922320', '922340', '922350', '922360', '922370', '922380', '922390', '932370', '932380', '932390', '942380', '942390', '942400', '942410', '942420', '952410']

#%% Thermal-hydraulics
if thermal_coupling:
    plot_thermal = True
    TH = {'solver': 'GeN-Foam', 'step_size':5, 'max_steps':1, 'nnodes':2, 'time_limit':1, 'positions_scale':100, 'fuel_mat': 'fuel',
    'fields_of_interest':['Q', 'Tfav.nuclearSteadyStatePebble', 'Tmav.nuclearSteadyStatePebble', 'T']}#    'convergence_criteria':{'Q':0.03, 'Tfav.nuclearSteadyStatePebble':0.01, 'Tmav.nuclearSteadyStatePebble':0.01, 'T':0.01, 'keff':30e-5}}

#%% Motion

# Discrete motion
if discrete_motion:
    motion_direction = -1
    Nrows_to_move = [20]*20+[10]*20+[5]*(Nsteps-20-20)
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
    DEM_start = 40
    DEM_end = 129
    if looping:
    # Looper

        looper_Nr = 5
        looper_Nz = 10
        looper_method = 'rz'

#%% Outputing
output_folder_name = 'HTR10_coupled2' # Name of the output folder
verbosity_mode = 1

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
    restart_discarded_write_every = 1

import numpy as np
base_reordering = np.loadtxt('/global/scratch/users/yvesrobert/HxF_dev/Models/HTR-10_long_det/base_reordering.txt', dtype=int)

#%% Restart read
if restart_calculation:
    restart_step = 1435
    restart_data = '/global/scratch/users/yvesrobert/HxF_dev/Cases/HTR10_restart/Data/core_1435.csv'
    restart_binary = '/global/scratch/users/yvesrobert/HxF_dev/Cases/HTR10_restart/wrk_Serpent/input.inp.wrk_1435'
    different_positions = True
elif read_first_compositions:
    # restart_binary = '/global/scratch/users/clementvayrondelamoureyre/HTR10_large/input.inp.wrk_750'
    # restart_data = '/global/scratch/users/clementvayrondelamoureyre/HTR10_large/core_750.csv'
    restart_binary ='/global/scratch/users/clementvayrondelamoureyre/HTR10_large/input.inp.wrk_750'
    restart_data ='/global/scratch/users/clementvayrondelamoureyre/HTR10_large/core_750.csv'


if saving:
    write_global = True
    write_incore = True
    write_reinserted = True
    write_discarded = True
    write_discharged = True