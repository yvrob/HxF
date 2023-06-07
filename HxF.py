# Known issues:
# - Restart, first step is completely wrong, but then it works fine (no impact seen)
# - Domain decomposition, if uncommented, takes too long to update domains -> changes to Serpent source code to do

print("""
============================================================
||                                                        ||
||    8 8888        8 `8.`8888.      ,8' 8 88888888888    ||
||    8 8888        8  `8.`8888.    ,8'  8 8888           ||
||    8 8888        8   `8.`8888.  ,8'   8 8888           ||
||    8 8888        8    `8.`c8888.,8'    8 8888           ||
||    8 8888        8     `8.`88888'     8 88888888888    ||
||    8 8888        8     .88.`8888.     8 8888           ||
||    8 8888888888888    .8'`8.`8888.    8 8888           ||
||    8 8888        8   .8'  `8.`8888.   8 8888           ||
||    8 8888        8  .8'    `8.`8888.  8 8888           ||
||    8 8888        8 .8'      `8.`8888. 8 8888           ||
||                                                        ||
============================================================
""")

from datetime import datetime
from dateutil.tz import tzlocal

#%% Printing time
tz = tzlocal() # time zone
now = datetime.now(tz) # date
date_str = now.strftime('%A, %B %d, %Y at %I:%M %p %Z')
print(f"Starting calculation on: {date_str}\n")

#%% Import pebble bed library and utilities
print('Importing modules\n')
from Source.Pebble_Bed import *
from Source.Utilities import *
from Source.Looping import *

import pandas as pd
import os
from glob import glob
import numpy as np
import importlib
import atexit
import sys
from copy import deepcopy
import inspect
import shutil

#%% Read input from file
print('Reading input\n')
filename_input = sys.argv[1]
sys.path.append('Source')
sys.path.append('Utils')
default_parameters = importlib.import_module('Default_Input.py'.replace(".py", "").replace('./', '')).__dict__
globals().update(default_parameters) # Will read all defaults parameters in input file and import them here
parameters = {k: v for k, v in default_parameters.items() if not (k.startswith('__') or k.endswith('__') or callable(v))}
sys.path.append(os.path.dirname(filename_input))
input_parameters = importlib.import_module(os.path.basename(filename_input).replace(".py", "").replace('./', '')).__dict__
globals().update(input_parameters) # Will read all parameters in input file and import them here (overwrites default)
parameters.update({k: v for k, v in input_parameters.items() if not (k.startswith('__') or k.endswith('__') or callable(v))})

try:
    ncores = int(sys.argv[2])
    parameters['ncores'] = ncores
    if len(sys.argv)>2:
        nnodes = int(sys.argv[3])
        parameters['nnodes'] = nnodes
except:
    pass

# Print
print("Simulation Input Parameters Summary:")
print("-" * 50)
bools = {}
others = {}
paths = {}
for key, value in sorted(parameters.items()):
    if isinstance(value, bool):
        bools[key] = value
    elif isinstance(value, str) and ('/' in value or '\\' in value):
        paths[key] = value
    else:
        others[key] = value
for key, value in sorted(bools.items(), key=lambda x: x[1], reverse=True):
    print("{:<40}{}".format(key, value))
print()
for key, value in sorted(others.items()):
    print("{:<40}{}".format(key, value))
print()
for key, value in sorted(paths.items()):
    print("{:<40}{}".format(key, value))
print("-" * 50)

# Constants
DEPLETION = 1
DECAY = 2
DAYS = 86400

#%% Prepare the case

original_path = str(os.getcwd())
init_case(case_name, filename_input, path_to_case, output_folder_name) # Copy input files in output folder
os.chdir(f'Cases/{output_folder_name}') # Go to output folder
atexit.register(os.chdir, original_path) # Come back to original folder when the script ends (even if error)

#%% Inventory
inventory_names  = process_inventory(inventory_names) # translate keywords to isotopes if needed
if inventory_names:
    with open(main_input_file, 'a') as f:
        f.write(f'\nset inventory {" ".join(inventory_names)}\n')
else:
    inventory_names = []

#%% Domain decomposition
if domain_decomposition:
    if nnodes==1: # for debugging
        ndomains=20
    else:
        ndomains=nnodes
    decomposition_types, decomposition_domains = nodes_to_dd(ndomains, allowed_decomposition_types, max_domains)
    print(f'Using Domain Decomposition with domains of types "{decomposition_types}" and breakout {decomposition_domains}')
    # Add line for domain decomposition to the input
    with open(main_input_file, 'a') as f:
        f.write(f'\nset dd 5 "initial_domains.txt"\n') # 5 is file decomposition

#%% Change positions, universes and radii based on input

# If restarting from a step, calculate equivalent step and positions
if transport:
    if restart_calculation:
        first_step = restart_step
        print(f'Restarting from step {restart_step}, binary data at "{restart_binary}" and reading data at "{restart_data}".')
        if not os.path.exists(restart_data):
            raise Exception(f'Restart mode selected but no restart data table found at {restart_data}')
        restart_files = glob(restart_binary+'*')
        if len(restart_files)==0:
            raise Exception(f'Restart mode selected but no restart binary found at {restart_binary}')
        else:
            nrestarts = len(restart_files) # count number of restart files
        with open(main_input_file, 'a') as f:
            f.write(f'\nset rfr idx 0 "{restart_binary}" {nrestarts}\n')
    else:
        first_step = 0
        if read_firt_compositions:
            if not os.path.exists(restart_data):
                raise Exception(f'First composition binary mode selected but no restart data table found at {restart_data}')
            restart_files = glob(restart_binary+'*')
            if len(restart_files)==0:
                raise Exception(f'First composition binary mode selected but no restart binary found at {restart_binary}')
            else:
                nrestarts = len(restart_files) # count number of restart files
            with open(main_input_file, 'a') as f:
                f.write(f'\nset rfr idx 0 "{restart_binary}" {nrestarts}\n')
else:
    first_step = 0

# DEM motion: import first DEM positions
if not discrete_motion:
    print(f'Using DEM motion from folder {positions_folder}')
    if not os.path.exists(positions_folder):
        raise Exception(f'DEM motion selected but positions folder {positions_folder} does not exist.')
    position_files = natural_sort(glob(positions_folder+'/step*.csv')) # List all positions found with DEM
    if looping:
        DEM_start, DEM_end, transition, error_rz = looper(position_files, DEM_start, DEM_end, looper_Nr, looper_Nz, looper_method)
        print(f'Looper ready: from DEM step {DEM_start} to {DEM_end} (err={error_rz:.2f})')
    position_files = position_files[DEM_start:DEM_end+1]

    if restart_calculation or read_firt_compositions:
        data = pd.read_csv(restart_data, index_col=0)
    else:
        data = pd.read_csv(position_files[0])[['x','y','z']]*positions_scale + np.array(positions_translation) # read new positions from DEM file
        data['r'] = r_pebbles # Change radii

    data['r_dist'] = np.linalg.norm(data[['x', 'y']], axis=1)
    time_per_pass = data.shape[0]/circulation_rate

# Discrete motion: import from pbed file, assign row/column ID and make motion matrix
else:
    print(f'Using discrete motion')
    if restart_calculation or read_firt_compositions:
        data = pd.read_csv(restart_data, index_col=0)
    else:
        data = pd.read_csv(pbed_file, header=None, names=['x', 'y', 'z', 'r', 'uni'], delim_whitespace=True)

    columns_group = data.groupby(["y", "x"], sort=False).ngroup().add(1) - 1
    data["column_id"] = columns_group
    rows_group = data.groupby("z")
    zrows = list(rows_group.groups.keys())
    data["row_id"] = rows_group.ngroup().add(1) - 1
    rows_group = data.groupby(data["row_id"] // 2)
    data["row_id"] = rows_group.ngroup().add(1) - 1

    motion_step = motion_direction * (np.diff(zrows).mean()*2)
    if discrete_motion == +1:
        zlim = zrows[-1] + motion_step/4
    else:
        zlim = zrows[0] + motion_step/4
print()

# Assign random universe, based on fuel fraction, if not restarting
if not restart_calculation:
    active_pebbles_dict = dict()
    threshold_pebbles_dict = dict()
    pebbles_fracs = getdict(pebbles_dict, 'pebbles_frac')
    pebbles_fracs /= pebbles_fracs.sum() # normalize fractions
    print(f'Assigning random universes {pebbles_dict.keys()} corresponding to {getdict(pebbles_dict, "mat_name")} with fractions of {pebbles_fracs*100}%\n')
    data['uni'] = assign_random_array(np.arange(data.shape[0]), getdict(pebbles_dict), list(pebbles_fracs)) # Split into fuel and non-fuel
    data['isactive'] = False
    for uni_id, (uni_name, uni) in enumerate(pebbles_dict.items()):
        data[f'pebble_type_{uni_id}'] = (data['uni'] == uni_name)
        if 'mat_name' in uni.keys():
            data.loc[data['uni'] == uni_name, 'mat_name'] = uni['mat_name']
            data.loc[data['uni'] == uni_name, 'isactive'] = True
            active_pebbles_dict[uni_name] = pebbles_dict[uni_name]
        #     if any(str(element) in ['922330', '942380', '902320', '922380', '922350', '942390', '942410'] for element in nuclides_list[uni['mat_name']]):
        #         data.loc[data['uni'] == uni_name, 'isactive'] = True
        #         active_pebbles_dict[uni_name] = pebbles_dict[uni_name]

        if ('threshold_type' in uni.keys() and uni['threshold_type']) or 'threshold' in uni.keys():
            # data.loc[data['uni'] == uni_name, 'isactive'] = True
            threshold_pebbles_dict[uni_name] = pebbles_dict[uni_name]



# Domain decomposition

if domain_decomposition:
    first_pbed = Pebble_bed()
    first_pbed.read_dataframe(data)
    first_pbed.decompose_in_domains(decomposition_domains, decomposition_types)
    first_pbed.data.loc[data['isactive'], 'domain_id'].astype(int).to_csv('initial_domains.txt', index=False, header=False)
    data.loc[data['isactive'], 'domain_id'] = first_pbed.data.loc[data['isactive'], 'domain_id']

# Save
data[['x','y','z','r','uni']].to_csv(pbed_file, header=False, index=False, sep='\t')
data['r_dist'] = np.linalg.norm(data[['x', 'y']], axis=1)
data.loc[list(data.sort_values(['z', 'r_dist']).index), 'id'] = np.arange(len(data)).astype(int) # sort by z and radial distance

# Count
print(f'{len(data)} pebbles in the core:')
for uni_id, (uni_name, uni) in enumerate(pebbles_dict.items()):
    print(f'\t- {(data["uni"]==uni_name).sum()} pebbles with universe {uni_name}')
print('')

# Prepare case
if not use_decnfy_lib: # remove link to decay/nfy library if not needed
    print('Not using decay/nfylibraries')
    os.environ.pop('SERPENT_DECLIB', None)    
    os.environ.pop('SERPENT_NFYLIB', None)    

serpent, serpent_input_files = start_Serpent(os.environ["SERPENT_EXE"], ncores, main_input_file, nnodes, verbosity_mode)
nplots = count_plots(serpent_input_files) # count number of plot commands in input files



#%% Time and filling step-dependent list with values if needed

# Create dictionnary of step-wise variables
step_wise_variables = dict()
for var in ["Nrows_to_move", "neutrons_per_cycle", "power_normalization_value", "DEM_step_increment"]:
    try: # Add only the ones in the input
        step_wise_variables[var] = eval(var)
        try: # If single value, make list of the same value at every step
            _ = len(step_wise_variables[var])
        except:
            step_wise_variables[var] = [step_wise_variables[var]] * Nsteps
    except:
        pass

# Process threshold for each pebble type
for uni_id, (uni_name, uni) in enumerate(threshold_pebbles_dict.items()):
    var = f'threshold_{uni["mat_name"]}'
    try: # Add only the ones in the input
        step_wise_variables[var] = uni['threshold']
        try: # If single value, make list of the same value at every step
            _ = len(step_wise_variables[var])
        except:
            step_wise_variables[var] = [step_wise_variables[var]] * Nsteps
    except:
        step_wise_variables[var] = None


# Add time steps
step_wise_variables["time_step"] = [0] # first time step is 0, only transport
for step in range(1, Nsteps):
    if discrete_motion:
        step_wise_variables["time_step"].append(step_wise_variables["Nrows_to_move"][step-1] * (time_per_pass/(data.row_id.max() + 1))) # days)
    else:
        step_wise_variables["time_step"].append(step_wise_variables["DEM_step_increment"][step-1] * DEM_circulation_step/circulation_rate) # days

# Test if right number of values
for var in step_wise_variables:
    if len(step_wise_variables[var]) != Nsteps:
        raise Exception(f'Error: {var} should have as many values as the number of steps ({len(step_wise_variables[var])} vs {Nsteps})')

# Initialize time
curr_time = np.sum(step_wise_variables["time_step"][:first_step+1]).round(6)*DAYS # Initial time step
print(f'Setting Serpent time to {curr_time/DAYS:.2f} days.')
serpent.set_current_time(curr_time) # Communicate to Serpent

#%% Transfer information

# Get from Serpent once

# Find one material zone per universe name, and get composition
serpent_materials = natural_sort(Serpent_get_values('materials', serpent))
nuclides_list = {mat_name: [] for mat_name in getdict(active_pebbles_dict, 'mat_name')}
fresh_fuel = {mat_name: [] for mat_name in getdict(active_pebbles_dict, 'mat_name')}
for mat_name in getdict(active_pebbles_dict, 'mat_name'):
    for serpent_mat in serpent_materials:
        if serpent_mat.startswith(f'{mat_name}z'):
            nuclides_list[mat_name] = list(Serpent_get_values(f'composition_{serpent_mat}_zai', serpent))
            fresh_fuel[mat_name] = list(Serpent_get_values(f'composition_{serpent_mat}_adens', serpent))
            break

#%% Create transferrables, which will be callable during calculation (Serpent <=> Python)

# Transferrables will be in the "tra" dictionnary
tra = dict()

##### Python => Serpent #####
tra['plot']          = Serpent_set_values("plot_geometry", 1, serpent, communicate=False) # Can be called to plot latest position/geometry
tra['write_restart'] = get_transferrable("write_restart", serpent, input_parameter=True) # Can be called to plot latest position/geometry
tra['xyzr_in']       = get_transferrable(f"pbed_{pbed_universe_name}_xyzr", serpent, input_parameter=True) # Can be called to change pebbles positions
tra['burnup_in'] = {mat_name: Serpent_get_material_wise(mat_name, 'burnup', serpent, input_parameter=True) for mat_name in getdict(active_pebbles_dict, 'mat_name')}
tra['reset_fuel'] = {mat_name: Serpent_get_material_wise(mat_name, 'reset', serpent, prefix='composition', input_parameter=True) for mat_name in getdict(active_pebbles_dict, 'mat_name')}
tra['burnable_fuel'] = {mat_name: Serpent_get_material_wise(mat_name, 'burnable', serpent, prefix='material', input_parameter=True) for mat_name in getdict(active_pebbles_dict, 'mat_name')}

if domain_decomposition:
    tra['domain_in'] = {mat_name: get_transferrable(f"material_div_{mat_name}_domain", serpent, input_parameter=True) for mat_name in getdict(active_pebbles_dict, 'mat_name')}

for uni_id, (uni_name, uni) in enumerate(threshold_pebbles_dict.items()):
    if uni['threshold_type'] != 'passes':
        tra[f"{uni['threshold_type']}_{uni['mat_name']}_in"] = get_transferrable(f"material_div_{uni['mat_name']}_{uni['threshold_type']}", serpent, input_parameter=True)

tra['switch_mode']   = get_transferrable('burn_step_type', serpent, input_parameter=True) # Can be called to switch between decay mode and depletion mode (useful for decay)
tra['time_in']       = get_transferrable('burn_time', serpent, input_parameter=True) # Can be called to control the simulation time (useful for decay)

# Set fuel volumes
for uni_id, (uni_name, uni) in enumerate(active_pebbles_dict.items()):
    V_per_pebble = 4/3 * np.pi * uni['r_fuel_kernel']**3 * uni['Ntrisos']
    N = data[f'pebble_type_{uni_id}'].sum()
    print(f'Calculated volume for {uni["mat_name"]}: {V_per_pebble:.2E} cm^3 (x{N} pebbles)')
    Serpent_set_values(f'material_div_{uni["mat_name"]}_volume', np.full(N, V_per_pebble), serpent) # allows overwriting volume definitions in Serpent

# Try to add step-wise variables, if given
try:
    power_normalization_field = f'norm_{power_normalization_field}'
    tra[power_normalization_field] = get_transferrable(power_normalization_field, serpent, input_parameter=True) # Normalizes the power in Serpent
except:
    pass
try:
    tra['neutrons_per_cycle'] = get_transferrable("neutrons_per_cycle", serpent, input_parameter=True) # Give number of neutrons to Serpent
except:
    pass

#### Serpent => Python #####
tra['keff']          = get_transferrable('ANA_KEFF', serpent) # Monitor multiplication factor
tra['keff_rel_unc']          = get_transferrable('ANA_KEFF_rel_unc', serpent) # Monitor multiplication factor
tra['bu_out'] = {mat_name: get_transferrable(f"material_div_{mat_name}_burnup", serpent) for mat_name in getdict(active_pebbles_dict, 'mat_name')} # Monitor burnup (MWd/kg)
tra['fima_out'] = {mat_name: get_transferrable(f"material_div_{mat_name}_fima", serpent) for mat_name in getdict(active_pebbles_dict, 'mat_name')} # Monitor burnup (%fima)
tra['decayheat'] = {mat_name: get_transferrable(f"material_div_{mat_name}_decayheat", serpent) for mat_name in getdict(active_pebbles_dict, 'mat_name')} # Monitor decay heat, values make no sense (Serpent code to change?), and does not work with DD
tra['activity'] = {mat_name: get_transferrable(f"material_div_{mat_name}_activity", serpent) for mat_name in getdict(active_pebbles_dict, 'mat_name')} # Monitor activity, values make no sense (Serpent code to change?), and does not work with DD

# Monitor isotopes in inventory
for name in inventory_names:
    tra[name] = {mat_name: get_transferrable(f'material_div_{mat_name}_{name}', serpent) for mat_name in getdict(active_pebbles_dict, 'mat_name')}

# Monitor tallies
for name in detector_names:
   tra[name] = get_transferrable(f'DET_{name}', serpent) # Tally value
   tra[f'{name}_rel_unc'] = get_transferrable(f'DET_{name}_rel_unc', serpent) # Tally uncertainty

#%% Create pebble bed object, which will be filled with data at each step

# Import with positions and pebble nature (fuel yes or no)
pbed = Pebble_bed(verbose=False)

pbed.read_dataframe(pd.DataFrame(data[['x','y','z','r']], columns=[*'xyzr'])) # Read first positions
pbed.data['id'] = data['id']

pbed.data['isactive'] = False
for uni_id, (uni_name, uni) in enumerate(pebbles_dict.items()):
    pbed.data[f'pebble_type_{uni_id}'] = (data['uni'] == uni_name)
    if 'mat_name' in uni.keys():
        pbed.data.loc[data['uni'] == uni_name, 'mat_name'] = uni['mat_name']
        pbed.data.loc[data['uni'] == uni_name, 'isactive'] = True

if domain_decomposition:
    pbed.data.loc[pbed.data['isactive'], 'domain_id'] = data['domain_id']

# Initialize columns with default values
pbed.data['initial'] = 1 # Initial pebbles (=1) should not be considered for discarded pebbles, for instance
pbed.data['passes'] = 1 # Monitor number of passes, first all at 1
pbed.data['recirculated'] = False # Monitor which pebbles just recirculated at each step
pbed.data['discarded'] = False # Monitor which pebbles should be discarded at each step
pbed.data['residence_time'] = 0.0 # Monitor total residence time
pbed.data['pass_residence_time'] = 0.0 # Monitor residence time for the current pass
# if not domain_decomposition:
#     pbed.data.loc[pbed.data['isactive'], 'decay_heat'] = 0.0 # Monitor decay heat
pbed.data.loc[pbed.data['isactive'], 'burnup'] = 0.0 # Monitor burnup
pbed.data.loc[pbed.data['isactive'], 'fima'] = 0.0 # Monitor burnup
pbed.data.loc[pbed.data['isactive'], 'pass_burnup'] = 0.0 # Monitor burnup for the current pass
pbed.data['insertion_step'] = 0 # Monitor when pebbles were inserted
pbed.data['avg_r_dist'] = np.array(pbed.data['r_dist']) # Monitor where pebble are radially, on average
pbed.data['pass_agg_r_dist'] = 0.0 # Monitor where pebble are radially, on average
pbed.data['pass_avg_r_dist'] = 0.0 # Monitor where pebble are radially, on average
pbed.data['pass_nsteps'] = 0
pass_dependent_names = ['pass_residence_time', 'pass_burnup', 'pass_agg_r_dist', 'pass_avg_r_dist', 'pass_nsteps']


# Initialize columns for each detector and create time-integrated parameters (fluences, energy)
for name in detector_names:
    pbed.data[name] = np.nan
    pbed.data[f'{name}_rel_unc'] = np.nan
    for subname in [f'integrated_{name}', f'integrated_{name}_unc', f'integrated_{name}_rel_unc', f'pass_integrated_{name}', f'pass_integrated_{name}_unc', f'pass_integrated_{name}_rel_unc']:
        pbed.data[subname] = 0.0
        if 'power' in name: # Only non-fuel should give energy
            pbed.data.loc[~pbed.data['isactive'], subname] = np.nan
    pass_dependent_names = pass_dependent_names + [f'pass_integrated_{name}', f'pass_integrated_{name}_unc', f'pass_integrated_{name}_rel_unc']

if restart_calculation:
    pbed.data[data.columns] = data.copy()
else:
    # Initialize fields corresponding to thresholds
    for uni_id, (uni_name, uni) in enumerate(threshold_pebbles_dict.items()):
        if uni['threshold_type'] == 'passes':
            pbed.data.loc[pbed.data[f'pebble_type_{uni_id}'], 'passes'] = np.random.randint(1, step_wise_variables[f"threshold_{uni['mat_name']}"][0] + 1, pbed.data[f'pebble_type_{uni_id}'].sum())
        else:
            Serpent_set_values(tra[f"{uni['threshold_type']}_{uni['mat_name']}_in"], np.random.uniform(0, step_wise_variables[f"threshold_{uni['mat_name']}"][0], pbed.data[f'pebble_type_{uni_id}'].sum()))

for uni_id, (uni_name, uni) in enumerate(active_pebbles_dict.items()):
    fuel_name = uni['mat_name']

    # Burnup
    pbed.data.loc[pbed.data[f'pebble_type_{uni_id}'], 'burnup'] = Serpent_get_values(tra['bu_out'][fuel_name]).astype(float)
    pbed.data.loc[pbed.data[f'pebble_type_{uni_id}'], 'fima'] = Serpent_get_values(tra['fima_out'][fuel_name]).astype(float)
    pbed.data.loc[pbed.data[f'pebble_type_{uni_id}'], 'decayheat'] = Serpent_get_values(tra['decayheat'][fuel_name]).astype(float)
    pbed.data.loc[pbed.data[f'pebble_type_{uni_id}'], 'activity'] = Serpent_get_values(tra['activity'][fuel_name]).astype(float)

    # Initialize columns for each isotope in inventory
    pbed.data[inventory_names] = np.nan
    for name in inventory_names:
        pbed.data.loc[pbed.data[f'pebble_type_{uni_id}'], name] = Serpent_get_values(tra[name][fuel_name])

#%% Initialize secondary information to get at each step

pbed.cycle_hist = pd.DataFrame(columns=['time', 'passes', 'recirculated', 'discarded', 'keff', 'keff_relative_uncertainty', 'keff_absolute_uncertainty']+list(step_wise_variables.keys())) # keff at each time step
pbed.discarded_data =  pd.DataFrame(columns=list(pbed.data.columns.drop(['x','y','z','r_dist','azim_angle','isactive','recirculated','discarded', 'pass_agg_r_dist']+[f'pebble_type_{i}' for i in range(len(pebbles_dict))]+list(detector_names)+[f'{name}_rel_unc' for name in detector_names]))+pass_dependent_names+['discard_step']) # Discarded pebbles data
if 'uni' in pbed.discarded_data.columns:
    pbed.discarded_data = pbed.discarded_data.drop(columns='uni')
# pbed.discharged_data = pd.DataFrame(columns=list(pbed.data.columns.drop(['x','y','z','r_dist','azim_angle','isactive','recirculated', 'pass_agg_r_dist']+[f'pebble_type_{i}' for i in range(len(pebbles_dict))]+list(detector_names)+[f'{name}_rel_unc' for name in detector_names]))) # Discharged pebbles data

pbed.discharged_fuel_data = pd.DataFrame(columns=list(pbed.data.columns.drop(['x','y','z','r_dist','azim_angle','isactive','recirculated', 'pass_agg_r_dist']+[f'pebble_type_{i}' for i in range(len(pebbles_dict))]+list(detector_names)+[f'{name}_rel_unc' for name in detector_names]))) # Discharged pebbles data, only fuel pebbles

# pbed.reinserted_data = pd.DataFrame(columns=list(pbed.data.columns.drop(['x','y','z','r_dist','azim_angle','isactive','recirculated', 'pass_agg_r_dist']+[f'pebble_type_{i}' for i in range(len(pebbles_dict))]+list(detector_names)+[f'{name}_rel_unc' for name in detector_names]))) # Re-inserted pebbles data
pbed.reinserted_fuel_data = pd.DataFrame(columns=list(pbed.data.columns.drop(['x','y','z','r_dist','azim_angle','isactive','recirculated', 'pass_agg_r_dist']+[f'pebble_type_{i}' for i in range(len(pebbles_dict))]+list(detector_names)+[f'{name}_rel_unc' for name in detector_names]))) # Re-inserted pebbles data, only fuel pebbles

#%% Main loop
for step in range(first_step, Nsteps):
    #### Set up step ####
    pbed.old_data = pbed.data.copy() # Keep old data
    pbed.data['recirculated'] = False # All pebbles are reset to non-recirculated
    pbed.data['discarded'] = False # All pebbles are reset to non-discarded

    #### Set step-wise parameters ####
    if step == 0:
        if "neutrons_per_cycle" in step_wise_variables:
            Serpent_set_values(tra['neutrons_per_cycle'], step_wise_variables["neutrons_per_cycle"][0])
        if "power_normalization_value" in step_wise_variables:
            Serpent_set_values(tra[power_normalization_field], step_wise_variables["power_normalization_value"][0])
    else:
        if discrete_motion:
            motion_step = motion_direction * (np.diff(zrows).mean()*2) * step_wise_variables["Nrows_to_move"][step-1]
        if "neutrons_per_cycle" in step_wise_variables:
            Serpent_set_values(tra['neutrons_per_cycle'], step_wise_variables["neutrons_per_cycle"][step-1])
        if "power_normalization_value":
            Serpent_set_values(tra[power_normalization_field], step_wise_variables["power_normalization_value"][step-1])

    if step == first_step:
        time_step = 0
    else:
        time_step = step_wise_variables["time_step"][step]
    next_time = curr_time + time_step*DAYS
    Npasses = next_time / (time_per_pass*DAYS)
    pbed.cycle_hist.loc[step, 'time'] = next_time/DAYS
    pbed.cycle_hist.loc[step, 'passes'] = Npasses

    print(f'Step {step}: t_ini={curr_time/DAYS:.1f} days, t_end={(next_time/DAYS):.1f} days ({Npasses:.2f} passes)')
    print('\tUsed variables:')
    for var in step_wise_variables:
        if step == 0:
            pbed.cycle_hist.loc[step, var] = step_wise_variables[var][0]
            print(f'\t- {var}: {step_wise_variables[var][0]}')
        else:
            pbed.cycle_hist.loc[step, var] = step_wise_variables[var][step-1]
            print(f'\t- {var}: {step_wise_variables[var][step-1]}')
    print()
    
    #### First step, just run domain decomposition and transport if needed ####
    if step == first_step:
        # Transport step. Not mandatory, only if using no-burnup or specifying to solve first transport
        if transport:
            if resolve_first:
                print('\tSending Serpent signal to run transport')
                serpent.solve()
                print('\tWaiting for Serpent...')
                keff = Serpent_get_values(tra['keff'])[0]
                keff_rel_unc = Serpent_get_values(tra['keff_rel_unc'])[0]
                keff_unc = keff*keff_rel_unc
                pbed.cycle_hist.loc[step, ['keff', 'keff_relative_uncertainty', 'keff_absolute_uncertainty']] = [keff, keff_rel_unc, keff_unc]
                print(f'\tDone. keff = {keff:.5f} +/- {keff_unc*1e5:.0f} pcm')

    #### Other than first step, move pebbles, recirculate, discard, replace with fresh fuel
    else:
        # Move pebbles (discrete motion case)
        if discrete_motion:
            new_positions = pbed.data[['x','y','z']] + np.array([0, 0, motion_step])  # apply motion
            zrows = np.unique(new_positions['z'])

            # Shuffle recirculating, line by line
            for i in range(0, len(zrows[::2])):
                if motion_direction * zrows[2*i] > motion_direction * zlim:
                    beyond_lim = (((new_positions['z']-zrows[2*i]).abs()<1e-3) | ((new_positions['z']-zrows[2*i+1]).abs()<1e-3))
                    new_positions[beyond_lim] = new_positions[beyond_lim].sample(frac=1).values

            # Shift pebbles that should recirculate
            new_positions.loc[motion_direction * new_positions['z'] >= motion_direction * zlim, 'z'] += -motion_direction*(np.abs(zrows[-1]-zrows[0]+(zrows[1]-zrows[0])))
        else:
            # Move pebbles (DEM case), loop if necessary
            nsteps_to_loop = (DEM_end-DEM_start)/DEM_step_increment # won't work if different step increments! after how many steps do we loop
            if step<=nsteps_to_loop: # if first loop (original), no modification
                nloops = 0
                equivalent_step = step
            else: # need to have an equivalent step, and apply transition indices as many times as there were loops
                nloops = int((step-1)//nsteps_to_loop)
                equivalent_step = int((step-1)%nsteps_to_loop)+1
                if nloops!=0 and equivalent_step==1:
                    print(f'Looping! Number of times looped: {nloops}')
            new_positions = pd.read_csv(position_files[equivalent_step])[['x','y','z']]*positions_scale + np.array(positions_translation) # read new positions from DEM file

            # Apply looping transition, if needed
            indices = np.arange(new_positions.shape[0], dtype=int)
            for i in range(nloops): # does not go here if nloops=0
                indices = indices[transition]
            new_positions[['x', 'y', 'z']] = new_positions.loc[indices][['x', 'y', 'z']].values

        new_positions['r'] = r_pebbles # assign r
        new_xyzr = new_positions.values.ravel() # flatten array
        Serpent_set_values(tra['xyzr_in'], new_xyzr) # communicate to Serpent
        pbed.data[[*'xyzr']] = xyzr_to_array(new_xyzr) # fill table with new positions
        pbed.data['r_dist'] = np.linalg.norm(pbed.data[['x', 'y']], axis=1) # calculate new radial positions

        # Decompose in domains and transfer pebbles which change domains
        # if domain_decomposition:
        #     pbed.decompose_in_domains(decomposition_domains, decomposition_types)
        #     changing_domain = (pbed.data[pbed.data['isactive']].domain_id != pbed.old_data[pbed.old_data['isactive']].domain_id)
        #     print(f'To change domain: {changing_domain.sum()}')
        #     Serpent_set_values(tra['domain_in'], pbed.data[pbed.data['isactive']].domain_id) # to correct with multiple fuel

        # Detect recirculated pebbles based on z increment between new and old data
        if discrete_motion:
            recirc_threshold=0
        else:
            # Warning: if the DEM step is huge, many pebbles could recirculate and this would not be valid anymore, then comment that and put the motion_direction in the input
            if (pbed.data['z'] > pbed.old_data['z']).sum()/pbed.data.shape[0] > 0.5: # if the majority went up
                motion_direction = 1
            else:
                motion_direction = -1

        print(f'\tRecirculation criterion: deltaZ {"< -" if motion_direction==+1 else ">"} {recirc_threshold} cm')
        pbed.data['recirculated'] = (-motion_direction*(pbed.data['z'] - pbed.old_data['z']) > recirc_threshold) # select recirculated pebbles (deltaZ > threshold)
        Nrecirculated = pbed.data["recirculated"].sum()
        print(f'\t\t{Nrecirculated} pebbles to recirculate\n')

        # Detect discarded pebbles based on set threshold
        pbed.data['discarded'] = False
        print(f'\tDiscard test for {len(threshold_pebbles_dict)} types of pebbles:')
        for i, (uni_name, uni) in enumerate(threshold_pebbles_dict.items()):
            threshold_type = uni['threshold_type']
            threshold_dir = uni['threshold_dir']
            var = f'threshold_{uni["mat_name"]}'
            threshold_val = step_wise_variables[var][step]
            if threshold_dir > 0:
                print(f'\t\t{uni["mat_name"]} criterion: {threshold_type} >= {threshold_val}')
                to_discard = (pbed.data[f'pebble_type_{i}']) & (pbed.data['recirculated']) & (pbed.data[threshold_type] >= threshold_val) # select recirculated pebbles with the right fuel which satisfy the discard criterion
            else:
                print(f'\t\t{uni["mat_name"]} criterion: {threshold_type} <= {threshold_val}')
                to_discard = (pbed.data[f'pebble_type_{i}']) & (pbed.data['recirculated']) & (pbed.data[threshold_type] <= threshold_val) # select recirculated pebbles with the right fuel which satisfy the discard criterion
            pbed.data.loc[to_discard, 'discarded'] = True

            Ndiscarded = pbed.data.loc[pbed.data[f'pebble_type_{i}'], "discarded"].sum()
            print(f'\t\t\t{Ndiscarded} pebbles to discard\n')
        Ndiscarded_tot = pbed.data["discarded"].sum()
        print(f'\t\t{Ndiscarded} pebbles to discard\n')

        # Add to pebble inventory table
        pbed.cycle_hist.loc[step, ['recirculated', 'discarded']] = [Nrecirculated, Ndiscarded]

        # Before anything, decay pebbles which were recirculated (if decay step>0)
        if transport and decay_step > 0:
            print(f'\tDecaying recirculated pebbles for {decay_step} days')

            # Switch from depletion mode to decay mode (no fission/absorption, just decay)
            Serpent_set_values(tra['switch_mode'], DECAY)

            # Pebbles which are not recirculated should not decay, make them "non-burnable"
            for i, (uni_name, uni) in enumerate(active_pebbles_dict.items()):
                not_decaying = ~pbed.data.loc[pbed.data[f'pebble_type_{i}'], "recirculated"].values # select non-recirculated pebbles
                Serpent_set_multiple_values(tra['burnable_fuel'][getdict(active_pebbles_dict, 'mat_name')[i]][not_decaying], np.zeros(not_decaying.sum(), dtype=int)) # make them non-burnable

            # Decay burnable pebbles
            serpent.advance_to_time(curr_time + decay_step*DAYS)

            # Come back to normal
            for i, (uni_name, uni) in enumerate(active_pebbles_dict.items()):
                not_decaying = ~pbed.data.loc[pbed.data[f'pebble_type_{i}'], "recirculated"].values
                Serpent_set_multiple_values(tra['burnable_fuel'][getdict(active_pebbles_dict, 'mat_name')[i]][not_decaying], np.ones(not_decaying.sum(), dtype=int))

            Serpent_set_values(tra['switch_mode'], DEPLETION) # back to depletion mode
            Serpent_set_values(tra['time_in'], curr_time) # we artificially applied decay for a given time, but time should not change: come back to time before decay

        # Store discharged information
        print(f'\tStoring {Nrecirculated} discharged pebbles information')
        # pbed.discharged_data = pbed.data.loc[pbed.data['recirculated'], pbed.discharged_data.columns].copy()
        pbed.discharged_fuel_data = pbed.data.loc[((pbed.data['isactive']) & pbed.data['recirculated']), pbed.discharged_fuel_data.columns].copy()

        # Store discarded information
        print(f'\tStoring {Ndiscarded} discarded pebbles information')
        pbed.discarded_data = pbed.data.loc[pbed.data['discarded'], pbed.discarded_data.columns.drop('discard_step')].copy() # Copy into discarded data (only non-initial are considered)
        pbed.discarded_data['discard_step'] = step-1 # Add discarded step

        # Insert fresh pebbles
        print(f'\tInserting {Ndiscarded} fresh pebbles')
        if transport:
            for uni_id, (uni_name, uni) in enumerate(threshold_pebbles_dict.items()):
                material = uni['mat_name']
                Serpent_set_multiple_values(tra['reset_fuel'][material][pbed.data.loc[pbed.data[f"pebble_type_{uni_id}"], "discarded"]], np.ones(Ndiscarded).astype(int)) # Reset fuel composition for discarded pebbles (=>fresh)
                Serpent_set_multiple_values(tra['burnup_in'][material][pbed.data.loc[pbed.data[f"pebble_type_{uni_id}"], "discarded"]], np.zeros(Ndiscarded)) # Reset burnup for discarded pebbles to 0

            # Reset integrated tallies to 0
            for name in detector_names:
                pbed.data.loc[pbed.data['discarded'], f'integrated_{name}'] = 0.0
                pbed.data.loc[pbed.data['discarded'], f'integrated_{name}_unc'] = 0.0

        pbed.data.loc[pbed.data['discarded'], 'insertion_step'] = step # set insertion step of fresh pebble to current step
        pbed.data.loc[pbed.data['discarded'], 'residence_time'] = 0.0 # reset residence time to 0
        pbed.data.loc[pbed.data['discarded'], 'passes'] = 1 # reset number of passes to 1
        pbed.data.loc[pbed.data['discarded'], 'initial'] = 0 # once replaced at least once, pebbles are not initial anymore

        # Store reinserted information (fresh + non-discarded pebbles)
        print(f'\tStoring {Nrecirculated} reinserted pebbles information')
        # pbed.reinserted_data = pbed.data.loc[pbed.data['recirculated'], pbed.reinserted_data.columns].copy()
        pbed.reinserted_fuel_data = pbed.data.loc[(pbed.data['isactive']) & (pbed.data['recirculated']), pbed.reinserted_fuel_data.columns].copy()
        
        for name in inventory_names:
            for uni_id, (uni_name, uni) in enumerate(threshold_pebbles_dict.items()):
                id_nuc = nuclides_list[uni['mat_name']].index(int(name))
                discarded = (pbed.data['discarded']) & (pbed.data[f'pebble_type_{uni_id}'])
                if discarded.sum() > 0:
                    pbed.reinserted_fuel_data.loc[discarded, name] = fresh_fuel[uni['mat_name']][id_nuc]

        # Reset pass-dependent variables for recirculated data
        for name in pass_dependent_names:
            pbed.data.loc[pbed.data['recirculated'], name] = 0.0
            if 'power' in name or 'burnup' in name:
                pbed.data.loc[(pbed.data['recirculated']) & (~pbed.data['isactive']), name] = np.nan

        curr_time = float(next_time) # Increment time
        pbed.data['pass_nsteps'] += 1
        pbed.data['pass_residence_time'] += time_step  # Increment pass residence times
        pbed.data['residence_time'] += time_step  # Increment residence times
        pbed.data.loc[(pbed.data['recirculated']) & (~pbed.data['discarded']), 'passes'] += 1 # Rercirculating (not fresh) pebbles have an incremented pass

        # Run Serpent for transport and depletion
        if transport:
            print('\tSending Serpent signal to run transport')
            serpent.advance_to_time(next_time) # Run Serpent transport + burn to next time (predictor if using predictor/corrector)
            if correct:
                serpent.correct() # Run corrector step if using (predictor/corrector)
            if write_restart and step%restart_write_every==0:
                print(f'\tWriting restart file for step {step}')
                Serpent_set_values(tra['write_restart'], step)
            print('\tWaiting for Serpent...')
            keff = Serpent_get_values(tra['keff'])[0]
            keff_rel_unc = Serpent_get_values(tra['keff_rel_unc'])[0]
            keff_unc = keff*keff_rel_unc
            pbed.cycle_hist.loc[step, ['keff', 'keff_relative_uncertainty', 'keff_absolute_uncertainty']] = [keff, keff_rel_unc, keff_unc]
            print(f'\tDone. keff = {keff:.5f} +/- {keff_unc*1e5:.0f} pcm')

    #### Monitor quantities ####
    print('\tExtracting burnups, detectors and compositions')

    if transport:
        # if not domain_decomposition:
            # # Extract decay heat
            # pbed.data.loc[pbed.data['isactive'], 'decay_heat'] = Serpent_get_values(tra['decay_heat']).astype(float)

        # Extract burnups and calculate pass burnup
        for uni_id, (uni_name, uni) in enumerate(active_pebbles_dict.items()):
            material = uni['mat_name']
            pbed.data.loc[pbed.data[f'pebble_type_{uni_id}'], 'burnup'] = Serpent_get_values(tra['bu_out'][material]).astype(float)
            pbed.data.loc[pbed.data[f'pebble_type_{uni_id}'], 'fima'] = Serpent_get_values(tra['fima_out'][material]).astype(float)
        if step > 0:
            # Pass burnup incremented for fuel pebbles
            pbed.data.loc[(pbed.data['isactive'] & (~pbed.data['recirculated'])), 'pass_burnup'] += pbed.data.loc[pbed.data['isactive'], 'burnup'] - pbed.old_data.loc[pbed.data['isactive'], 'burnup']

            # Special case for pebbles which recirculated
            pbed.data.loc[(pbed.data['isactive'] & (pbed.data['recirculated'])), 'pass_burnup'] = pbed.data.loc[(pbed.data['isactive'] & (pbed.data['recirculated'])), 'burnup'] - pbed.old_data.loc[(pbed.data['isactive'] & (pbed.data['recirculated'])), 'burnup']

        # Extract isotopic inventory
        for uni_id, (uni_name, uni) in enumerate(active_pebbles_dict.items()):
            material = uni['mat_name']
            for name in inventory_names:
                pbed.data.loc[pbed.data[f'pebble_type_{uni_id}'], name] = Serpent_get_values(tra[name][material])

        # Extract detector values and uncertainties and calculate integrated values
        for name in list(detector_names):
            pbed.data.to_csv('test.csv')
            # Tallies
            pbed.data[name] = Serpent_get_values(tra[name])
            pbed.data[f'{name}_rel_unc'] =  Serpent_get_values(tra[f'{name}_rel_unc'])

            # Integrated tallies
            if step > 0:
                added = time_step * DAYS * pbed.data[name]
                added_unc = added * pbed.data[f'{name}_rel_unc']
                pbed.data[f'integrated_{name}'] += added
                pbed.data[f'integrated_{name}_unc'] += added_unc
                pbed.data[f'integrated_{name}_rel_unc'] = pbed.data[f'integrated_{name}_unc']/pbed.data[f'integrated_{name}'] # to modify
                pbed.data[f'pass_integrated_{name}'] += added
                pbed.data[f'pass_integrated_{name}_unc'] += added_unc
                pbed.data[f'pass_integrated_{name}_rel_unc'] = pbed.data[f'pass_integrated_{name}_unc']/pbed.data[f'pass_integrated_{name}']  # to modify
                if 'power' in name: # Only fuel has power, rest is "nan"
                    for subname in [f'integrated_{name}', f'integrated_{name}_unc', f'integrated_{name}_rel_unc', f'pass_integrated_{name}', f'pass_integrated_{name}_unc', f'pass_integrated_{name}_rel_unc']:
                        pbed.data.loc[~pbed.data['isactive'], subname] = np.nan

    # Calculate new average distance
    pbed.data['avg_r_dist'] = (pbed.data['avg_r_dist']*(step-pbed.data['insertion_step'])+pbed.data['r_dist'])/(step-pbed.data['insertion_step']+1)
    pbed.data['pass_agg_r_dist'] += pbed.data['r_dist']*time_step
    pbed.data['pass_avg_r_dist'] = pbed.data['pass_agg_r_dist']/pbed.data['pass_residence_time']

    #### Plots ###
    if plotting:
        print('\tPlotting')

        # Plot latest geometry and save it to folder
        plot_Serpent(tra['plot'], plots_folder_path='wrk_Serpent', output_suffix=f'_{step}', output_folder='Plots', nplots=nplots, delay=60)

        # Create copy of table and fill non-fuel data with "nans"
        pbed_fuel = deepcopy(pbed)
        pbed_fuel.data.loc[~pbed.data['isactive'], pbed.data.columns.drop(['x', 'y', 'z', 'r', 'isactive'] + [f'pebble_type_{i}' for i in range(len(pebbles_dict))])] = np.nan
        cmap = cm.get_cmap('turbo')
        cmap.set_bad([0.8, 0.8, 0.8], alpha = 1.) # nans will appear gray

        # Slices of the core, with fuel-only or whole core
        plt.close('all')
        if transport:
            plt.figure(figsize=(11.5, 15))
            plt.subplot(3,4,1)
            pbed.plot2D('id', field_title='Pebble ID', plot_title=f'Step {step}: {curr_time/DAYS:.1f} days', colormap=cmap, new_fig=False, shrink_cbar=0.9, pad_cbar=0.12, clim=[0, pbed.data['id'].max()], verbose=False)
            plt.subplot(3,4,2)
            pbed.plot2D('pass_residence_time', field_title='Pass residence time [days]', plot_title='', colormap=cmap, new_fig=False, shrink_cbar=0.9, pad_cbar=0.12, clim=[0, pbed.data['pass_residence_time'].max()], verbose=False)
            plt.subplot(3,4,3)
            pbed_fuel.plot2D('residence_time', field_title='Total residence time [days]', plot_title='', colormap=cmap, new_fig=False, shrink_cbar=0.9, pad_cbar=0.12, clim=[0, pbed_fuel.data['residence_time'].max()], verbose=False)
            plt.subplot(3,4,4)
            pbed_fuel.plot2D('passes', field_title='Number of passes', plot_title='',  colormap=cmap, new_fig=False, shrink_cbar=0.9, pad_cbar=0.12, clim=[1, pbed_fuel.data['passes'].max()], verbose=False)
            plt.subplot(3,4,5)
            pbed.plot2D('flux_pebbles_thermal', field_title='Thermal flux [n/cm$^2$.s]', plot_title='',  colormap=cmap, new_fig=False, shrink_cbar=0.9, pad_cbar=0.12, clim=[0, pbed.data['flux_pebbles_thermal'].max()], verbose=False)
            plt.subplot(3,4,6)
            pbed.plot2D('flux_pebbles_fast', field_title='Fast flux [n/cm$^2$.s]', plot_title='',  colormap=cmap, new_fig=False, shrink_cbar=0.9, pad_cbar=0.12, clim=[0, pbed.data['flux_pebbles_fast'].max()], verbose=False)
            plt.subplot(3,4,7)
            pbed_fuel.plot2D('power_pebbles', field_title='Power [W]', plot_title='',  colormap=cmap, new_fig=False, shrink_cbar=0.9, pad_cbar=0.12, clim=[0, pbed_fuel.data['power_pebbles'].max()], verbose=False)
            plt.subplot(3,4,8)
            pbed_fuel.plot2D('burnup', field_title='Burnup [MWd/kg]', plot_title='',  colormap=cmap, new_fig=False, shrink_cbar=0.9, pad_cbar=0.12, clim=[0, pbed_fuel.data['burnup'].max()], verbose=False)
            plt.subplot(3,4,9)
            pbed_fuel.plot2D('integrated_flux_pebbles_thermal', field_title='Thermal fluence [n/cm$^2$]', plot_title='',  colormap=cmap, new_fig=False, shrink_cbar=0.9, pad_cbar=0.12, clim=[0, pbed_fuel.data['integrated_flux_pebbles_thermal'].max()], verbose=False)
            plt.subplot(3,4,10)
            pbed_fuel.plot2D('integrated_flux_pebbles_fast', field_title='Fast fluence [n/cm$^2$]', plot_title='',  colormap=cmap, new_fig=False, shrink_cbar=0.9, pad_cbar=0.12, clim=[0, pbed_fuel.data['integrated_flux_pebbles_fast'].max()], verbose=False)
            plt.subplot(3,4,11)
            pbed_fuel.plot2D('integrated_power_pebbles', field_title='Energy [Wd]', plot_title='',  colormap=cmap, new_fig=False, shrink_cbar=0.9, pad_cbar=0.12, clim=[0, pbed_fuel.data['integrated_power_pebbles'].max()], verbose=False)
            plt.subplot(3,4,12)
            pbed_fuel.plot2D('pass_burnup', field_title='$\Delta$BU [MWd/kg]', plot_title='',  colormap=cmap, new_fig=False, shrink_cbar=0.9, pad_cbar=0.12, clim=[0, pbed_fuel.data['pass_burnup'].max()], verbose=False)
        else:
            plt.figure(figsize=(11.5, 15))
            plt.subplot(1,4,1)
            pbed.plot2D('id', field_title='Pebble ID', plot_title=f'Step {step}: {curr_time/DAYS:.1f} days', colormap=cmap, new_fig=False, shrink_cbar=0.9, pad_cbar=0.12, clim=[0, pbed.data['id'].max()], verbose=False)
            plt.subplot(1,4,2)
            pbed.plot2D('pass_residence_time', field_title='Pass residence time [days]', plot_title='', colormap=cmap, new_fig=False, shrink_cbar=0.9, pad_cbar=0.12, clim=[0, pbed.data['pass_residence_time'].max()], verbose=False)
            plt.subplot(1,4,3)
            pbed.plot2D('residence_time', field_title='Total residence time [days]', plot_title='', colormap=cmap, new_fig=False, shrink_cbar=0.9, pad_cbar=0.12, clim=[0, pbed.data['residence_time'].max()], verbose=False)
            plt.subplot(1,4,4)
            pbed_fuel.plot2D('passes', field_title='Number of passes', plot_title='',  colormap=cmap, new_fig=False, shrink_cbar=0.9, pad_cbar=0.12, clim=[1, pbed.data['passes'].max()], verbose=False)
        plt.tight_layout()
        plt.savefig(f'Plots/core_{step}.png', dpi=400, bbox_inches='tight')# save to output

        # # Plot discarded pebbles
        # if len(pbed.discarded_data)>0:
        #     plt.figure(figsize=(11.5, 15))
        #     pbed.discarded_data.drop(columns=[i for i in pbed.discarded_data if '_unc' in i]).astype(float).hist() # Plot distributions for each quantity
        #     plt.tight_layout()
        #     plt.savefig(f'Plots/discarded_{step}.png', dpi=400, bbox_inches='tight')

        # Plot keff
        if transport:
            plt.figure()
            errorbar(pbed.cycle_hist['passes'].astype(float), pbed.cycle_hist['keff'].astype(float), pbed.cycle_hist['keff_absolute_uncertainty'].astype(float), label='k$_{eff}$', labelize_error='Uncertainty', linewidth=1.5)
            plt.xlabel('Passes')
            plt.ylabel('Multiplication factor')
            plt.grid()
            plt.gca().set_axisbelow(True)
            plt.savefig(f'Plots/keff_{step}.png', dpi=400, bbox_inches='tight')

        plt.show()

    #### Save data ####
    if saving:
        if write_global:
            pbed.cycle_hist.to_csv(f'Data/cycle_{step}.csv') # global inventory data
        if write_incore:
            pbed.data.to_csv(f'Data/core_{step}.csv') # whole core data
        if write_discharged:
            # pbed.discharged_data.to_csv(f'Data/discharged_{step}.csv') # discarded pebbles data
            pbed.discharged_fuel_data.to_csv(f'Data/discharged_fuel_{step}.csv') # discarded pebbles data
        if write_discarded:
            pbed.discarded_data.to_csv(f'Data/discarded_{step}.csv') # discarded pebbles data
        if write_reinserted:
            # pbed.reinserted_data.to_csv(f'Data/reinserted_{step}.csv') # re-inserted pebbles data
            pbed.reinserted_fuel_data.to_csv(f'Data/reinserted_fuel_{step}.csv') # re-inserted pebbles data, only fuel pebbles
# %%
