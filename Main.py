print("""
============================================================
||                                                        ||
||    8 8888        8 `8.`8888.      ,8' 8 88888888888    ||
||    8 8888        8  `8.`8888.    ,8'  8 8888           ||
||    8 8888        8   `8.`8888.  ,8'   8 8888           ||
||    8 8888        8    `8.`8888.,8'    8 8888           ||
||    8 8888        8     `8.`88888'     8 88888888888    ||
||    8 8888        8     .88.`8888.     8 8888           ||
||    8 8888888888888    .8'`8.`8888.    8 8888           ||
||    8 8888        8   .8'  `8.`8888.   8 8888           ||
||    8 8888        8  .8'    `8.`8888.  8 8888           ||
||    8 8888        8 .8'      `8.`8888. 8 8888           ||
||                                                        ||
============================================================
""")


#%% Import pebble bed library and utilities

from Pebble_Bed import *
from Utilities import *
from Looping import *

import pandas as pd
import os
from glob import glob
import numpy as np
import importlib
import atexit
import sys
from copy import deepcopy
from IPython import get_ipython

#%% Constants

DEPLETION = 1
DECAY = 2
DAYS = 86400

#%% Read nodes and cores from arguments and set domain decomposition if needed (need for explicit name in Jupyter)
if get_ipython().__class__.__name__ != 'ZMQInteractiveShell':
    filename_input = sys.argv[1]
else:
    filename_input = 'Input'

#%% Read input from file
globals().update(importlib.import_module('./Default_Input.py'.replace(".py", "").replace('./', '')).__dict__) # Will read all defaults parameters in input file and import them here
globals().update(importlib.import_module(filename_input.replace(".py", "").replace('./', '')).__dict__) # Will read all parameters in input file and import them here (overwrites default)

try:
    ncores = int(sys.argv[2])
    if len(sys.argv)>2:
        nnodes = int(sys.argv[3])
except:
    pass

#%% Domain decomposition
if domain_decomposition:
    decomposition_types, decomposition_domains = nodes_to_dd(nnodes, allowed_decomposition_types, max_domains)
    print(f'Using Domain Decomposition with domains of types "{decomposition_types}" and breakout {decomposition_domains}')

#%% Prepare the case

os.chdir(os.path.dirname(__file__)) # Get to the current path
original_path = str(os.getcwd())
init_case(case_name, filename_input, path_to_case, output_folder_name) # Copy input files in output folder
os.chdir(f'Cases/{output_folder_name}') # Go to output folder
atexit.register(os.chdir, original_path) # Come back to original folder when the script ends (even if error)

#%% Change positions, universes and radii based on input

# If restarting from a step, calculate equivalent step and positions
if restart_calculation:
    first_step = restart_step
    print(f'Restarting from step {restart_step}, binary data at "{restart_binary}" and reading data at "{restart_data}".')
    if not os.path.exists(restart_data):
        raise Exception(f'Restart mode selected but no restart data table found at {restart_data}')
    if not os.path.exists(restart_binary):
        raise Exception(f'Restart mode selected but no restart binary found at {restart_binary}')
    else:
        nrestarts = len(glob(f'{restart_binary}*')) # count number of restart files
    with open(main_input_file, 'a') as f:
        f.write(f'\nset rfr idx 0 "{restart_binary}" {nrestarts}\n')
else:
    first_step = 0
    if read_firt_compositions:
        if not os.path.exists(restart_data):
            raise Exception(f'First composition binary mode selected but no restart data table found at {restart_data}')
        if not os.path.exists(restart_binary):
            raise Exception(f'First composition binary mode selected but no restart binary found at {restart_binary}')
        else:
            nrestarts = len(glob(f'{restart_binary}*')) # count number of restart files
        with open(main_input_file, 'a') as f:
            f.write(f'\nset rfr idx 0 "{restart_binary}" {nrestarts}\n')

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
    print(f'Assigning random universes with fuel fraction of {fuel_frac*100:.2f}%\n')
    if fuel_frac != 1:
        data['uni'] = assign_random_array(np.arange(len(data)), [fuel_pebble_universe_name, graph_pebble_universe_name], [fuel_frac, 1-fuel_frac]) # Spl
    else:
        data['uni'] = fuel_pebble_universe_name
else:
    data.loc[data['fuel'], 'uni']  = fuel_pebble_universe_name
    data.loc[~data['fuel'], 'uni'] = graph_pebble_universe_name
# Save
data[['x','y','z','r','uni']].to_csv(pbed_file, header=False, index=False, sep='\t')
data.loc[list(data.sort_values(['z', 'r_dist']).index), 'id'] = np.arange(len(data)).astype(int) # sort by z and radial distance

# Count
Npebbles = len(data['r'])
Nfuel = sum(['fuel' in u for u in data['uni']])

# Prepare case
serpent, serpent_input_files = start_Serpent(sssexe, ncores, main_input_file, nnodes, verbosity_mode)
nplots = count_plots(serpent_input_files) # count number of plot commands in input files
#erase_working_plots() # erase the first plots

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
fuel_indices = np.array([int(name.split(f'{fuel_mat_name}z')[-1]) for name in natural_sort(Serpent_get_values('materials', serpent)) if f'{fuel_mat_name}z' in name]) # Indices for fuel pebbles (=row numbers in pbed file)
nuclides_list = list(Serpent_get_values(f'composition_{fuel_mat_name}z{fuel_indices[0]}_zai', serpent)) # List of nuclides in fuel
fresh_fuel = Serpent_get_values(f'composition_{fuel_mat_name}z{fuel_indices[0]}_adens', serpent) # Fresh fuel composition
if domain_decomposition:
    initial_domains = Serpent_get_values(f"material_div_{fuel_mat_name}_domain", serpent) # Initial domains for each pebble

#%% Create transferrables, which will be callable during calculation (Serpent <=> Python)

# Transferrables will be in the "tra" dictionnary
tra = dict()

##### Python => Serpent #####
tra['plot']          = Serpent_set_values("plot_geometry", 1, serpent, communicate=False) # Can be called to plot latest position/geometry
tra['write_restart'] = get_transferrable("write_restart", serpent, input_parameter=True) # Can be called to plot latest position/geometry
tra['xyzr_in']       = get_transferrable(f"pbed_{pbed_universe_name}_xyzr", serpent, input_parameter=True) # Can be called to change pebbles positions
tra['bu_in']         = Serpent_get_material_wise(fuel_mat_name, 'burnup', serpent, input_parameter=True) # Can be called to change pebble-wise burnup
tra['reset_fuel']    = Serpent_get_material_wise(fuel_mat_name, 'reset', serpent, prefix='composition', input_parameter=True) # Can be called to reset the fuel composition to the original one (discarded -> fresh)
#tra['compo_in']    = Serpent_get_material_wise(fuel_mat_name, 'adens', serpent, prefix='composition', input_parameter=True) # Can be called to reset the fuel composition to the original one (discarded -> fresh)
tra['burnable_fuel'] = Serpent_get_material_wise(fuel_mat_name, 'burnable', serpent, prefix='material', input_parameter=True) # Can be called to switch between burnable and non-burnable pebbles (useful for decay)
tra['switch_mode']   = get_transferrable('burn_step_type', serpent, input_parameter=True) # Can be called to switch between decay mode and depletion mode (useful for decay)
tra['time_in']       = get_transferrable('burn_time', serpent, input_parameter=True) # Can be called to control the simulation time (useful for decay)
if domain_decomposition:
    tra['domain_in']   = get_transferrable(f"material_div_{fuel_mat_name}_domain", serpent, input_parameter=True) # Can be called to control which domain pebbles are transfered to
if threshold_type not in ['passes', 'burnup']:
    tra[f'{threshold_type}_in'] = Serpent_get_material_wise(fuel_mat_name, threshold_type, serpent, prefix='composition', input_parameter=True)

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
tra['bu_out']      = get_transferrable(f"material_div_{fuel_mat_name}_burnup", serpent) # Monitor burnup
tra['keff']          = get_transferrable('ANA_KEFF', serpent) # Monitor multiplication factor
tra['keff_rel_unc']          = get_transferrable('ANA_KEFF_rel_unc', serpent) # Monitor multiplication factor

# Monitor isotopes in inventory
inventory_names  = list(Serpent_get_values('inventory', serpent))
for name in inventory_names:
    tra[name] = get_transferrable(f'material_div_{fuel_mat_name}_{name}', serpent)

# Monitor tallies
for name in detector_names:
   tra[name] = get_transferrable(f'DET_{name}', serpent) # Tally value
   tra[f'{name}_rel_unc'] = get_transferrable(f'DET_{name}_rel_unc', serpent) # Tally uncertainty

#%% Create pebble bed object, which will be filled with data at each step

# Import with positions and pebble nature (fuel yes or no)
pbed = Pebble_bed(verbose=False)

pbed.read_dataframe(pd.DataFrame(data[['x','y','z','r']], columns=[*'xyzr'])) # Read first positions
pbed.data['id'] = data['id']
pbed.data['fuel'] = ['fuel' in u for u in data['uni']]

# Initialize columns with default values
pbed.data['initial'] = 1 # Initial pebbles (=1) should not be considered for discarded pebbles, for instance
pbed.data['passes'] = 1 # Monitor number of passes, first all at 1
pbed.data['recirculated'] = False # Monitor which pebbles just recirculated at each step
pbed.data['discarded'] = False # Monitor which pebbles should be discarded at each step
pbed.data['residence_time'] = 0.0 # Monitor total residence time
pbed.data['pass_residence_time'] = 0.0 # Monitor residence time for the current pass
pbed.data.loc[pbed.data['fuel'], 'burnup'] = 0 # Monitor burnup
pbed.data.loc[pbed.data['fuel'], 'pass_burnup'] = 0.0 # Monitor burnup for the current pass
pbed.data['insertion_step'] = 0 # Monitor when pebbles were inserted
pbed.data['avg_r_dist'] = np.array(pbed.data['r_dist']) # Monitor where pebble are radially, on average
pbed.data['agg_r_dist_pass'] = 0.0 # Monitor where pebble are radially, on average

# Initialize columns for each detector and create time-integrated parameters (fluences, energy)
for name in detector_names:
    pbed.data[name] = np.nan
    pbed.data[f'{name}_rel_unc'] = np.nan
    pbed.data[f'integrated_{name}'] = 0.0
    pbed.data[f'integrated_{name}_unc'] = 0.0
    pbed.data[f'integrated_{name}_rel_unc'] = 0.0
    if 'power' in name: # Only non-fuel should give energy
        pbed.data.loc[~pbed.data['fuel'], f'integrated_{name}'] = np.nan
        pbed.data.loc[~pbed.data['fuel'], f'integrated_{name}_unc'] = np.nan
        pbed.data.loc[~pbed.data['fuel'], f'integrated_{name}_rel_unc'] = np.nan

if restart_calculation:
    pbed.data[data.columns] = data.copy()
else:
    # Initialize field corresponding to threshold
    if threshold_type == 'passes':
        pbed.data.loc[pbed.data['fuel'], threshold_type] = np.random.randint(1, threshold+1, Nfuel) # Uniform integer distribution between 0 and threshold value, in Python
    else:
        Serpent_set_multiple_values(tra[f'{threshold_type}_in'], np.random.uniform(0, threshold, Nfuel)) # Uniform distribution between 0 and threshold value, in Serpent

#%% Communicate some data from Serpent (checking purpose, threshold, and filling data for restart)

# Burnup
pbed.data.loc[pbed.data['fuel'], 'burnup'] = Serpent_get_values(tra['bu_out'])

# Initialize columns for each isotope in inventory
for name in inventory_names:
    pbed.data.loc[pbed.data['fuel'], name] = Serpent_get_values(tra[name])

# Domain decomposition
if domain_decomposition:
    pbed.data.loc[pbed.data['fuel'], 'domain_id'] = initial_domains # Monitor which domain pebbles are in

#%% Initialize secondary information to get at each step

pbed.cycle_hist = pd.DataFrame(columns=['time', 'passes', 'recirculated', 'discarded', 'keff', 'keff_relative_uncertainty', 'keff_absolute_uncertainty']) # keff at each time step
pbed.discarded_data =  pd.DataFrame(columns=list(pbed.data.columns.drop(['x','y','z','r_dist','azim_angle','fuel','recirculated','discarded','pass_residence_time','pass_burnup']+list(detector_names)+[f'{name}_rel_unc' for name in detector_names]))+['discard_step']) # Discarded pebbles data
if 'uni' in pbed.discarded_data.columns:
    pbed.discarded_data = pbed.discarded_data.drop(columns='uni')
pbed.reinserted_data = pd.DataFrame(columns=list(pbed.data.columns.drop(['x','y','z','r_dist','azim_angle','fuel','recirculated','discarded']+list(detector_names)+[f'{name}_rel_unc' for name in detector_names]))) # Re-inserted pebbles data
pbed.reinserted_fuel_data = pd.DataFrame(columns=list(pbed.data.columns.drop(['x','y','z','r_dist','azim_angle','fuel','recirculated','discarded']+list(detector_names)+[f'{name}_rel_unc' for name in detector_names]))) # Re-inserted pebbles data, only fuel pebbles
pbed.tracked_reinserted_nuclide = pd.DataFrame() # added for ML project

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
            print(f'\t- {var}: {step_wise_variables[var][0]}')
        else:
            print(f'\t- {var}: {step_wise_variables[var][step-1]}')
    print()

    #### First step, just run domain decomposition and transport if needed ####
    if step == first_step:
        # Decompose in domains and transfer pebbles which change domains
        if domain_decomposition:
            pbed.decompose_in_domains(decomposition_domains, decomposition_types)
            changing_domain = (pbed.data[pbed.data['fuel']].domain_id != initial_domains)
            print(f'\tTo change domain: {changing_domain.sum()}')
            Serpent_set_values(tra['domain_in'], pbed.data[pbed.data['fuel']].domain_id)

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
            for i in range(len(zrows[::2])):
                if motion_direction * zrows[i] > motion_direction * zlim:
                    new_positions.loc[(new_positions['z']-zrows[i]).abs()<1e-3] = new_positions.loc[(new_positions['z']-zrows[i]).abs()<1e-3].sample(frac=1).values
                if motion_direction * zrows[i+1] > motion_direction * zlim:
                    new_positions.loc[(new_positions['z']-zrows[i+1]).abs()<1e-3] = new_positions.loc[(new_positions['z']-zrows[i+1]).abs()<1e-3].sample(frac=1).values

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
        if domain_decomposition:
            pbed.decompose_in_domains(decomposition_domains, decomposition_types)
            changing_domain = (pbed.data[pbed.data['fuel']].domain_id != pbed.old_data[pbed.old_data['fuel']].domain_id)
            print(f'To change domain: {changing_domain.sum()}')
            Serpent_set_values(tra['domain_in'], pbed.data[pbed.data['fuel']].domain_id)

        # Detect recirculated pebbles based on z increment between new and old data
        print(f'\tRecirculation criterion: deltaZ >= {recirc_threshold} cm')
        pbed.data['recirculated'] = (-motion_direction*(pbed.data['z'] - pbed.old_data['z']) >= recirc_threshold) # select recirculated pebbles (deltaZ > threshold)
        Nrecirculated = pbed.data["recirculated"].sum()
        print(f'\t\t{Nrecirculated} pebbles to recirculate\n')

        # Detect discarded pebbles based on set threshold
        print(f'\tDiscard criterion: {threshold_type} > {threshold}')
        pbed.data.loc[pbed.data['recirculated'] & (pbed.data[threshold_type] > threshold), 'discarded'] = True
        Ndiscarded = pbed.data["discarded"].sum()
        print(f'\t\t{Ndiscarded} pebbles to discard\n')

        # Add to pebble inventory table
        pbed.cycle_hist.loc[step, ['recirculated', 'discarded']] = [Nrecirculated, Ndiscarded]

        # Before anything, decay pebbles which were recirculated (if decay step>0)
        if transport and decay_step > 0:
            print(f'\tDecaying recirculated pebbles for {decay_step} days')

            # Switch from depletion mode to decay mode (no fission/absorption, just decay)
            Serpent_set_values(tra['switch_mode'], DECAY)

            # Pebbles which are not recirculated should not decay, make them "non-burnable"
            not_decaying = ~pbed.data.loc[pbed.data['fuel'], "recirculated"].values # select non-recirculated pebbles
            Serpent_set_multiple_values(tra['burnable_fuel'][not_decaying], np.zeros(not_decaying.sum()).astype(int)) # make them non-burnable

            # Decay burnable pebbles
            serpent.advance_to_time(curr_time + decay_step*DAYS)

            # Come back to normal
            Serpent_set_multiple_values(tra['burnable_fuel'][not_decaying], np.ones(not_decaying.sum()).astype(int)) # all pebbles are burnable
            Serpent_set_values(tra['switch_mode'], DEPLETION) # back to depletion mode
            Serpent_set_values(tra['time_in'], curr_time) # we artificially applied decay for a given time, but time should not change: come back to time before decay

        # Store discarded information
        print(f'\tStoring {Ndiscarded} discarded pebbles information')
        pbed.discarded_data = pbed.data.loc[pbed.data['discarded'], pbed.discarded_data.columns.drop('discard_step')].copy() # Copy into discarded data (only non-initial are considered)
        pbed.discarded_data['discard_step'] = step-1 # Add discarded step

        # Insert fresh pebbles
        print(f'\tInserting {Ndiscarded} fresh pebbles')
        if transport:
            Serpent_set_multiple_values(tra['reset_fuel'][pbed.data.loc[pbed.data["fuel"],"discarded"].values], np.ones(Ndiscarded).astype(int)) # Reset fuel composition for discarded pebbles (=>fresh)
            Serpent_set_multiple_values(tra['bu_in'][pbed.data.loc[pbed.data["fuel"],"discarded"].values], np.zeros(Ndiscarded)) # Reset burnup for discarded pebbles to 0

            # Reset integrated tallies to 0
            for name in detector_names:
                pbed.data.loc[pbed.data['discarded'], f'integrated_{name}'] = 0.0
                pbed.data.loc[pbed.data['discarded'], f'integrated_{name}_unc'] = 0.0

        pbed.data.loc[pbed.data['discarded'], 'insertion_step'] = step # set insertion step of fresh pebble to current step
        pbed.data.loc[pbed.data['discarded'], 'residence_time'] = 0 # reset residence time to 0
        pbed.data.loc[pbed.data['discarded'], 'passes'] = 1 # reset number of passes to 1
        pbed.data.loc[pbed.data['discarded'], 'initial'] = 0 # once replaced at least once, pebbles are not initial anymore

        # Store reinserted information (fresh + non-discarded pebbles)
        print(f'\tStoring {Nrecirculated} reinserted pebbles information')
        pbed.reinserted_data = pbed.data.loc[pbed.data['recirculated'], pbed.reinserted_data.columns].copy()
        pbed.reinserted_fuel_data = pbed.data.loc[((pbed.data['fuel']) & pbed.data['recirculated']), pbed.reinserted_data.columns].copy()
        for i in inventory_names:
            pbed.reinserted_data.loc[pbed.data['discarded'], i] = fresh_fuel[nuclides_list.index(int(i))]
            pbed.reinserted_fuel_data.loc[pbed.data['discarded'], i] = fresh_fuel[nuclides_list.index(int(i))]

        # Reset recirculated data
        pbed.data.loc[pbed.data['recirculated'], 'passes'] += 1 # add one pass
        pbed.data.loc[pbed.data['recirculated'], 'pass_residence_time'] = 0.0 # reset pass residence time
        pbed.data.loc[pbed.data['recirculated'], 'pass_burnup'] = 0.0 # reset pass burnup
        pbed.data.loc[pbed.data['recirculated'], 'agg_r_dist_pass'] = 0.0

        curr_time = float(next_time) # Increment time
        pbed.data['pass_residence_time'] += time_step  # Increment pass residence times
        pbed.data['residence_time'] += time_step  # Increment residence times

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
        # Extract burnups and calculate pass burnup
        pbed.data.loc[pbed.data['fuel'], 'burnup'] = Serpent_get_values(tra['bu_out'])
        if step > 0:
            # Pass burnup incremented for fuel pebbles
            pbed.data.loc[(pbed.data['fuel'] & (~pbed.data['recirculated'])), 'pass_burnup'] += pbed.data.loc[pbed.data['fuel'], 'burnup'] - pbed.old_data.loc[pbed.data['fuel'], 'burnup']

            # Special case for pebbles which recirculated
            pbed.data.loc[(pbed.data['fuel'] & (~pbed.data['recirculated'])), 'pass_burnup'] += pbed.data.loc[pbed.data['fuel'], 'burnup'] - pbed.old_data.loc[pbed.data['fuel'], 'burnup']

        # Extract isotopic inventory
        for name in list(inventory_names):
            pbed.data.loc[pbed.data['fuel'], name] = Serpent_get_values(tra[name])

        # Extract detector values and uncertainties and calculate integrated values
        for name in list(detector_names):
            # Tallies
            pbed.data[name] = Serpent_get_values(tra[name])
            pbed.data[f'{name}_rel_unc'] =  Serpent_get_values(tra[f'{name}_rel_unc'])

            # Integrated tallies
            if step > 0:
                pbed.data[f'integrated_{name}'] += time_step * DAYS * pbed.data[name]
                pbed.data[f'integrated_{name}_unc'] += (time_step * DAYS * pbed.data[name])*pbed.data[f'{name}_rel_unc']
                pbed.data[f'integrated_{name}_rel_unc'] = pbed.data[f'integrated_{name}_unc']/pbed.data[f'integrated_{name}']
                if 'power' in name: # Only fuel has power, rest is "nan"
                    pbed.data.loc[~pbed.data['fuel'], f'integrated_{name}'] = np.nan
                    pbed.data.loc[~pbed.data['fuel'], f'integrated_{name}_unc'] = np.nan
                    pbed.data.loc[~pbed.data['fuel'], f'integrated_{name}_rel_unc'] = np.nan

    # Calculate new average distance
    pbed.data['avg_r_dist'] = (pbed.data['avg_r_dist']*(step-pbed.data['insertion_step'])+pbed.data['r_dist'])/(step-pbed.data['insertion_step']+1)
    pbed.data['agg_r_dist_pass'] += pbed.data['r_dist']*time_step

    #### Plots ###
    if plotting:
        print('\tPlotting')

        # Plot latest geometry and save it to folder
        plot_Serpent(tra['plot'], plots_folder_path='wrk_Serpent', output_suffix=f'_{step}', output_folder='Plots', nplots=nplots, delay=60)

        # Create copy of table and fill non-fuel data with "nans"
        pbed_fuel = deepcopy(pbed)
        pbed_fuel.data.loc[~pbed.data['fuel'], pbed.data.columns.drop(['x', 'y', 'z', 'r', 'fuel'])] = np.nan
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
        plt.savefig(f'Plots/core_{step}.png', dpi=400, bbox_inches='tight')# save to output

        # Plot discarded pebbles
        if len(pbed.discarded_data)>0:
            plt.figure(figsize=(11.5, 15))
            pbed.discarded_data.drop(columns=[i for i in pbed.discarded_data if '_unc' in i]).astype(float).hist() # Plot distributions for each quantity
            plt.tight_layout()
            plt.savefig(f'Plots/discarded_{step}.png', dpi=400, bbox_inches='tight')

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
        pbed.cycle_hist.to_csv(f'Data/cycle_{step}.csv') # global inventory data
        pbed.data.to_csv(f'Data/data_{step}.csv') # whole core data
        pbed.discarded_data.to_csv(f'Data/discarded_{step}.csv') # discarded pebbles data
        pbed.reinserted_data.to_csv(f'Data/reinserted_{step}.csv') # re-inserted pebbles data
        pbed.reinserted_fuel_data.to_csv(f'Data/reinserted_fuel_{step}.csv') # re-inserted pebbles data, only fuel pebbles

        # Added for ML project: Cs137 vs keff
        nuclide = '551370'
        if transport and nuclide in pbed.reinserted_fuel_data.columns:
            pbed.tracked_reinserted_nuclide.loc[step, ['keff', 'keff_relative_uncertainty', 'keff_absolute_uncertainty']] = [keff, keff_rel_unc, keff_unc]
            tracked = pbed.reinserted_fuel_data[nuclide]
            tracked.index = [f'{ZAI_to_name(nuclide)}_{i} [at/b.cm]' for i in range(len(tracked))]
            pbed.tracked_reinserted_nuclide.loc[step, tracked.index] = tracked.T
            pbed.tracked_reinserted_nuclide.to_csv(f'Data/tracked_{ZAI_to_name(nuclide)}_{step}.csv')