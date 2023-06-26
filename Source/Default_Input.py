##### Default Input parameters #####
# Input file
#filename_input = 'Input'

# Pebbles
fuel_frac = 1

#%% Depletion steps
power_normalization_field = 'power'
decay_step = 0 # days
transport = True
resolve_first = True
correct = False
use_decnfy_lib = True

#%% Outputing
verbosity_mode = 0
plotting = True
saving = True
write_global = True
write_incore = True
write_reinserted = True
write_discarded = True
write_discharged = True
inventory_names = ""

#%% Domain decomposition
domain_decomposition = False
allowed_decomposition_types = 'ars'
max_domains = [6, 4, 10]

#%% Configuration
from multiprocessing import cpu_count
ncores = cpu_count()
nnodes = 1 

#%% Restart read/write
restart_calculation = False
read_firt_compositions = False
write_restart = False
write_restart_discarded = False
write_restart_discharged = False
restart_write_every = 1
restart_discharged_write_every = 1
restart_discarded_write_every = 1