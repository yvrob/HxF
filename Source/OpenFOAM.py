import os
import numpy as np
import time
import subprocess

def create_sbatch_file(step, substep, current_TH_step, TH_step_size, nnodes, partition, cpus_per_task, qos, account, time_lim=1, solver='GeN-Foam'):
    ntasks = int(cpus_per_task)*int(nnodes)
    formatted_time = '{:02d}:{:02d}:{:02d}'.format(int(time_lim), int((time_lim * 60) % 60), int((time_lim * 3600) % 60))
    job_name = f'HxF_TH_step{step}.{substep}'
    remove_mpi_cmd = 'unset "${!OMPI_@}" "${!MPI_@}"'
    first_block = ''
    if current_TH_step == 0:
        first_block += f'''source {os.environ["OF_BASHRC"]}
foamDictionary ./OF/system/decomposeParDict -entry numberOfSubdomains -set {ntasks}
foamDictionary ./OF/system/thermoMechanicalRegion/decomposeParDict -entry numberOfSubdomains -set {ntasks}
foamDictionary ./OF/system/neutroRegion/decomposeParDict -entry numberOfSubdomains -set {ntasks}
foamDictionary ./OF/system/fluidRegion/decomposeParDict -entry numberOfSubdomains -set {ntasks}
foamDictionary ./OF/system/controlDict -entry writeInterval -set {TH_step_size}
decomposePar -case ./OF -allRegions -latestTime -force
postProcess -func writeCellVolumes -case ./OF -latestTime 
postProcess -func writeCellCentres -case ./OF -latestTime 
'''

    config = f'''#!/bin/bash
#SBATCH --job-name={job_name}
#SBATCH --output={job_name}.o
#SBATCH --error={job_name}.error

#SBATCH --account={account}
#SBATCH --qos={qos}
#SBATCH --partition={partition}
#SBATCH --time={formatted_time}

#SBATCH --nodes={nnodes}
#SBATCH --cpus-per-task={cpus_per_task}
#SBATCH --ntasks-per-node=1
'''
    content = f'''
###### Script ######
# Initial actions
cd ./OF
{remove_mpi_cmd}

# Add block to power density
block=$(cat <<ETH
boundaryField
{{
    Inlet
    {{
        type            zeroGradient;
    }}
    Outlet
    {{
        type            zeroGradient;
    }}
    Wall
    {{
        type            zeroGradient;
    }}
}}
ETH
)

if ! grep -q 'boundaryField' "{current_TH_step}/fluidRegion/powerDensityNeutronics"; then
    cat "{current_TH_step}/fluidRegion/powerDensityNeutronics" > "{current_TH_step}/fluidRegion/powerDensityNeutronics.tmp"
    echo "$block" >> "{current_TH_step}/fluidRegion/powerDensityNeutronics.tmp"
    mv "{current_TH_step}/fluidRegion/powerDensityNeutronics.tmp" "{current_TH_step}/fluidRegion/powerDensityNeutronics"
fi

# Prepare TH case

mpirun -np {ntasks} -oversubscribe --bind-to none {solver} -parallel
wait
reconstructPar -allRegions -latestTime
cd ..
'''
    return config + first_block + content

def execute_GeN_Foam(of_step, of_iteration, of_substep, TH):
    with open('execute_TH.sub', 'w') as fw:
        fw.write(create_sbatch_file(of_iteration, of_substep, of_step, TH['step_size'], TH['nnodes'], os.environ['PARTITION'], os.environ['PARTITION_CPUS_PER_NODE'], os.environ['QOS'], os.environ['GROUP_ACCOUNT'], TH['time_limit'], TH['solver']))
    output = subprocess.check_output('sbatch execute_TH.sub', shell=True, universal_newlines=True)
    job_id = output.strip().split()[-1]
    while True:
        output = subprocess.check_output(f'squeue -j {job_id}', shell=True, universal_newlines=True)
        if job_id not in output:
            break
        time.sleep(5) 
    print(f"Job {job_id} completed.")
    

def read_GeN_Foam(name_field, folder='0/', path_case='./', region='', dtype=float):
    field_path = f'{path_case}/{folder}/{region}/{name_field}'

    with open(field_path, 'r') as f:
        lines = f.readlines()

    for i, line in enumerate(lines):
        if 'internalField' in line and 'nonuniform' not in line:
            line = line.strip()[:-1].split('uniform')[-1].strip()
            return float(line)

        if 'internalField' in line:
            start_index = i + 3
            break

    nvalues = int(lines[start_index - 2])
    dim = len(lines[start_index].replace('(', '').replace(')', '').split())
    field = np.empty((nvalues, dim), dtype=dtype)

    for i, line in enumerate(lines[start_index:start_index+nvalues]):
        field[i] = np.array(line.replace('(', '').replace(')', '').split()).astype(dtype)

    if dim == 1:
        field = field.ravel()

    return field

# def init_TH(TH):
#     source {os.environ["OF_BASHRC"]}
