#%% Modules
print('Importing modules')
import math
import numpy as np
import struct
import os
import concurrent.futures
import shutil
from glob import glob
import pandas as pd
from time import time
from itertools import repeat
import pickle
import sys
import multiprocessing

def PolarAngle(position):
    x = position[0]
    y = position[1]
    if x == 0.0 and y > 0.0:
        t = np.pi / 2.0
    elif x == 0.0 and y < 0.0:
        t = 3.0 * np.pi / 2.0
    elif x > 0.0 and y == 0.0:
        t = 0.0
    elif x < 0.0 and y == 0.0:
        t = np.pi
    elif x > 0.0 and y > 0.0:
        t = math.atan(y / x)
    elif x < 0.0 and y > 0.0:
        t = math.atan(y / x) + np.pi
    elif x < 0.0 and y < 0.0:
        t = math.atan(y / x) + np.pi
    elif x > 0.0 and y < 0.0:
        t = math.atan(y / x) + 2.0 * np.pi
    else:
        return 0.0
    return t


def replicate_dd(data, mpitasks):
    print("Replicating Serpent domain decomposition")
    positions = np.array(data[["x", "y", "z"]])

    subn = 100000
    nr = np.zeros((subn), int)
    nmat = len(positions)

    # Assign angular segment
    for i in positions:
        f = PolarAngle(i) / (2 * np.pi)
        f_truncated = f - float(int(f))
        n = int((float(subn) - 1.0) * f_truncated)
        nr[n] += 1

    # Assign MPI ID to segments
    domain_id = 0
    s = 0
    for n in range(subn):
        s += nr[n]
        nr[n] = domain_id
        if s > int(float(nmat) / float(mpitasks)) - 1:
            s = 0
            if domain_id < mpitasks - 1:
                domain_id += 1

    # Assign MPI ID to materials
    domain_id_list = []
    domains_indices = [[] for i in range(mpitasks)]
    for i, pos in enumerate(positions):
        f = PolarAngle(pos) / (2 * np.pi)
        f_truncated = f - float(int(f))
        n = int((float(subn) - 1.0) * f_truncated)
        if nr[n] < 0 or nr[n] > mpitasks - 1:
            raise Exception("indexing error: {} {}".format(nr[n], n))
        else:
            domain_id_list.append(nr[n])
            domains_indices[nr[n]].append(i)

    # Fill domains and check number of materials per domain
    s = 0
    MPI_cnt = [0 for i in range(mpitasks)]
    for domain_id in range(mpitasks):
        for i in domain_id_list:
            if i == domain_id:
                MPI_cnt[domain_id] += 1
        if MPI_cnt[domain_id] == 0:
            raise Exception(
                "Decomposition failed: no materials assigned to task {}".format(
                    domain_id + 1
                )
            )
        s += MPI_cnt[domain_id]
        print(
            "\tDomain {}: {} materials ({:.3f}%)".format(
                domain_id + 1, MPI_cnt[domain_id], 100.0 * (MPI_cnt[domain_id] / nmat)
            )
        )

    if s < nmat:
        print(0, "{} materials not decomposed".format(nmat - s))
    elif s > nmat:
        raise Exception("Not possible")

    # Assign to pebble
    data["domain_id"] = domain_id_list
    print("\tDone")
    return data


def write_domain_restart(data, domain_id, fuel_material, dep_output, files_out):
    pc = 0
    t0 = time()
    data = data[data["domain_id"] == domain_id].sort_index(ascending=False)
    with open(dep_output, "rb") as f:
        interpolator, list_iso = pickle.load(f)
    with open(files_out[domain_id], "wb") as fo:
        # Write actual materials
        cnt = 0
        pc_inc = 10

        # Zones
        for i in range(len(data)):
            mat = data.iloc[i]
            name = "{}z{}".format(fuel_material, mat.name + 1)
            adens = interpolator(mat["burnup"])
            nnuc = adens.shape[0] - 1
            content = b""
            content += struct.pack("q", len(name))
            content += struct.pack("{}s".format(len(name)), str.encode(name))
            content += struct.pack("d", 0.0)  # BU global
            content += struct.pack("d", 0.0)  # BU days
            content += struct.pack("q", nnuc)
            content += struct.pack("d", adens[-1])
            content += struct.pack("d", 0.0)  # mdens
            content += struct.pack("d", mat["burnup"])
            content += struct.pack("q", -1)
            content += struct.pack("d", adens[-2])

            for j in range(len(adens) - 2):
                content += struct.pack("q", list_iso[j])
                content += struct.pack("d", adens[j])
            fo.write(content)

            # Timer
            current_pc = (cnt + 1) / len(data) * 100
            if current_pc >= pc:
                if pc == 0:
                    print("\t[{}]\t{}%".format(domain_id, pc))
                else:
                    dt = int(time() - t0)
                    estimation = int(dt / (cnt + 1) * len(data))
                    print(
                        "\t[{}]\t{}%. Elapsed time: {}s. Estimation: {}s. Remaining time: {}s".format(
                            domain_id, pc, dt, estimation, estimation - dt
                        )
                    )
                while pc <= current_pc:
                    pc += pc_inc

            cnt += 1
        dt = int(time() - t0)
        print("\t[{}]\t100%. Elapsed time: {}s.".format(domain_id, dt))


def create_decomposed_restart_files(data, mpitasks, parallel, dep_output, path="./", name_out="first_compos", fuel_material="fuel"):
    path = path + "/restart"
    if os.path.exists(path):
        shutil.rmtree(path)
    os.makedirs(path)

    files_out = []
    for i in range(mpitasks):
        files_out.append("{}/{}.wrk_dd{}".format(path, name_out, i))
        try:
            os.remove(files_out[-1])
        except:
            pass
    print(files_out)
    print("Creating decomposed restart files")
    if parallel:
        MPItasks = [i for i in range(mpitasks)]
        with concurrent.futures.ProcessPoolExecutor() as executor:
            for domain_id in zip(MPItasks, executor.map(write_domain_restart, repeat(data), MPItasks, repeat(fuel_material), repeat(dep_output), repeat(files_out))):
                pass
    else:
        for domain_id in range(mpitasks):
            write_domain_restart(data, domain_id, fuel_material, dep_output, files_out)


#%% Input
print('Reading input')
parallel_creation = True
fuel_material_name = "fuel_1"
interpolator_data = "gFHR_equilibrium.interpolator"
name_out = "first_compos"
positions_file = '../../Inputs/gFHR/fpb_pos'
max_burnup = 100  # MWd/kgHM
npasses = 8
residence_time = 522  # days
fuel_frac = 1
fuel_pebble_universe_name = 'u_fuel_pebble'
motion_direction = +1
output_path = './'

#%% Process input
print('Processing input')
time_per_pass = residence_time / npasses
if parallel_creation:
    nfiles = multiprocessing.cpu_count()

#%% Create dataframe

# Read positions
print('Reading positions file')
data = pd.read_csv(positions_file, sep="\t", header=None, names=["x", "y", "z", "r_pebbles", "uni"])
data["r_dist"] = np.linalg.norm(data[["x", "y"]], axis=1)

# Assign universes
print('Assigning initial universes')
if fuel_frac != 1:
    data['uni'] = assign_random_array(np.arange(len(data)), [fuel_pebble_universe_name, graph_pebble_universe_name], [fuel_frac, 1-fuel_frac]) # Spl
else:
    data['uni'] = fuel_pebble_universe_name

data['fuel'] = ['fuel' in u for u in data['uni']]
Nfuel = sum(['fuel' in u for u in data['uni']])

# Initialize passes
print('Initializing passes')
data.loc[data['fuel'], 'passes'] = np.random.randint(1, npasses+1, Nfuel)

# Passes and z to residence time estimate
print('Estimating residence times from passes and elevation')
previous_passes_time = time_per_pass * (data.passes - 1)
H = data.z.max() - data.z.min()
if motion_direction == +1:
    dz = data.z - data.z.min()
elif motion_direction == -1:
    dz = data.z.max() - data.z
data['pass_residence_time'] = dz / H * time_per_pass
data['residence_time'] = previous_passes_time + data['pass_residence_time']

# Residence time to burnup estimate, assuming linear burnup (rough assumption)
print('Estimating burnups from residence times')
data["burnup"] = data["residence_time"] / residence_time * max_burnup

# Make domains for parallel processing only
print(f'Decomposing into {nfiles} sector domains')
data = replicate_dd(data, nfiles)

# Save csv
print(data.describe(percentiles=[]).T,'\n')
print('Saving csv file')
data.to_csv(f'data_{name_out}.csv')

# Make restart files
print('Making restart files file')
create_decomposed_restart_files(data, nfiles, parallel_creation, interpolator_data, path=output_path, name_out=name_out, fuel_material=fuel_material_name)

