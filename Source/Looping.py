import pandas as pd
import numpy  as np
from joblib import Parallel, delayed, cpu_count
import os
import re

#%% Functions
def natural_sort(l):
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
    return sorted(l, key=alphanum_key)

def make_equal_zones(data_target, data_source, Nr, Nz):
    grouped_data_target = [group for _, group in data_target.groupby(pd.qcut(data_target['r_dist'], Nr, labels=False))]
    grouped_data_source = [group for _, group in data_source.groupby(pd.qcut(data_source['r_dist'], Nr, labels=False))]
    
    # Radial match
    for i in range(Nr-1):
        if len(grouped_data_target[i]) != len(grouped_data_source[i]):
            if len(grouped_data_target[i]) > len(grouped_data_source[i]):
                while len(grouped_data_target[i]) > len(grouped_data_source[i]):
                    maxr_id = grouped_data_target[i]['r_dist'].idxmax()
                    grouped_data_target[i+1] = pd.concat([grouped_data_target[i+1], pd.DataFrame(grouped_data_target[i].loc[maxr_id]).T],axis=0)
                    grouped_data_target[i] = grouped_data_target[i].drop(maxr_id)
            else:
                while len(grouped_data_target[i]) < len(grouped_data_source[i]):
                    maxr_id = grouped_data_source[i]['r_dist'].idxmax()
                    grouped_data_source[i+1] = pd.concat([grouped_data_source[i+1], pd.DataFrame(grouped_data_source[i].loc[maxr_id]).T],axis=0)
                    grouped_data_source[i] = grouped_data_source[i].drop(maxr_id)

        grouped_data_target[i] = grouped_data_target[i].sort_index()
        grouped_data_target[i+1] = grouped_data_target[i+1].sort_index()
        grouped_data_source[i] = grouped_data_source[i].sort_index()
        grouped_data_source[i+1] = grouped_data_source[i+1].sort_index()

    for i in range(Nr):
        if (len(grouped_data_target[i]) != len(grouped_data_source[i])):
            raise Exception('Problem with radial zones assignment')
    
    # Create axial groups within each radial group
    grouped_data_target = [group.groupby(pd.qcut(group['z'], Nz, labels=False)) for group in grouped_data_target]
    grouped_data_source = [group.groupby(pd.qcut(group['z'], Nz, labels=False)) for group in grouped_data_source]
    for i in range(len(grouped_data_target)):
        grouped_data_target[i] = [group for _, group in grouped_data_target[i]]
        grouped_data_source[i] = [group for _, group in grouped_data_source[i]]

    # Axial match
    for i in range(Nr):
        for j in range(Nz-1):
            if len(grouped_data_target[i][j]) != len(grouped_data_source[i][j]):
                if len(grouped_data_target[i][j]) > len(grouped_data_source[i][j]):
                    while len(grouped_data_target[i][j]) > len(grouped_data_source[i][j]):
                        maxz_id = grouped_data_target[i][j]['z'].idxmax()
                        grouped_data_target[i][j+1] = pd.concat([grouped_data_target[i][j+1], pd.DataFrame(grouped_data_target[i][j].loc[maxz_id]).T],axis=0)
                        grouped_data_target[i][j] = grouped_data_target[i][j].drop(maxz_id)
                else:
                    while len(grouped_data_target[i][j]) < len(grouped_data_source[i][j]):
                        maxz_id = grouped_data_source[i][j]['z'].idxmax()
                        grouped_data_source[i][j+1] = pd.concat([grouped_data_source[i][j+1], pd.DataFrame(grouped_data_source[i][j].loc[maxz_id]).T],axis=0)
                        grouped_data_source[i][j] = grouped_data_source[i][j].drop(maxz_id)

            grouped_data_target[i][j] = grouped_data_target[i][j].sort_index()
            grouped_data_target[i][j+1] = grouped_data_target[i][j+1].sort_index()
            grouped_data_source[i][j] = grouped_data_source[i][j].sort_index()
            grouped_data_source[i][j+1] = grouped_data_source[i][j+1].sort_index()

    zones1 = np.empty((Nz, Nr), dtype=object)
    zones2 = np.empty((Nz, Nr), dtype=object)
    zones_indices1 = np.empty((Nz, Nr), dtype=np.ndarray)
    zones_indices2 = np.empty((Nz, Nr), dtype=np.ndarray)
    for i in range(zones1.shape[0]):
        for j in range(zones1.shape[1]):
            zones1[i, j] = grouped_data_target[j][i]
            zones_indices1[i, j] = np.array(zones1[i, j].index.tolist())
            zones2[i, j] = grouped_data_source[j][i]
            zones_indices2[i, j] = np.array(zones2[i, j].index.tolist())
    for i in range(zones1.shape[0]):
        for j in range(zones1.shape[1]):
            if zones1[i, j].shape != zones2[i,j].shape:
                raise Exception('Problem with axial zones assignment')
    return zones1, zones2, zones_indices1, zones_indices2

def match_points_parallel(i_z, i_r, zones1, zones2, method, exponent):
    message = f'\t[{(os.getpid() - 1) % cpu_count()}] Matching positions for radial zone {i_r}, axial zone {i_z}'
    print(message, flush=True)
    return match_points(zones1[i_z, i_r], zones2[i_z, i_r], method=method, exponent=exponent)

def match_points(points1, points2, method='rz', exponent=1):
    from scipy.optimize import linear_sum_assignment
    # Two methods: minimize (x,y,z) or (r_dist, z) distance
    if method == 'xyz':
        cost = np.linalg.norm(points2[['x', 'y', 'z']].values[:, np.newaxis, :] - points1[['x', 'y', 'z']].values, axis=2)
    elif method == 'rz':
        rz1 = points1[['r_dist', 'z']]
        rz2 = points2[['r_dist', 'z']]
        cost = np.linalg.norm(rz2.values[:, np.newaxis, :] - rz1.values, axis=2)
    cost = cost**exponent
    _, indices = linear_sum_assignment(cost)
    points3 = points1.iloc[indices]
    return points3, indices, cost

def make_transition(zones1, zones2, method='rz', exponent=1, workers=-1):
    Nz, Nr = zones1.shape
    transition = np.zeros(np.sum([zones1[i][j].shape[0] for i in range(Nz) for j in range(Nr)]), dtype=int) # initialize with zeros

    tasks = [(i_z, i_r) for i_z in range(Nz) for i_r in range(Nr)]
    results = Parallel(n_jobs=cpu_count())(
        delayed(match_points_parallel)(*task, zones1, zones2, method, exponent) for task in tasks
    )

    for i, (points, indices, _) in enumerate(results):
        transition[zones2[i // Nr, i % Nr].index] = points.index

    if len(np.unique(transition)) != len(transition):
        raise Exception("Wrong transition!")
    print('\tDone.')
    return transition

def looper(file_start, file_end, Nr, Nz, method, workers=-1):
    print('Setting up loop:')
    # Load positions
    print(f'\tFirst DEM step: {file_start}')
    data_start = pd.read_csv(file_start)
    data_start['r_dist'] = np.linalg.norm(data_start[['x', 'y']], axis=1)
    data_start['id'] = np.arange(data_start.shape[0])
    print(f'\tLast DEM step: {file_end}')
    data_end = pd.read_csv(file_end)
    data_end['r_dist'] = np.linalg.norm(data_end[['x', 'y']], axis=1)
    data_end['id'] = np.arange(data_end.shape[0])
    print(f'\tNumber of zones: {Nr} radial and {Nz} axial zones')
    print(f'\tOptimization method: {method}')    
    zones1, zones2, zones_indices1, zones_indices2 = make_equal_zones(data_start, data_end, Nr, Nz)
    transition = make_transition(zones1, zones2, method=method, exponent=2, workers=workers)
    print((data_start.loc[transition].reset_index(drop=True)-data_end).abs()[['z', 'r_dist']].describe().T)
    print("Done.")
    return transition

# def make_equal_zones(data, looper_Nz, looper_Nr):
#     cut_data = pd.qcut(data['r_dist'].rank(method = 'first'), looper_Nr)
#     grouped_data = data.groupby(cut_data)
#     radial_groups = [group for _, group in grouped_data]
#     zones = np.empty((looper_Nz, looper_Nr), dtype=object)
#     zones_indices = np.empty((looper_Nz, looper_Nr), dtype=list)
#     for i_r, rad_group in enumerate(radial_groups):
#         cut_data = pd.qcut(rad_group['z'], looper_Nz)
#         grouped_data = rad_group.groupby(cut_data)
#         axial_groups = [group for _, group in grouped_data]
#         for i_z in range(looper_Nz):
#             zones[i_z, i_r] = axial_groups[i_z]
#             zones_indices[i_z, i_r] = axial_groups[i_z].index.values
#     return zones_indices, zones

# def match_equal_zones(zones1, zones2, zones_indices2):
#     looper_Nz, looper_Nr = zones1.shape
#     not_match = True
#     ntrys_max = 1000
#     ntrys = 0
#     while not_match:
#         if ntrys == ntrys_max:
#             raise(Exception('Problem to make equal zones'))
#         not_match = False
#         not_matching = []
#         for i in range(looper_Nz):
#             for j in range(looper_Nr):
#                 if zones1[i,j].shape[0] != zones2[i,j].shape[0]:
#                     not_matching.append((i, j))
#                     not_match = True
#         #print(not_matching)
#         for i in not_matching:
#             for j in not_matching:
#                 if i[0] > j[0] and i[1]==j[1]:
#                     dif = zones1[i].shape[0] - zones2[i].shape[0]
#                     k = (i[0]-1, i[1])
#                     if dif > 0:
#                         # More in zone 1 than zone 2, take from zone below the highest point and transfer it
#                         max_z_id = zones2[k]['z'].argmax()
#                         zones2[i].loc[zones2[k].iloc[max_z_id].name] = zones2[k].iloc[max_z_id]
#                         zones_indices2[i] = np.hstack([zones_indices2[i], zones_indices2[k][max_z_id]])
#                         # Remove from original zone
#                         zones2[k].drop(zones2[k].iloc[max_z_id].name)
#                         zones_indices2[k] = np.delete(zones_indices2[k], max_z_id, axis=0)
#                     elif dif < 0:
#                         # More in zone 2 than zone 1, get the lowest point in zone and transfer it to zone below
#                         min_z_id = zones2[k]['z'].argmin()
#                         zones2[k].loc[zones2[i].iloc[min_z_id].name] = zones2[i].iloc[min_z_id]
#                         zones_indices2[k] = np.hstack([zones_indices2[k], zones_indices2[i][min_z_id]])
#                         # Remove from original zone
#                         zones2[i].drop(zones2[i].iloc[min_z_id].name)
#                         zones_indices2[i] = np.delete(zones_indices2[i], min_z_id, axis=0)
#         ntrys += 1
#     return zones1, zones2, zones_indices2

# def match_points(points1, points2, method='rz'):
#     import scipy.optimize
#     if method == 'xyz':
#         cost = np.linalg.norm(points2[['x', 'y', 'z']].values[:, np.newaxis, :] - points1[['x', 'y', 'z']].values, axis=2)
#     elif method == 'rz':
#         rz1 = points1[['r_dist', 'z']]
#         rz2 = points2[['r_dist', 'z']]
#         cost = np.linalg.norm(rz2.values[:, np.newaxis, :] - rz1.values, axis=2)
#     _, indices = scipy.optimize.linear_sum_assignment(cost)
#     #print(indices.shape, points1.shape)
#     points3 = points1.iloc[indices]
#     return points3, indices, cost

# def make_transition(zones1, zones2, zones_indices2, data_start, method='rz'):
#     zones3 = np.empty_like(zones1)
#     looper_Nz, looper_Nr = zones1.shape
#     for i_z in range(looper_Nz):
#         for i_r in range(looper_Nr):
#             points, indices, cost = match_points(zones1[i_z, i_r], zones2[i_z, i_r], method=method)
#             zones3[i_z, i_r] = points
#     indices = np.concatenate(np.concatenate(zones_indices2))
#     points = pd.DataFrame(np.concatenate(np.concatenate(zones3)), columns=data_start.columns)
#     transition = points.loc[indices.argsort(), 'id'].astype(int).values
#     return transition

# def looper(files, DEM_start, DEM_end, Nr, Nz, method):
#     first_start_step = int(DEM_start)
#     first_end_step = int(DEM_end)
#     found = False
#     for DEM_end in range(first_end_step, first_start_step, -1):
#         for DEM_start in range(first_start_step, DEM_end):
#             print(DEM_start, DEM_end)
#             data_start = pd.read_csv(files[DEM_start])
#             data_start['r_dist'] = np.linalg.norm(data_start[['x', 'y']], axis=1)
#             data_start['id'] = np.arange(data_start.shape[0])

#             data_end = pd.read_csv(files[DEM_end])
#             data_end['r_dist'] = np.linalg.norm(data_end[['x', 'y']], axis=1)
#             data_end['id'] = np.arange(data_end.shape[0])
#             try:
#                 zones_indices1, zones1 = make_equal_zones(data_start, Nz, Nr)
#                 zones_indices2, zones2 = make_equal_zones(data_end, Nz, Nr)
#                 zones1, zones2, zones_indices2 = match_equal_zones(zones1, zones2, zones_indices2)
#                 transition = make_transition(zones1, zones2, zones_indices2, data_start, method=method)
#                 found = True
#                 break
#             except:
#                 pass
#         if found:
#             break
#     data_matched = data_start.loc[transition].reset_index(drop=True)
#     data_matched['id'] = data_matched.index

#     errors_rz = data_matched[['r_dist', 'z']] - data_end[['r_dist', 'z']]
#     errors = np.linalg.norm(data_matched[['r_dist', 'z']] - data_end[['r_dist', 'z']], axis=1)

#     error_rz = np.linalg.norm(errors_rz)
#     return DEM_start, DEM_end, transition, error_rz

if __name__ == '__main__':
    from Utilities import natural_sort
    import os
    from glob import glob

    step = 750
    positions_folder = '/global/scratch/users/co_nuclear/pebble_positions_larger/'

    print(f'Using DEM motion from folder {positions_folder}')
    if not os.path.exists(positions_folder):
        raise Exception(f'DEM motion selected but positions folder {positions_folder} does not exist.')
    position_files = natural_sort(glob(positions_folder+'/step*.csv')) # List all positions found with DEM
    DEM_start, DEM_end, transition, error_rz = looper(position_files, DEM_start, DEM_end, looper_Nr, looper_Nz, looper_method)
    print(f'Looper ready: from DEM step {DEM_start} to {DEM_end} (err={error_rz:.2f})')

    # Just for checking, check that positions in data correspond to the right positions
    nsteps_to_loop = (DEM_end-DEM_start)/DEM_step_increment # won't work if different step increments! after how many steps do we loop
    if step<=nsteps_to_loop: # if first loop (original), no modification
        nloops = 0
        equivalent_step = step
    else: # need to have an equivalent step, and apply transition indices as many times as there were loops
        nloops = int((step-1)//nsteps_to_loop)
        equivalent_step = int((step-1)%nsteps_to_loop)+1
    
    print(f'\tFirst check: Using positions at file: {position_files[equivalent_step]} (loop #{nloops}, equivalent step #{equivalent_step})')
    positions = pd.read_csv(position_files[equivalent_step])[['x','y','z']]*positions_scale + np.array(positions_translation) # read new positions from DEM file
    
    # Apply looping transition, if needed
    indices = np.arange(positions.shape[0], dtype=int)
    for i in range(nloops): # does not go here if nloops=0
        indices = indices[transition]
    positions[['x', 'y', 'z']] = positions.loc[indices][['x', 'y', 'z']].values
    print(positions)