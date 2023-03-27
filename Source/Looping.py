import numpy as np
import pandas as pd
import re

#%% Functions
def natural_sort(l):
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
    return sorted(l, key=alphanum_key)

def make_equal_zones(data, looper_Nz, looper_Nr):
    cut_data = pd.qcut(data['r_dist'].rank(method = 'first'), looper_Nr)
    grouped_data = data.groupby(cut_data)
    radial_groups = [group for _, group in grouped_data]
    zones = np.empty((looper_Nz, looper_Nr), dtype=object)
    zones_indices = np.empty((looper_Nz, looper_Nr), dtype=list)
    for i_r, rad_group in enumerate(radial_groups):
        cut_data = pd.qcut(rad_group['z'], looper_Nz)
        grouped_data = rad_group.groupby(cut_data)
        axial_groups = [group for _, group in grouped_data]
        for i_z in range(looper_Nz):
            zones[i_z, i_r] = axial_groups[i_z]
            zones_indices[i_z, i_r] = axial_groups[i_z].index.values
    return zones_indices, zones

def match_equal_zones(zones1, zones2, zones_indices2):
    looper_Nz, looper_Nr = zones1.shape
    not_match = True
    ntrys_max = 1000
    ntrys = 0
    while not_match:
        if ntrys == ntrys_max:
            raise(Exception('Problem to make equal zones'))
        not_match = False
        not_matching = []
        for i in range(looper_Nz):
            for j in range(looper_Nr):
                if zones1[i,j].shape[0] != zones2[i,j].shape[0]:
                    not_matching.append((i, j))
                    not_match = True
        #print(not_matching)
        for i in not_matching:
            for j in not_matching:
                if i[0] > j[0] and i[1]==j[1]:
                    dif = zones1[i].shape[0] - zones2[i].shape[0]
                    k = (i[0]-1, i[1])
                    if dif > 0:
                        # More in zone 1 than zone 2, take from zone below the highest point and transfer it
                        max_z_id = zones2[k]['z'].argmax()
                        zones2[i].loc[zones2[k].iloc[max_z_id].name] = zones2[k].iloc[max_z_id]
                        zones_indices2[i] = np.hstack([zones_indices2[i], zones_indices2[k][max_z_id]])
                        # Remove from original zone
                        zones2[k].drop(zones2[k].iloc[max_z_id].name)
                        zones_indices2[k] = np.delete(zones_indices2[k], max_z_id, axis=0)
                    elif dif < 0:
                        # More in zone 2 than zone 1, get the lowest point in zone and transfer it to zone below
                        min_z_id = zones2[k]['z'].argmin()
                        zones2[k].loc[zones2[i].iloc[min_z_id].name] = zones2[i].iloc[min_z_id]
                        zones_indices2[k] = np.hstack([zones_indices2[k], zones_indices2[i][min_z_id]])
                        # Remove from original zone
                        zones2[i].drop(zones2[i].iloc[min_z_id].name)
                        zones_indices2[i] = np.delete(zones_indices2[i], min_z_id, axis=0)
        ntrys += 1
    return zones1, zones2, zones_indices2

def match_points(points1, points2, method='rz'):
    import scipy.optimize
    if method == 'xyz':
        cost = np.linalg.norm(points2[['x', 'y', 'z']].values[:, np.newaxis, :] - points1[['x', 'y', 'z']].values, axis=2)
    elif method == 'rz':
        rz1 = points1[['r_dist', 'z']]
        rz2 = points2[['r_dist', 'z']]
        cost = np.linalg.norm(rz2.values[:, np.newaxis, :] - rz1.values, axis=2)
    _, indices = scipy.optimize.linear_sum_assignment(cost)
    #print(indices.shape, points1.shape)
    points3 = points1.iloc[indices]
    return points3, indices, cost

def make_transition(zones1, zones2, zones_indices2, data_start, method='rz'):
    zones3 = np.empty_like(zones1)
    looper_Nz, looper_Nr = zones1.shape
    for i_z in range(looper_Nz):
        for i_r in range(looper_Nr):
            points, indices, cost = match_points(zones1[i_z, i_r], zones2[i_z, i_r], method=method)
            zones3[i_z, i_r] = points
    indices = np.concatenate(np.concatenate(zones_indices2))
    points = pd.DataFrame(np.concatenate(np.concatenate(zones3)), columns=data_start.columns)
    transition = points.loc[indices.argsort(), 'id'].astype(int).values
    return transition

def looper(files, DEM_start, DEM_end, Nr, Nz, method):
    first_start_step = int(DEM_start)
    first_end_step = int(DEM_end)
    found = False
    for DEM_end in range(first_end_step, first_start_step, -1):
        for DEM_start in range(first_start_step, DEM_end):
            print(DEM_start, DEM_end)
            data_start = pd.read_csv(files[DEM_start])
            data_start['r_dist'] = np.linalg.norm(data_start[['x', 'y']], axis=1)
            data_start['id'] = np.arange(data_start.shape[0])

            data_end = pd.read_csv(files[DEM_end])
            data_end['r_dist'] = np.linalg.norm(data_end[['x', 'y']], axis=1)
            data_end['id'] = np.arange(data_end.shape[0])
            try:
                zones_indices1, zones1 = make_equal_zones(data_start, Nz, Nr)
                zones_indices2, zones2 = make_equal_zones(data_end, Nz, Nr)
                zones1, zones2, zones_indices2 = match_equal_zones(zones1, zones2, zones_indices2)
                transition = make_transition(zones1, zones2, zones_indices2, data_start, method=method)
                found = True
                break
            except:
                pass
        if found:
            break
    data_matched = data_start.loc[transition].reset_index(drop=True)
    data_matched['id'] = data_matched.index

    errors_rz = data_matched[['r_dist', 'z']] - data_end[['r_dist', 'z']]
    errors = np.linalg.norm(data_matched[['r_dist', 'z']] - data_end[['r_dist', 'z']], axis=1)

    error_rz = np.linalg.norm(errors_rz)
    return DEM_start, DEM_end, transition, error_rz

