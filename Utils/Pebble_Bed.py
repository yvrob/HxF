import serpentTools
import pandas as pd
import numpy as np
from matplotlib.patches import Circle
from matplotlib.collections import PatchCollection
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import Normalize, LogNorm
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import ticker
from copy import deepcopy
from glob import glob
import matplotlib.image as mpimg
import os
import difflib
#%matplotlib widget

class Pebble_bed:

    #### GENERAL ####

    def __init__(self, verbose=True, level=0, log_delim='  '):
        log_print(f'Creating empty Pebble_bed object', verbose, level, log_delim)
        self.read_files = []

    def __repr__(self) -> str:
        if hasattr(self, 'data'):
            return self.data.__repr__()
        else:
            return 'Empty Pebble_bed object'

    #### READING ####

    def read_pbed_file(self, pbed_file_path, *, calculate_dist=True, calculate_radial_dist=True, radial_center=[0,0], verbose=True, level=0, log_delim='  '):
        log_print(f'Reading pbed file from {pbed_file_path}', verbose, level, log_delim)
        data = pd.read_csv(pbed_file_path, delim_whitespace=True, header=None, names=["x", "y", "z", "r", "uni"])
        self.pbed_path = pbed_file_path
        self.data = data
        self.data['id'] = np.array(self.data.index)
        self.process_data(calculate_dist, calculate_radial_dist, radial_center, verbose=verbose, level=level+1, log_delim=log_delim)
        self.read_files.append(pbed_file_path)
        log_print(f'Done.', verbose, level, log_delim, end_block=True)
        return data

    def read_datafile(self, data_file_path, sep=',', names='infer', first_column_index=True, r=None, calculate_dist=True, calculate_radial_dist=True, calculate_center=True, radial_center=[0,0], verbose=True, level=0, log_delim='  '):
        log_print(f'Reading data from data file {data_file_path}', verbose, level, log_delim)
        if first_column_index:
            index_col = 0
        else:
            index_col = None
        if not sep:
            self.data = pd.read_csv(data_file_path, delim_whitespace=True, header=names, index_col=index_col)
        else:
            self.data = pd.read_csv(data_file_path, sep=sep, header=names, index_col=index_col)

        if not hasattr(self.data, 'r'):
            self.data['r'] = r
        if not hasattr(self.data, 'id'):
            self.data['id'] = np.array(self.data.index)
        self.process_data(calculate_dist, calculate_radial_dist, radial_center, calculate_center=calculate_center, verbose=verbose, level=level+1, log_delim=log_delim)
        log_print(f'Done.', verbose, level, log_delim, end_block=True)

    def read_dataframe(self, dataframe, r=None, calculate_dist=True, calculate_radial_dist=True, calculate_center=True, calculate_angles=True, radial_center=[0,0], verbose=True, level=0, log_delim='  '):
        log_print(f'Reading data from data table', verbose, level, log_delim)
        self.data = dataframe
        if not hasattr(self.data, 'r'):
            self.data['r'] = r
        if not hasattr(self.data, 'id'):
            self.data['id'] = np.array(self.data.index)
        self.process_data(calculate_dist=calculate_dist, calculate_radial_dist=calculate_radial_dist, radial_center=radial_center, calculate_center=calculate_center, calculate_angles=calculate_angles, verbose=verbose, level=level+1, log_delim=log_delim)
        log_print(f'Done.', verbose, level, log_delim, end_block=True)

    def read_pos_table(self, table, universes_included=False, calculate_dist=True, calculate_radial_dist=True, radial_center=[0,0], verbose=True, level=0, log_delim='  '):
        log_print(f'Generating data from table containing xyzr', verbose, level, log_delim)
        if isinstance(table, pd.DataFrame):
            self.data = table
        else:
            if universes_included:
                self.data = pd.DataFrame(table, columns=['x', 'y', 'z', 'r', 'uni'])
            else:
                self.data = pd.DataFrame(table, columns=['x', 'y', 'z', 'r'])
        self.data['id'] = np.array(self.data.index)
        self.process_data(calculate_dist, calculate_radial_dist, radial_center, verbose=verbose, level=level+1, log_delim=log_delim)
        log_print(f'Done.', verbose, level, log_delim, end_block=True)


    def read_xyzr(self, sss_ov_xyzr, calculate_dist=True, calculate_radial_dist=True, radial_center=[0,0], verbose=True, level=0, log_delim='  '):
        log_print(f'Generating data from Serpent sss_ov_xyzr field', verbose, level, log_delim)
        data = np.reshape(sss_ov_xyzr, (-1, 4))
        self.data = pd.DataFrame(data, columns=['x', 'y', 'z', 'r'])
        self.data['id'] = np.array(self.data.index)
        self.process_data(calculate_dist, calculate_radial_dist, radial_center, verbose=verbose, level=level+1, log_delim=log_delim)
        log_print(f'Done.', verbose, level, log_delim, end_block=True)

    #### PROCESSING RAW DATA ####

    def reset(self, keep_xyzr_if_exists=False, process_data=True, verbose=True, level=0, log_delim='  '):
        log_print(f'Resetting pbed data', verbose, level, log_delim)
        if keep_xyzr_if_exists:
            kept = list(self.data.columns.intersection(['x', 'y', 'z', 'r', 'id', 'uni']))
            log_print(f"Keeping {kept}", verbose, level+1, log_delim)
            if hasattr(self, 'data'):
                self.data = self.data[kept]
            if process_data:
                self.process_data(verbose, level+2, log_delim)
        else:
            self = Pebble_bed(verbose, level+1, log_delim)
        log_print(f'Done.', verbose, level, log_delim, end_block=True)

    def get_positions(self, indices=[0, 1, 2]):
        return self.data[np.array(['x', 'y', 'z'])[indices]]

    def process_data(self, calculate_dist=True, calculate_radial_dist=True, calculate_angles=True, calculate_center=True, calculate_ghosts=False, radial_center=[0,0], verbose=True, level=0, log_delim='  '):
        log_print(f'Processing pbed data', verbose, level, log_delim)
        if calculate_ghosts:
            self.filter_ghosts()
        if calculate_center:
            self.center = np.array(radial_center+[np.nanmean(self.filter_ghosts().data.z)])
        if calculate_dist:
            self.data["dist"] = np.linalg.norm(self.get_positions() - self.center, axis=1)
        if calculate_radial_dist:
            self.data["r_dist"] = np.linalg.norm(self.get_positions([0, 1]) - self.center[:2], axis=1)
        if calculate_angles:
            self.data["azim_angle"] = np.arctan2(-(self.data.y-self.center[0]), -(self.data.x-self.center[1]))*180/np.pi+180 # made to match Serpent angles (y=0 is reference, reversed clock)
        self.N_elements = len(self.data)
        log_print(f'Summary:', verbose, level+1, log_delim)
        log_print(f'Number of elements: {self.N_elements}', verbose, level+2, log_delim)
        log_print(f'Fields: {list(self.data.columns)}', verbose, level+2, log_delim)
        self.radii_list = self.data.r.unique()
        log_print(f'Radii: {self.radii_list}', verbose, level+2, log_delim)
        self.box = np.array([[self.data.x.min(), self.data.x.max()],
                             [self.data.y.min(), self.data.y.max()],
                             [self.data.z.min(), self.data.z.max()]])
        log_print(f'[x,y,z] limits: [{self.box[0,0]:.2f}, {self.box[0,1]:.2f}], [{self.box[1,0]:.2f}, {self.box[1,1]:.2f}],[{self.box[2,0]:.2f}, {self.box[2,1]:.2f}]', verbose, level+2, log_delim)

        if 'uni' in self.data.columns:
            self.universes_list = self.data.uni.unique()
            log_print(f'Universes: {list(self.universes_list)}', verbose, level+2, log_delim)
        log_print(f'Done.', verbose, level, log_delim, end_block=True)

    def filter_ghosts(self, filter=True):
        if filter:
            subpbed = deepcopy(self)
            are_ghosts = self.is_ghost()
            self.N_ghosts = sum(are_ghosts)
            subpbed.data.loc[are_ghosts, :] = np.nan
        return subpbed

    def is_ghost(self, index='all'):
        if isinstance(index, str) and index=='all':
            return np.logical_or(self.data['r'] == 0, np.isnan(self.data['r']))
        else:
            return np.logical_or(self.data.loc[index, 'r'] == 0, np.isnan(self.data.loc[index, 'r']))

    #### PROCESSING DATA ####

    def to_xyzr(self):
        return np.array(self.data[['x','y','z','r']]).flatten()

    def extract(self, fields='all', indices='all'):
        subpbed = deepcopy(self)

        if isinstance(fields, (list, tuple, np.ndarray)) or (isinstance(fields, str) and fields != 'all'):
            subpbed.data = subpbed.data[fields]

        if isinstance(indices, (list, tuple, np.ndarray)) or isinstance(indices, int):
           subpbed.data = subpbed.data.loc[indices]

        return subpbed

    #### POST PROCESSING ####

    def read_results(self, res_file_path):
        if not hasattr(self, 'results'):
            self.results = dict()
        res_reader = serpentTools.read(res_file_path)
        for parameter in res_reader.resdata.keys():
            self.results[parameter] = res_reader.resdata[parameter]

    def extract_res(self, parameter, res_file_path=None, infer_if_incorrect=True, which_rows='all', which_columns='all'):
        if not isinstance(res_file_path, type(None)):
            self.read_results(res_file_path)
        parameter = find_matching(parameter, self.results.keys(), infer_if_incorrect=infer_if_incorrect, lower_ok=True, cutoff=0.4)
        array = self.results[parameter]
        if isinstance(which_rows, str) and which_rows == 'all':
            which_rows = [i for i in range(array.shape[0])]
        array = array[which_rows, :]
        if isinstance(which_columns, str) and which_columns == 'all':
            which_columns = [i for i in range(array.shape[1])]
        array = array[:, which_columns]
        return array

    def extract_time(self, res_file_path=None, which_steps='all', burnup=False):
        if burnup:
            parameter = 'burnup'
        else:
            parameter = 'burnDays'
        return self.extract_res(parameter, res_file_path, which_rows=which_steps, which_columns=0)

    def distribution(self, name, bins=10, ignore_ghosts=False, normalized=False, weights=None, over_max=False, cumulative=False):
        if ignore_ghosts:
            array = self.data.loc[~self.is_ghost(), name]
        else:
            array = self.data[name]

        if isinstance(weights, type(None)):
            weights=np.ones(len(array))
        elif isinstance(weights, (int,float)):
            weights= np.ones(len(array))*weights

        cnt, vals = np.histogram(array, bins=bins, weights=weights)
        if over_max:
            cnt /= np.max(cnt)
        elif normalized:
            cnt /= np.sum(cnt)
        if cumulative:
            cnt = np.cumsum(cnt)
        return cnt, vals

    def profile1D(self, field, fieldx='r_dist', bins=20, error_field=None, absolute_error=False, ignore_ghosts=False, xcut=None, ycut=None, verbose=True, level=0, log_delim='  '):
        profile = pd.Series(dtype=object)
        if isinstance(field, (tuple, list, np.ndarray, pd.DataFrame, pd.Series)):
            log_print(f'1D profile for pbed of user-defined field vs {fieldx}', verbose, level, log_delim)
            self.data['user-defined'] = np.array(field)
        else:
            log_print(f'1D profile for pbed of field "{field}" vs {fieldx}', verbose, level, log_delim)

        if ignore_ghosts:
            data = self.data[~self.is_ghost()]
        else:
            data = self.data
        if not isinstance(xcut, type(None)):
            data = data.loc[data[fieldx].between(xcut[0], xcut[1])]
        if not isinstance(ycut, type(None)):
            data = data.loc[data[field].between(ycut[0], ycut[1])]
        cut_data, bin_limits = pd.cut(data[fieldx], bins=bins, retbins=True)
        grouped_data = data.groupby(cut_data)
        counts = grouped_data.size()
        mean_array = grouped_data[field].mean()
        if not isinstance(error_field, type(None)):
            error_array = grouped_data[error_field].mean()
            if not absolute_error:
                error_array *= mean_array.values
        else:
            error_array = grouped_data[field].std()
        profile['values_array'] = np.array(mean_array.values)
        profile['errors_array'] = np.array(error_array.values)
        profile['data'] = grouped_data
        profile['bins'] = np.array(bin_limits)
        profile['counts_array'] = np.array(counts)
        profile['total_count'] = np.sum(profile['counts_array'])
        profile['weighted_values'] = np.multiply(profile['values_array'], profile['counts_array']/profile['total_count'])

        log_print(f'Done.', verbose, level, log_delim, end_block=True)

        return profile

    def summarize_field(self, names='all', ignore_ghosts=False):
        if ignore_ghosts:
            data = self.data[~self.is_ghost()]
        else:
            data = self.data
        if isinstance(names, str) and names == 'all':
            names = data.columns
        return pd.DataFrame(data[names].describe())

    def add_field(self, name, array, err=False, verbose=True, level=0, log_delim='  '):
        log_print(f'Adding field {name} manually to data (N={len(array)})', verbose, level, log_delim)
        if not hasattr(self, 'fields'):
            self.fields = dict()
        if not err:
            self.fields[name] = array
            if len(array) == len(self.data):
                self.data[name] = self.fields[name]
        else:
            self.fields_err[name] = array
            if len(array) == len(self.data):
                self.data[name+'_err'] = self.fields_err[name]
        log_print(f'Done.', verbose, level, log_delim, end_block=True)

    def read_detector(self, det_file_path, which_dets='all', suffix='', verbose=True, level=0, log_delim='  '):
        log_print(f'Reading det file from {det_file_path}', verbose, level, log_delim)
        det_reader = serpentTools.read(det_file_path, reader='det')
        if not hasattr(self, 'detector_names'):
            self.detector_names = []
        if not hasattr(self, 'fields'):
            self.fields = dict()
        if not hasattr(self, 'fields_err'):
            self.fields_err = dict()
        if not hasattr(self, 'fields_grids'):
            self.fields_grids = dict()

        if which_dets=='all':
            which_dets = list(det_reader.detectors.keys())
        elif isinstance(which_dets, str):
            which_dets = [which_dets]
        for det_name in which_dets:
            tallies = det_reader[det_name].tallies
            errors = det_reader[det_name].errors
            grids = det_reader[det_name].grids
            det_name += suffix

            self.detector_names.append(det_name)
            self.fields[det_name] = tallies
            self.fields_err[det_name] = errors
            self.fields_grids[det_name] = grids

            # Handles up to 1 extra grid level
            tallies = np.atleast_2d(tallies)
            err = np.atleast_2d(self.fields_err[det_name])
            if tallies.shape[1] == len(self.data) and len(grids) <=1:
                if len(grids) == 0:
                    self.data[det_name] = self.fields[det_name]
                    self.data[det_name + '_err'] = self.fields_err[det_name]
                else:
                    name_grid = list(grids)[0]
                    for i_grid in range(grids[name_grid].shape[0]):
                        self.data[f'{det_name}_{name_grid}{i_grid}'] = tallies[i_grid, :]
                        self.data[f'{det_name}_{name_grid}{i_grid}_err'] = err[i_grid, :]
        self.read_files.append(det_file_path)
        log_print(f'Added following detectors: {which_dets}', verbose, level+1, log_delim)
        log_print(f'Done.', verbose, level, log_delim, end_block=True)

    def read_depletion(self, dep_file_path, material_name, dd=False, global_concentration=True, steps='all', fields='all', verbose=True, level=0, log_delim='  '):
        if dd:
            if dep_file_path[-2:] == '.m':
                dep_file_path = dep_file_path[:-2]
            dep_file_path = glob(dep_file_path + '_dd*.m')
        else:
            if dep_file_path[-2:] != '.m':
                dep_file_path += '.m'
            dep_file_path = [dep_file_path]


        log_print(f'Reading dep file from {dep_file_path}', verbose, level, log_delim)
        dep_reader = serpentTools.read(dep_file_path[0], reader='dep')
        for i in range(1, len(dep_file_path)):
            additional_reader = serpentTools.read(dep_file_path[i], reader='dep')
            for mat_name, mat in additional_reader.materials.items():
                dep_reader.materials[mat_name] = mat

        if not hasattr(self, 'fields'):
            self.fields = dict()

        if fields=='all':
            fields = [i for i in list(dep_reader.materials.values())[0].data.keys()]
        if steps=='all':
            steps = [i for i in range(len(list(dep_reader.materials.values())[0].burnup))]
        if material_name=='infer':
            candidates = []
            for name in [mat.name for mat in dep_reader.materials.values()]:
                if 'z' in name and name.split('z')[-1].isdigit():
                    candidates.append(name)
            if len(candidates) == 0:
                raise Exception('No divided material found')
            if len(candidates) > 1:
                raise Exception('More than one divided material found, please choose between: ', candidates)
            else:
                material_name = candidates[0]

        nuc_wise = []
        nuc_names = dep_reader.materials[material_name].names

        if global_concentration:
            for field_name in fields:
                if isinstance(dep_reader.materials[material_name].data[field_name][0], (tuple, list, np.ndarray)):
                    for i_name, name in enumerate(nuc_names):
                        self.fields[f'{name}_{field_name}'] = dep_reader.materials[material_name].data[field_name][i_name, :]
                else:
                    self.fields[f'{field_name}'] = dep_reader.materials[material_name].data[field_name]

        else:
            for step in steps:
                for field_name in fields:
                    if isinstance(dep_reader.materials[material_name].data[field_name][0], (tuple, list, np.ndarray)):
                        for name in nuc_names:
                            self.fields[f'{name}_{field_name}_{step}'] = np.empty(len(self.data))
                        nuc_wise.append(True)
                    else:
                        self.fields[f'{field_name}_{step}'] = np.empty(len(self.data))
                        nuc_wise.append(False)
            for step in steps:
                for i in range(1, len(self.data)+1):
                    for i_field, field_name in enumerate(fields):
                        if nuc_wise[i_field]:
                            for i_name, name in enumerate(nuc_names):
                                self.fields[f'{name}_{field_name}_{step}'][i-1] = dep_reader.materials[f'{material_name}z{i}'].data[field_name][step][i_name]
                        else:
                            self.fields[f'{field_name}_{step}'][i-1] = dep_reader.materials[f'{material_name}z{i}'].data[field_name][step]

                for i_field, field_name in enumerate(fields):
                    if nuc_wise[i_field]:
                        for name in nuc_names:
                            self.data[f'{name}_{field_name}_{step}'] = self.fields[f'{name}_{field_name}_{step}']
                    else:
                        self.data[f'{field_name}_{step}'] = self.fields[f'{field_name}_{step}']

            for file in dep_file_path:
                self.read_files.append(file)

        log_print(f'Added following fields for step(s) {steps}: {fields}', verbose, level+1, log_delim)
        log_print(f'Done.', verbose, level, log_delim, end_block=True)

    #### DISCARDED ####

    def read_discarded_files(self, discard_folder, list_indices='all', ignore_init=False, verbose=True, level=0, log_delim='  '):
        discard_files = glob(f'{discard_folder}/*')
        prefix = os.path.commonprefix(discard_files)
        if isinstance(list_indices, str) and list_indices == 'all':
            list_indices = np.sort([int(file.split(prefix)[-1]) for file in discard_files])

        log_print(f'Reading {len(list_indices)} discarded pebbles files from {discard_folder}.', verbose, level, log_delim)
        self.discarded = pd.DataFrame()
        for index in list_indices:
            file = f'{prefix}{index}'
            self.discarded = pd.concat((self.discarded, self.read_discarded_file(file, ignore_init, verbose=verbose, level=level+1, log_delim=log_delim)))
        log_print(f'Found {len(self.discarded)} discarded pebbles.', verbose, level+1, log_delim)
        log_print(f'Done.', verbose, level, log_delim, end_block=True)

    def read_discarded_file(self, discard_file_path, ignore_init=False, verbose=True, level=0, log_delim='  '):
        log_print(f'Reading discarded pebbles file {discard_file_path}.', verbose, level, log_delim)
        discarded = pd.read_csv(discard_file_path)
        if ignore_init:
            discarded = discarded.loc[~discarded.init_pebble.astype(bool)]
        log_print(f'Done.', verbose, level, log_delim, end_block=True)
        return discarded

    def group_discarded(self, field_groupby='discarded_step', which_fields='all', verbose=True, level=0, log_delim='  '):
        log_print(f'Grouping discarded data by {field_groupby} with following fields: {which_fields}', verbose, level, log_delim)
        if isinstance(which_fields, str) and which_fields == 'all':
            which_fields = self.discarded.columns
        elif isinstance(which_fields, str):
            which_fields = [which_fields]
        groups = self.discarded.groupby(field_groupby)[which_fields]
        log_print(f'Found {len(groups)} groups', verbose, level, log_delim)
        log_print(f'Done.', verbose, level, log_delim, end_block=True)
        return groups

    #### SORTING DATA ####

    def slice(self, dir_id='x', val='middle', *, tol=None, verbose=True, level=0, log_delim='  '):
        dir_id = str(dir_id)
        # if dir_id.lower() in ["x", "y", "z"]:
        #     dir_id = ["x", "y", "z"].index(dir_id.lower())
        # dir_id = int(dir_id)

        if isinstance(val, str) and val == 'middle':
            val = np.nanmean(self.filter_ghosts().data[dir_id])

        log_print(f'Slicing pbed in direction {dir_id} at value {val:.4E}', verbose, level, log_delim)
        if isinstance(tol, type(None)):
            log_print('Tolerance set to radii values', verbose, level=level+1, log_delim=log_delim)
        else:
            log_print(f'Tolerance set to {tol}', verbose, level=level+1, log_delim=log_delim)

        sub_pbed = deepcopy(self)

        if isinstance(tol, type(None)):
            sub_pbed.data = self.data[np.abs(self.data[dir_id] - val) <= self.data["r"]]
        else:
            sub_pbed.data =  self.data[np.abs(self.data[dir_id] - val) <= tol]
        sub_pbed.N_elements

        log_print(f'Done.', verbose, level, log_delim, end_block=True)
        return sub_pbed


    def clip(self, dir_id=0, val='middle', direction=+1, verbose=True, level=0, log_delim='  ', logical_and=True):
        sdir = '>=' if direction==1 else '<='
        sub_pbed = deepcopy(self)
        if isinstance(dir_id, (list, np.ndarray, tuple)):
            dir_id_list = np.array(dir_id)
        else:
            dir_id_list = np.array([dir_id])
        if not isinstance(val, (list, np.ndarray, tuple)):
            val_list = np.array([val for _ in range(len(dir_id_list))])
        else:
            val_list = np.array(val)
        if not isinstance(direction, (list, np.ndarray, tuple)):
            direction_list = np.array([direction for _ in range(len(dir_id_list))])
        else:
            direction_list = np.array(direction)

        conditions = []
        for i, dir_id in enumerate(dir_id_list):
            dir_id = str(dir_id)
            val = val_list[i]
            direction = direction_list[i]
            if dir_id.lower() in ["x", "y", "z"]:
                if isinstance(val, str) and val == 'middle':
                    val = np.nanmean(self.filter_ghosts().data[dir_id])
                conditions.append(list(self.data[dir_id]*direction >= val*direction))
            elif dir_id == '4' or dir_id.lower() in ['dist', 'd']:
                if isinstance(val, str) and val == 'middle':
                    val = np.nanmean(self.filter_ghosts().data.dist)
                conditions.append(list(self.data.dist*direction >= val*direction))
            elif dir_id == '5' or dir_id.lower() in ['r_dist', 'rdist', 'r']:
                if isinstance(val, str) and val == 'middle':
                    val = np.nanmean(self.filter_ghosts().data.r_dist)
                conditions.append(list(self.data.r_dist*direction >= val*direction))
            else:
                if isinstance(val, str) and val == 'middle':
                    val = np.nanmean(self.filter_ghosts().data[dir_id])
                conditions.append(self.data[dir_id]*direction >= val*direction)
        if len(conditions) == 1:
            gathered_conditions = np.array(conditions[0])
        else:
            if logical_and:
                gathered_conditions = np.prod((np.array(conditions)), axis=0).astype(bool)
            else:
                gathered_conditions = np.sum((np.array(conditions)), axis=0).astype(bool)
        sub_pbed.data =  self.data[gathered_conditions]
        log_print(f'Done.', verbose, level, log_delim, end_block=True)
        return sub_pbed

    #### PLOTTING ####

    def projection(self, dir_id, val, verbose=True, level=0, log_delim='  '):
        dir_id = str(dir_id)

        log_print(f'Projecting pbed in direction {dir_id} at value {val:.4E}', verbose, level, log_delim)
        rel_pos = val - self.data[dir_id]
        tmp = self.data.r**2 - rel_pos**2
        r_projected = np.ones_like(tmp)*np.nan
        r_projected[tmp >= 0] = tmp[tmp >= 0]**0.5

        log_print(f'Done.', verbose, level, log_delim, end_block=True)
        return r_projected

    def plot_global_concentration(self, parameter, which_steps='all', dep_file_path=None, material_name=None, use_time=True, res_file_path=None, ylabel=None, plot_labels=None, new_fig=True, xlim=None, ylim=None, equal=False, fig_size=None, save_fig=False, fig_folder='./', fig_name=None, fig_suffix='', fig_dpi=400, plot_title=None, verbose=True, level=0, log_delim='  '):
        if not hasattr(self, 'fields'):
            self.fields = dict()

        tmp = find_matching(parameter, self.fields.keys())
        if isinstance(tmp, type(None)):
            if isinstance(dep_file_path, type(None)):
                raise Exception(f'{parameter} not found in Pebble Bed fields, please put a dep_file_path')
            self.read_depletion(dep_file_path, material_name, global_concentration=True)
            tmp = find_matching(parameter, self.fields.keys())

        array = self.fields[tmp]

        if use_time:
            x = self.extract_res('burnDays', res_file_path, which_rows=which_steps, which_columns=0)
            xlabel = 'Time [EPFD]'
        else:
            x = self.extract_res('burnStep', res_file_path, which_rows=which_steps, which_columns=0)
            xlabel = 'Step #'
        if isinstance(which_steps, str) and which_steps=='all':
            which_steps = np.arange(len(array))[:-1]
        array = array[which_steps]
        if new_fig:
            if not isinstance(fig_size, type(None)):
                plt.figure(figsize=fig_size)
            else:
                plt.figure()

        ax = plt.gca()
        plt.step(x, array, label=plot_labels)
        plt.xlabel(xlabel)
        if isinstance(ylabel, type(None)):
            plt.ylabel(parameter)
        else:
            plt.ylabel(ylabel)
        plt.legend()

        if not isinstance(xlim, type(None)):
            ax.set_xlim(xlim)
        if not isinstance(ylim, type(None)):
            ax.set_ylim(ylim)
        if isinstance(plot_title, type(None)):
            plt.title(f"Evolution of {parameter}")
        else:
            plt.title(plot_title)
        ax.autoscale_view()
        if equal:
            ax.set_aspect("equal", adjustable="box")
        plt.grid(visible=True)
        ax.set_axisbelow(True)

        plt.tight_layout()

        if save_fig:
            if isinstance(fig_name, type(None)):
                fig_path = f'{fig_folder}/global_conc_plot_{parameter}{fig_suffix}.png'
            else:
                fig_path = f'{fig_folder}/{fig_name}{fig_suffix}.png'
            log_print(f'Saving figure in {fig_path}', verbose, level+1, log_delim)
            plt.savefig(fig_path, dpi=fig_dpi, bbox_inches='tight')

        log_print(f'Done.', verbose, level, log_delim, end_block=True)
        return ax


    def plot_res(self, parameter, infer_if_incorrect=True, which_rows='all', which_columns='all', which_columns_errors=None, res_file_path=None, use_time=True, ylabel=None, plot_labels=None, new_fig=True, xlim=None, ylim=None, equal=False, fig_size=None, save_fig=False, fig_folder='./', fig_name=None, fig_suffix='', fig_dpi=400, plot_title=None, verbose=True, level=0, log_delim='  '):
        if not hasattr(self, 'results'):
            self.results = dict()

        tmp = find_matching(parameter, self.results.keys())
        if isinstance(tmp, type(None)):
            if isinstance(res_file_path, type(None)):
                raise Exception(f'{parameter} not found in Pebble Bed results, please put a res_file_path')
            array = self.extract_res(parameter, res_file_path, infer_if_incorrect, which_rows, which_columns)
        else:
            array = self.results[tmp]
        if use_time:
            x = self.extract_res('burnDays', res_file_path, which_rows=which_rows, which_columns=0)
            xlabel = 'Time [EPFD]'
        else:
            x = self.extract_res('burnStep', res_file_path, which_rows=which_rows, which_columns=0)
            xlabel = 'Step #'

        if new_fig:
            if not isinstance(fig_size, type(None)):
                plt.figure(figsize=fig_size)
            else:
                plt.figure()

        ax = plt.gca()
        if isinstance(which_columns_errors, type(None)):
            plt.plot(x, array, label=plot_labels)
        else:
            array_err = self.extract_res(parameter, res_file_path, infer_if_incorrect, which_rows=which_rows, which_columns=which_columns_errors)
            errorbar(x, array, yerr=array_err, label=plot_labels, steps=True)
        plt.xlabel(xlabel)
        if isinstance(ylabel, type(None)):
            plt.ylabel(parameter)
        else:
            plt.ylabel(ylabel)
        plt.legend()

        if not isinstance(xlim, type(None)):
            ax.set_xlim(xlim)
        if not isinstance(ylim, type(None)):
            ax.set_ylim(ylim)
        if isinstance(plot_title, type(None)):
            plt.title(f"Evolution of {parameter}")
        else:
            plt.title(plot_title)
        ax.autoscale_view()
        if equal:
            ax.set_aspect("equal", adjustable="box")
        plt.grid(visible=True)
        ax.set_axisbelow(True)

        plt.tight_layout()

        if save_fig:
            if isinstance(fig_name, type(None)):
                fig_path = f'{fig_folder}/res_plot_{parameter}{fig_suffix}.png'
            else:
                fig_path = f'{fig_folder}/{fig_name}{fig_suffix}.png'
            log_print(f'Saving figure in {fig_path}', verbose, level+1, log_delim)
            plt.savefig(fig_path, dpi=fig_dpi, bbox_inches='tight')

        log_print(f'Done.', verbose, level, log_delim, end_block=True)
        return ax

    def plot_distribution(self, parameter, bins=10, show_mean=False, show_quartiles=False, ignore_ghosts=False, normalized=False, norm_pc=True, over_max=False, fancy=True, weights=None, cumulative=False, xlabel=None, ylabel=None, plot_labels=None, alpha=0.6, new_fig=True, xlim=None, ylim=None, equal=False, lab=None, colormap='turbo', fig_size=None, save_fig=False, fig_folder='./', fig_name=None, fig_suffix='', fig_dpi=400, plot_title=None, linewidth=2, color=None, verbose=True, level=0, log_delim='  '):

        if new_fig:
            if not isinstance(fig_size, type(None)):
                plt.figure(figsize=fig_size)
            else:
                plt.figure()

        ax = plt.gca()
        cnt, vals = self.distribution(parameter, bins=bins, ignore_ghosts=ignore_ghosts, normalized=normalized, weights=weights, cumulative=cumulative, over_max=over_max)
        if normalized and norm_pc:
            cnt*=100
        if cumulative:
            if fancy:
                slopes = np.divide(np.diff(np.concatenate(([0],cnt))), 1) # np.diff(vals))
                cmap = cm.get_cmap(colormap)
                norm = Normalize(vmin = np.min(slopes), vmax = np.max(slopes))
                normalized_values = norm(slopes)
                colors = cmap(normalized_values)
                sm = cm.ScalarMappable(cmap=cmap, norm=norm)
                sm.set_array([])
                plt.colorbar(sm, label=lab)
            else:
                colors=color
            if over_max:
                lab='Distribution [a.u]'
            elif normalized:
                lab='Fraction of the core [%]'
            else:
                lab='Number of pebbles'

        else:
            colors=color
        if fancy:
            plt.bar(vals[:-1], cnt, width=np.diff(vals), color=colors, alpha=alpha, label=plot_labels)
        else:
            plt.plot(vals, np.concatenate(([cnt[0]], cnt)), drawstyle='steps-pre', color=color, label=plot_labels, linewidth=linewidth)
        if show_mean:
            plt.axvline(float(self.data[parameter].mean()), linestyle='--', color='k', label='Average', linewidth=linewidth)
        if show_quartiles:
            summary = self.summarize_field(parameter, ignore_ghosts=ignore_ghosts)
            plt.axvline(float(summary.loc['25%']), linestyle=':', color='r', label='Quartiles', linewidth=linewidth)
            plt.axvline(float(summary.loc['50%']), linestyle=':', color='r', linewidth=linewidth)
            plt.axvline(float(summary.loc['75%']), linestyle=':', color='r', linewidth=linewidth)
        if isinstance(xlabel, type(None)):
            plt.xlabel(parameter)
        else:
            plt.xlabel(xlabel)
        if isinstance(ylabel, type(None)):
            ylabel = 'Number of pebbles'
        if cumulative:
            ylabel = f'Cumulative number of pebbles'
        if over_max:
            ylabel ='Distribution [a.u]'
        if normalized:
            ylabel += ' [%]'

        plt.ylabel(ylabel)

        if not isinstance(xlim, type(None)):
            ax.set_xlim(xlim)
        if not isinstance(ylim, type(None)):
            ax.set_ylim(ylim)
        if isinstance(plot_title, type(None)):
            plt.title(f"Distribution of {parameter}")
        else:
            plt.title(plot_title)
        ax.autoscale_view()
        if equal:
            ax.set_aspect("equal", adjustable="box")
        plt.grid(visible=True)
        ax.set_axisbelow(True)
        ax.ticklabel_format(axis='x', style='', scilimits=(-2,3))
        plt.tight_layout()

        if save_fig:
            if isinstance(fig_name, type(None)):
                fig_path = f'{fig_folder}/dist_plot_{parameter}{fig_suffix}.png'
            else:
                fig_path = f'{fig_folder}/{fig_name}{fig_suffix}.png'
            log_print(f'Saving figure in {fig_path}', verbose, level+1, log_delim)
            plt.savefig(fig_path, dpi=fig_dpi, bbox_inches='tight')

        log_print(f'Done.', verbose, level, log_delim, end_block=True)
        return ax

    def plot1D_sum(self, field, fieldx='r_dist', bins=20, error_field=None, absolute_error=False, cumulative=False, ignore_ghosts=False, xlim=None, ylim=None, xlabel=None, ylabel=None, label=None, color=None, xcut=None, ycut=None, equal=False, field_title=None, fig_size=None, new_fig=True, save_fig=False, fig_folder='./', fig_name=None, fig_suffix='', fig_dpi=400, plot_title=None, verbose=True, level=0, log_delim='  '):
        if isinstance(field, (tuple, list, np.ndarray, pd.DataFrame, pd.Series)):
            log_print(f'1D plotting pbed of user-defined field vs {fieldx}', verbose, level, log_delim)
            if isinstance(field_title, type(None)):
                self.data['user-defined'] = np.array(field)
                field = 'user-defined'
                field_title = field
            else:
                self.data[field_title] = np.array(field)
                field = field_title

        else:
            if isinstance(field_title, type(None)):
                field_title = field
            log_print(f'1D plotting pbed of field "{field}" vs {fieldx}', verbose, level, log_delim)


        if ignore_ghosts:
            data = self.data[~self.is_ghost()]
        else:
            data = self.data
        if not isinstance(xcut, type(None)):
            data = data.loc[data[fieldx].between(xcut[0], xcut[1])]
        if not isinstance(ycut, type(None)):
            data = data.loc[data[field].between(ycut[0], ycut[1])]
        cut_data, bin_limits = pd.cut(data[fieldx], bins=bins, retbins=True)
        bin_middle = (bin_limits[1:]+bin_limits[:-1])/2
        grouped_data = data.groupby(cut_data)
        sum_array = grouped_data[field].sum()
        if cumulative:
            sum_array = pd.Series(np.cumsum(list(sum_array)))
        if new_fig:
            if not isinstance(fig_size, type(None)):
                plt.figure(figsize=fig_size)
            else:
                plt.figure()

        ax = plt.gca()
        if not isinstance(error_field, type(None)):
            err_array = [np.sqrt(np.sum(np.power(i, 2))) for i in np.array(grouped_data[error_field].apply(list))]
            if not absolute_error:
                err_array *= sum_array.values
            errorbar(bin_limits, np.concatenate(([sum_array.iloc[0]], sum_array)), np.concatenate(([err_array[0]], err_array)), steps=True, color=color, label=label)
        else:
            plt.plot(bin_limits, np.concatenate(([sum_array.iloc[0]], sum_array)), drawstyle='steps-pre', color=color, label=label)

        if not isinstance(xlabel, type(None)):
            plt.xlabel(xlabel)
        else:
            plt.xlabel(fieldx)
        if not isinstance(ylabel, type(None)):
            plt.ylabel(ylabel)
        else:
            plt.ylabel(field)
        if not isinstance(xlim, type(None)):
            ax.set_xlim(xlim)
        if not isinstance(ylim, type(None)):
            ax.set_ylim(ylim)
        if isinstance(plot_title, type(None)):
            plt.title(f"{field} vs {fieldx}")
        else:
            plt.title(plot_title)
        ax.autoscale_view()
        if equal:
            ax.set_aspect("equal", adjustable="box")
        plt.grid(visible=True)
        ax.set_axisbelow(True)

        plt.tight_layout()

        if save_fig:
            if isinstance(fig_name, type(None)):
                fig_path = f'{fig_folder}/1D_sum_plot_{field}_vs_{fieldx}{fig_suffix}.png'
            else:
                fig_path = f'{fig_folder}/{fig_name}{fig_suffix}.png'
            log_print(f'Saving figure in {fig_path}', verbose, level+1, log_delim)
            plt.savefig(fig_path, dpi=fig_dpi, bbox_inches='tight')

        log_print(f'Done.', verbose, level, log_delim, end_block=True)
        return ax

    def plot1D_slice(self, field, fieldx='r_dist', dir_id='z', val='middle', tol='r', bins=20, error_field=None, absolute_error=False, plot_std=True, plot_extremes=False, linewidth=1, ignore_ghosts=False, steps=True, xlim=None, ylim=None, xlabel=None, ylabel=None, xcut=None, ycut=None, equal=False, field_title=None, fig_size=None, new_fig=True, label=None, save_fig=False, fig_folder='./', alpha_error=0.3, color=None, fig_name=None, fig_suffix='', fig_dpi=400, plot_title=None, linestyle='-', verbose=True, level=0, log_delim='  '):
        if isinstance(val, str) and val == 'middle':
            val = np.nanmean(self.filter_ghosts().data[dir_id])

        if isinstance(field, (tuple, list, np.ndarray, pd.DataFrame, pd.Series)):
            log_print(f'1D plotting pbed of user-defined field vs {fieldx} in direction={dir_id} and values={val}', verbose, level, log_delim)
            if isinstance(field_title, type(None)):
                self.data['user-defined'] = np.array(field)
                field = 'user-defined'
                field_title = field
            else:
                self.data[field_title] = np.array(field)
                field = field_title

        else:
            if isinstance(field_title, type(None)):
                field_title = field
            log_print(f'1D plotting pbed of field "{field}" vs {fieldx} in direction={dir_id} and values={val}', verbose, level, log_delim)

        if isinstance(tol, str) and tol=='r':
            tol = self.data.r
        data = self.data.loc[np.abs(self.data[dir_id]-val)<=tol]
        if ignore_ghosts:
            data = data[~self.is_ghost()]
        if not isinstance(xcut, type(None)):
            data = data.loc[data[fieldx].between(xcut[0], xcut[1])]
        if not isinstance(ycut, type(None)):
            data = data.loc[data[field].between(ycut[0], ycut[1])]
        cut_data, bin_limits = pd.cut(data[fieldx], bins=bins, retbins=True)
        bin_middle = (bin_limits[1:]+bin_limits[:-1])/2
        grouped_data = data.groupby(cut_data)
        mean_array = grouped_data[field].mean()

        if new_fig:
            if not isinstance(fig_size, type(None)):
                plt.figure(figsize=fig_size)
            else:
                plt.figure()

        ax = plt.gca()
        if not isinstance(error_field, type(None)):
            mean_err_array = grouped_data[error_field].mean()
            if not absolute_error:
                mean_err_array *= mean_array.values
            if steps:
                x = bin_limits
                y = np.concatenate(([mean_array.iloc[0]], mean_array))
                yerr = np.concatenate(([mean_err_array.iloc[0]], mean_err_array))
            else:
                x = bin_middle
                y = mean_array
                yerr = mean_err_array
            errorbar(x, y, yerr, steps=steps, label=label, linewidth=linewidth, alpha=alpha_error, color=color, linestyle=linestyle)
        elif plot_std and not plot_extremes:
            std_array = grouped_data[field].std()
            if steps:
                x = bin_limits
                y = np.concatenate(([mean_array.iloc[0]], mean_array))
                yerr = np.concatenate(([mean_err_array.iloc[0]], mean_err_array))
            else:
                x = bin_middle
                y = mean_array
                yerr = std_array
            errorbar(x, y, yerr, steps=steps, label=label, linewidth=linewidth, alpha=alpha_error, color=color, linestyle=linestyle)
        elif plot_extremes:
            extreme_array = np.vstack((grouped_data[field].min(), grouped_data[field].max()))
            y = np.concatenate(([mean_array.iloc[0]], mean_array))
            yerr = np.hstack((np.atleast_2d(extreme_array[:,0]).T, extreme_array)).T
            if steps:
                x = bin_limits
                y = np.concatenate(([mean_array.iloc[0]], mean_array))
                yerr = np.hstack((np.atleast_2d(extreme_array[:,0]).T, extreme_array)).T
            else:
                x = bin_middle
                y = mean_array
                yerr = extreme_array
            yerr[:, 0] = y-yerr[:,0]
            yerr[:, 1] = yerr[:,1] - y
            errorbar(x, y, yerr, steps=steps, label=label, linewidth=linewidth, color=color, linestyle=linestyle)
        else:
            if steps:
                plt.plot(bin_limits, np.concatenate(([mean_array.iloc[0]], mean_array)), drawstyle='steps-pre', label=label, linewidth=linewidth, color=color, linestyle=linestyle)
            else:
                plt.plot(bin_middle, mean_array, label=label, linewidth=linewidth, color=color, linestyle=linestyle)

        if not isinstance(xlabel, type(None)):
            plt.xlabel(xlabel)
        else:
            plt.xlabel(fieldx)
        if not isinstance(ylabel, type(None)):
            plt.ylabel(ylabel)
        else:
            plt.ylabel(field)
        if not isinstance(xlim, type(None)):
            ax.set_xlim(xlim)
        if not isinstance(ylim, type(None)):
            ax.set_ylim(ylim)
        if isinstance(plot_title, type(None)):
            plt.title(f"{field} vs {fieldx}")
        else:
            plt.title(plot_title)
        ax.autoscale_view()
        if equal:
            ax.set_aspect("equal", adjustable="box")

        plt.grid(visible=True)
        ax.set_axisbelow(True)
        plt.tight_layout()

        if save_fig:
            if isinstance(fig_name, type(None)):
                fig_path = f'{fig_folder}/1D_mean_plot_{field}_vs_{fieldx}{fig_suffix}.png'
            else:
                fig_path = f'{fig_folder}/{fig_name}{fig_suffix}.png'
            log_print(f'Saving figure in {fig_path}', verbose, level+1, log_delim)
            plt.savefig(fig_path, dpi=fig_dpi, bbox_inches='tight')

        log_print(f'Done.', verbose, level, log_delim, end_block=True)
        return ax

    def plot1D_mean(self, field, fieldx='r_dist', bins=20, error_field=None, absolute_error=False, plot_std=True, plot_extremes=False, labelize_error=True, linewidth=1, ignore_ghosts=False, steps=True, xlim=None, ylim=None, xlabel=None, ylabel=None, xcut=None, ycut=None, equal=False, field_title=None, fig_size=None, new_fig=True, label=None, save_fig=False, fig_folder='./', alpha_error=0.3, fig_name=None, fig_suffix='', fig_dpi=400, plot_title=None, linestyle='-', color=None, verbose=True, level=0, log_delim='  '):
        if isinstance(field, (tuple, list, np.ndarray, pd.DataFrame, pd.Series)):
            log_print(f'1D plotting pbed of user-defined field vs {fieldx}', verbose, level, log_delim)
            if isinstance(field_title, type(None)):
                self.data['user-defined'] = np.array(field)
                field = 'user-defined'
                field_title = field
            else:
                self.data[field_title] = np.array(field)
                field = field_title

        else:
            if isinstance(field_title, type(None)):
                field_title = field
            log_print(f'1D plotting pbed of field "{field}" vs {fieldx}', verbose, level, log_delim)


        if ignore_ghosts:
            data = self.data[~self.is_ghost()]
        else:
            data = self.data
        if not isinstance(xcut, type(None)):
            data = data.loc[data[fieldx].between(xcut[0], xcut[1])]
        if not isinstance(ycut, type(None)):
            data = data.loc[data[field].between(ycut[0], ycut[1])]
        cut_data, bin_limits = pd.cut(data[fieldx], bins=bins, retbins=True)
        bin_middle = (bin_limits[1:]+bin_limits[:-1])/2
        grouped_data = data.groupby(cut_data)
        mean_array = grouped_data[field].mean()

        if new_fig:
            if not isinstance(fig_size, type(None)):
                plt.figure(figsize=fig_size)
            else:
                plt.figure()

        ax = plt.gca()
        if not isinstance(error_field, type(None)):
            mean_err_array = grouped_data[error_field].mean()
            if not absolute_error:
                mean_err_array *= mean_array.values
            if steps:
                x = bin_limits
                y = np.concatenate(([mean_array.iloc[0]], mean_array))
                yerr = np.concatenate(([mean_err_array.iloc[0]], mean_err_array))
            else:
                x = bin_middle
                y = mean_array
                yerr = mean_err_array
            errorbar(x, y, yerr, steps=steps, label=label, linewidth=linewidth, linestyle=linestyle, alpha=alpha_error, color=color, labelize_error=labelize_error)
        elif plot_std and not plot_extremes:
            std_array = grouped_data[field].std()
            if steps:
                x = bin_limits
                y = np.concatenate(([mean_array.iloc[0]], mean_array))
                yerr = np.concatenate(([std_array.iloc[0]], std_array))
            else:
                x = bin_middle
                y = mean_array
                yerr = std_array
            errorbar(x, y, yerr, steps=steps, label=label, linewidth=linewidth, linestyle=linestyle, alpha=alpha_error, color=color, labelize_error=labelize_error)
        elif plot_extremes:
            extreme_array = np.vstack((grouped_data[field].min(), grouped_data[field].max()))
            if steps:
                x = bin_limits
                y = np.concatenate(([mean_array.iloc[0]], mean_array))
                yerr = np.hstack((np.atleast_2d(extreme_array[:,0]).T, extreme_array)).T
            else:
                x = bin_middle
                y = mean_array
                yerr = extreme_array.T
            yerr[:, 0] = y-yerr[:,0]
            yerr[:, 1] = yerr[:,1] - y
            errorbar(x, y, yerr, steps=steps, label=label, linewidth=linewidth, color=color, linestyle=linestyle, labelize_error=labelize_error)
        else:
            if steps:
                plt.plot(bin_limits, np.concatenate(([mean_array.iloc[0]], mean_array)), drawstyle='steps-pre', label=label, linewidth=linewidth, color=color, linestyle=linestyle)
            else:
                plt.plot(bin_middle, mean_array, label=label, linewidth=linewidth, color=color, linestyle=linestyle)

        if not isinstance(xlabel, type(None)):
            plt.xlabel(xlabel)
        else:
            plt.xlabel(fieldx)
        if not isinstance(ylabel, type(None)):
            plt.ylabel(ylabel)
        else:
            plt.ylabel(field)
        if not isinstance(xlim, type(None)):
            ax.set_xlim(xlim)
        if not isinstance(ylim, type(None)):
            ax.set_ylim(ylim)
        if isinstance(plot_title, type(None)):
            plt.title(f"{field} vs {fieldx}")
        else:
            plt.title(plot_title)
        ax.autoscale_view()
        if equal:
            ax.set_aspect("equal", adjustable="box")

        plt.grid(visible=True)
        ax.set_axisbelow(True)
        plt.tight_layout()

        if save_fig:
            if isinstance(fig_name, type(None)):
                fig_path = f'{fig_folder}/1D_mean_plot_{field}_vs_{fieldx}{fig_suffix}.png'
            else:
                fig_path = f'{fig_folder}/{fig_name}{fig_suffix}.png'
            log_print(f'Saving figure in {fig_path}', verbose, level+1, log_delim)
            plt.savefig(fig_path, dpi=fig_dpi, bbox_inches='tight')

        log_print(f'Done.', verbose, level, log_delim, end_block=True)
        return ax


    def plot1D_extreme(self, field, fieldx='r_dist', bins=20, what='max', linewidth=1, ignore_ghosts=False, steps=True, xlim=None, ylim=None, xlabel=None, ylabel=None, xcut=None, ycut=None, equal=False, field_title=None, fig_size=None, new_fig=True, label=None, save_fig=False, fig_folder='./', alpha_error=0.3, fig_name=None, fig_suffix='', fig_dpi=400, plot_title=None, linestyle='-', color=None, verbose=True, level=0, log_delim='  '):
        if isinstance(field, (tuple, list, np.ndarray, pd.DataFrame, pd.Series)):
            log_print(f'1D plotting pbed of user-defined field vs {fieldx}', verbose, level, log_delim)
            if isinstance(field_title, type(None)):
                self.data['user-defined'] = np.array(field)
                field = 'user-defined'
                field_title = field
            else:
                self.data[field_title] = np.array(field)
                field = field_title

        else:
            if isinstance(field_title, type(None)):
                field_title = field
            log_print(f'1D plotting pbed of field "{field}" vs {fieldx}', verbose, level, log_delim)


        if ignore_ghosts:
            data = self.data[~self.is_ghost()]
        else:
            data = self.data
        if not isinstance(xcut, type(None)):
            data = data.loc[data[fieldx].between(xcut[0], xcut[1])]
        if not isinstance(ycut, type(None)):
            data = data.loc[data[field].between(ycut[0], ycut[1])]
        cut_data, bin_limits = pd.cut(data[fieldx], bins=bins, retbins=True)
        bin_middle = (bin_limits[1:]+bin_limits[:-1])/2
        grouped_data = data.groupby(cut_data)
        if what=='max':
            extreme_array = grouped_data[field].max()
        else:
            extreme_array = grouped_data[field].min()
        if new_fig:
            if not isinstance(fig_size, type(None)):
                plt.figure(figsize=fig_size)
            else:
                plt.figure()

        ax = plt.gca()
        if steps:
            x = bin_limits
            y = np.concatenate(([extreme_array.iloc[0]], extreme_array))
            plt.plot(x, y, drawstyle='steps-pre', label=label, linewidth=linewidth, color=color, linestyle=linestyle)
        else:
            x = bin_middle
            y = extreme_array
            plt.plot(x, extreme_array, label=label, linewidth=linewidth, color=color, linestyle=linestyle)


        if not isinstance(xlabel, type(None)):
            plt.xlabel(xlabel)
        else:
            plt.xlabel(fieldx)
        if not isinstance(ylabel, type(None)):
            plt.ylabel(ylabel)
        else:
            plt.ylabel(field)
        if not isinstance(xlim, type(None)):
            ax.set_xlim(xlim)
        if not isinstance(ylim, type(None)):
            ax.set_ylim(ylim)
        if isinstance(plot_title, type(None)):
            plt.title(f"{field} vs {fieldx}")
        else:
            plt.title(plot_title)
        ax.autoscale_view()
        if equal:
            ax.set_aspect("equal", adjustable="box")

        plt.grid(visible=True)
        ax.set_axisbelow(True)
        plt.tight_layout()

        if save_fig:
            if isinstance(fig_name, type(None)):
                fig_path = f'{fig_folder}/1D_mean_plot_{field}_vs_{fieldx}{fig_suffix}.png'
            else:
                fig_path = f'{fig_folder}/{fig_name}{fig_suffix}.png'
            log_print(f'Saving figure in {fig_path}', verbose, level+1, log_delim)
            plt.savefig(fig_path, dpi=fig_dpi, bbox_inches='tight')

        log_print(f'Done.', verbose, level, log_delim, end_block=True)
        return ax


    def plot2D(self, field='id', dir_id=0, val='middle', colormap='turbo', xlim=None, ylim=None, tol=None, equal=True, field_title=None, clim=None, log_color=False, vertical_cbar=False, shrink_cbar=1, pad_cbar=0.15, translation=None, no_ax=False, superimpose_Serpent=False, Serpent_xsize=None, Serpent_ysize=None, Serpent_geom_path=None, fig_size=None, new_fig=True, save_fig=False, fig_folder='./', fig_name=None, fig_suffix='', fig_dpi=400, plot_title=None, patch_clip=None, unit='cm', verbose=True, level=0, log_delim='  '):
        dir_id = str(dir_id)
        if dir_id == '0' or dir_id.lower() == 'x':
            xdir = 'y'
            ydir = 'z'
            dir_id = 'x'
        elif dir_id == '1' or dir_id.lower() == 'y':
            xdir = 'x'
            ydir = 'z'
            dir_id = 'y'
        elif dir_id == '2' or dir_id.lower() == 'z':
            xdir = 'x'
            ydir = 'y'
            dir_id = 'z'

        if isinstance(val, str) and val == 'middle':
            val = np.nanmean(self.filter_ghosts().data[dir_id])

        if isinstance(field, (tuple, list, np.ndarray, pd.DataFrame, pd.Series)):
            log_print(f'2D plotting pbed in direction {dir_id} at value {val:.4E}, showing user-defined array', verbose, level, log_delim)
            if isinstance(field_title, type(None)):
                self.data['user-defined'] = np.array(field)
                field = 'user-defined'
                field_title = field
            else:
                self.data[field_title] = np.array(field)
                field = field_title

        else:
            if isinstance(field_title, type(None)):
                field_title = field
            log_print(f'2D plotting pbed in direction {dir_id} at value {val:.4E}, showing field {field}', verbose, level, log_delim)

        # dir_id = int(dir_id)
        if isinstance(tol, type(None)):
            sub_pbed = self.slice(dir_id, val, verbose=verbose, level=level+1, log_delim=log_delim)
            r = sub_pbed.projection(dir_id, val, verbose=verbose, level=level+1, log_delim=log_delim)
        else:
            sub_pbed = self.slice(dir_id, val, verbose=verbose, tol=tol, level=level+1, log_delim=log_delim)
            r = np.array(sub_pbed.data.r)

        data = sub_pbed.data[sub_pbed.data.r != 0]

        x = np.array(data[xdir])
        y = np.array(data[ydir])

        if not isinstance(translation, type(None)):
            x += translation[0]
            y += translation[1]

        patches = []
        for i in range(len(data)):
            circle = Circle((x[i], y[i]), r[i])
            patches.append(circle)

        colors = np.array(data[field])

        if isinstance(clim, type(None)):
            clim = [data[field].min(), data[field].max()]
        if log_color:
            #sfmt = ticker.ScalarFormatter(useMathText=True)
            #sfmt.set_powerlimits((0, 0))
            p = PatchCollection(patches, cmap=colormap, norm=LogNorm())
        else:
            p = PatchCollection(patches, cmap=colormap)
        p.set_array(colors)
        p.set_clim(clim)

        if new_fig:
            if not isinstance(fig_size, type(None)):
                plt.figure(figsize=fig_size, facecolor = "white")
            else:
                plt.figure()

        if superimpose_Serpent:
            log_print(f'Superimposing Serpent geometry (x={Serpent_xsize}, y={Serpent_ysize}) from {Serpent_geom_path}', verbose, level+1, log_delim)
            self.show_Serpent_plot(Serpent_geom_path, new_fig=False, xlim=Serpent_xsize, ylim=Serpent_ysize, verbose=verbose, level=level+1, patch_clip=patch_clip, log_delim=log_delim)

        ax = plt.gca()
        ax.add_collection(p)
        if vertical_cbar:
            cbar = plt.colorbar(p, label=field_title, shrink=shrink_cbar, pad=pad_cbar)
        else:
            cbar = plt.colorbar(p, label=field_title, orientation="horizontal", shrink=shrink_cbar, pad=pad_cbar)
        cbar.formatter.set_powerlimits((-2, 4))
        if not isinstance(xlim, type(None)):
            ax.set_xlim(xlim)
        if not isinstance(ylim, type(None)):
            ax.set_ylim(ylim)
        xlab = xdir
        ylab = ydir
        if unit != '':
            xlab += f' [{unit}]'
            ylab += f' [{unit}]'
        plt.xlabel(xlab)
        plt.ylabel(ylab)
        if isinstance(plot_title, type(None)):
            plt.title("{}={:.3f}".format([dir_id], val))
        else:
            plt.title(plot_title)
        ax.autoscale_view()
        if equal:
            ax.set_aspect("equal", adjustable="box")
        plt.tight_layout()
        if no_ax:
            plt.axis('off')
        if save_fig:
            if isinstance(fig_name, type(None)):
                fig_path = f'{fig_folder}/2D_plot_{[dir_id]}{val:.2E}_{field}{fig_suffix}.png'
            else:
                fig_path = f'{fig_folder}/{fig_name}{fig_suffix}.png'
            log_print(f'Saving figure in {fig_path}', verbose, level+1, log_delim)
            plt.savefig(fig_path, dpi=fig_dpi, bbox_inches='tight')

        log_print(f'Done.', verbose, level, log_delim, end_block=True)
        return ax

    def plot3D(self, field='id', colormap='turbo', view=None, view_dist=None, xlim=None, ylim=None, zlim=None, transparent=False, translation=None, sample_fraction=None, fast=False, force_slow=False, scatter_size=10, alpha=1, show_ghosts=False, field_title=None, clim=None, shrink_cbar=1, pad_cbar=0.15, vertical_cbar=False, no_ax=False, plot_title=None, equal=True, fig_size=None, new_fig=True, save_fig=False, fig_folder='./', fig_name=None, fig_suffix='', fig_dpi=400, verbose=True, level=0, log_delim='  '):
        lim_fast = 1000

        data = self.data[self.data.r != 0]

        if not isinstance(sample_fraction, type(None)):
            data = data.sample(int(sample_fraction*len(data)))

        if isinstance(field, (tuple, list, np.ndarray, pd.DataFrame, pd.Series)):
            log_print(f'3D plotting pbed, showing user-defined array', verbose, level, log_delim)
            if isinstance(field_title, type(None)):
                field_title = 'user-defined'
                data[str(field_title)] = np.array(field)
                field = field_title
            else:
                data[field_title] = np.array(field)
                field = field_title
        else:
            if isinstance(field_title, type(None)):
                field_title = field
            log_print(f'3D plotting pbed, showing field {field}', verbose, level, log_delim)

        if not isinstance(sample_fraction, type(None)):
            log_print(f'Only showing {sample_fraction*100}% ({len(data)}) of the elements', verbose, level+1, log_delim)

        if isinstance(clim, type(None)):
            clim = [data[field].min(), data[field].max()]
        if new_fig:
            if not isinstance(fig_size, type(None)):
                plt.figure(figsize=fig_size, facecolor = "white")
            else:
                plt.figure()
            plt.gcf().add_subplot(111, projection="3d")
            #plt.gcf().add_subplot(111, projection="3d")  # , proj_type = 'ortho')

        ax = plt.gca()
        if len(data) > lim_fast and not force_slow:
            log_print(f'WARNING: Too many elements ({len(data)}>{lim_fast}), switched automatically to scatter plot, needing manual parameter for scatter size', verbose, level+1, log_delim)
            fast = True
            if isinstance(scatter_size, type(None)):
                raise Exception('Since the number of elements is too high, plot3D switched to scatter. Please put a size to elements, with scatter_size')

        if force_slow and len(data) > lim_fast:
            log_print(f'WARNING: Many elements to plot ({len(data)}>{lim_fast}), might take a while', verbose, level+1, log_delim)

        if not isinstance(translation, type(None)):
            data.x += translation[0]
            data.y += translation[1]
            data.z += translation[2]

        if not fast:
            values =  data[field]
            cmap = cm.get_cmap(colormap)
            norm = Normalize(vmin = clim[0], vmax = clim[1])
            normalized_values = norm(values)
            colors = cmap(normalized_values)

            i_row = 0
            for _, row in data.iterrows():
                # draw sphere
                u, v = np.mgrid[0:2*np.pi:10j, 0:np.pi:10j]
                x = row.r*np.cos(u)*np.sin(v)
                y = row.r*np.sin(u)*np.sin(v)
                z = row.r*np.cos(v)
                p = ax.plot_surface(x+row.x, y+row.y, z+row.z, color=colors[i_row], shade=False, alpha=alpha, zorder=1)
                i_row += 1
            if vertical_cbar:
                cb = plt.colorbar(p, label=field_title, shrink=shrink_cbar, pad=pad_cbar)
            else:
                cb = plt.colorbar(p, label=field_title, orientation="horizontal", shrink=shrink_cbar, pad=pad_cbar)
            cb.formatter.set_powerlimits((-2, 4))

        else:
            if isinstance(scatter_size, type(None)):
                raise Exception(f'Fast mode, needing manual parameter for scatter size')
            p = ax.scatter3D(
                data.x,
                data.y,
                data.z,
                s=scatter_size,
                c=data[field],
                alpha=alpha,
                zorder=1,
                vmin=clim[0],
                vmax=clim[1],
                cmap=colormap,
                edgecolor = None #'none'
            )
            if vertical_cbar:
                cb = plt.colorbar(p, label=field_title, shrink=shrink_cbar, pad=pad_cbar)
            else:
                cb = plt.colorbar(p, label=field_title, orientation="horizontal", shrink=shrink_cbar, pad=pad_cbar)
            cb.set_alpha(1)
            cb.formatter.set_powerlimits((-2, 4))
            cb.draw_all()

        ax.set_xlabel("x [cm]")
        ax.set_ylabel("y [cm]")
        ax.set_zlabel("z [cm]")

        if not isinstance(view, type(None)):
            ax.view_init(view[0], view[1])
        if not isinstance(view_dist, type(None)):
            ax.dist = view_dist
        if isinstance(plot_title, type(None)):
            plt.title(f"{field}")
        else:
            plt.title(plot_title)

        if equal:
            try:
                try:
                    ax.set_box_aspect((np.diff(ax.get_xlim())[0],np.diff(ax.get_ylim())[0],np.diff(ax.get_zlim())[0]))
                except:
                    ax.set_aspect('equal')
            except:
                pass
        ax.relim()      # make sure all the data fits
        ax.autoscale()  # auto-scale
        if not isinstance(xlim, type(None)):
            ax.set_xlim(xlim)
        if not isinstance(ylim, type(None)):
            ax.set_ylim(ylim)
        if not isinstance(zlim, type(None)):
            ax.set_zlim([zlim[0], zlim[1]])
        #ax.invert_zaxis()

        plt.tight_layout()
        if no_ax:
            plt.axis('off')
        if save_fig:
            if isinstance(fig_name, type(None)):
                fig_path = f'{fig_folder}/3D_plot_{field}{fig_suffix}.png'
            else:
                fig_path = f'{fig_folder}/{fig_name}{fig_suffix}.png'
            log_print(f'Saving figure in {fig_path}', verbose, level+1, log_delim)
            plt.savefig(fig_path, dpi=fig_dpi, bbox_inches='tight', transparent=transparent)

        log_print(f'Done.', verbose, level, log_delim, end_block=True)

        return ax

    def plot_summary(self, field='id', colormap='turbo', view=None, view_dist=10.5, xlim=None, ylim=None, zlim=None, translation=None, sample_fraction=None, fast=False, clipping_3D=True, shrink_cbar=1, pad_cbar=None, force_slow=False, scatter_size=10, alpha=1, equal=True, field_title=None, clim=None, superimpose_Serpent=False, Serpent_xsize=None, Serpent_ysize=None, Serpent_zsize=None, Serpent_paths_dirxyz = [None, None, None], patches_clip_2D=[None, None, None], save_fig=False, fig_size=None, fig_folder='./', fig_name=None, fig_suffix='', fig_dpi=400, verbose=True, level=0, log_delim='  '):
        if not isinstance(fig_size, type(None)):
            fig = plt.figure(figsize=fig_size)
        else:
            fig = plt.figure()

        log_print(f'Plotting of Pebble Bed', verbose, level, log_delim)
        if superimpose_Serpent:
            valx = np.mean(Serpent_xsize)
            valy = np.mean(Serpent_ysize)
            valz = np.mean(Serpent_zsize)
        else:
            valx = 'middle'
            valy = 'middle'
            valz = 'middle'
        if isinstance(clim, type(None)):
            clim = [np.min(self.data[field]), np.max(self.data[field])]

        ax1=fig.add_subplot(2,2,1)
        self.plot2D(field, dir_id=2, val=valz, plot_title='(X,Y)', translation=translation, superimpose_Serpent=superimpose_Serpent, Serpent_xsize=Serpent_xsize, Serpent_ysize=Serpent_ysize, Serpent_geom_path=Serpent_paths_dirxyz[2], new_fig=False, equal=equal, xlim=ylim, ylim=zlim, clim=clim, field_title=field_title, colormap=colormap, shrink_cbar=shrink_cbar, pad_cbar=pad_cbar, patch_clip=patches_clip_2D[2], verbose=verbose, level=level+1, log_delim=log_delim)
        ax2=fig.add_subplot(2,2,2)
        self.plot2D(field, dir_id=1, val=valy, plot_title='(X,Z)', translation=translation, superimpose_Serpent=superimpose_Serpent, Serpent_xsize=Serpent_xsize, Serpent_ysize=Serpent_zsize, Serpent_geom_path=Serpent_paths_dirxyz[1], new_fig=False, equal=equal, xlim=xlim, ylim=zlim, clim=clim, field_title=field_title, colormap=colormap, shrink_cbar=shrink_cbar, pad_cbar=pad_cbar, patch_clip=patches_clip_2D[1], verbose=verbose, level=level+1, log_delim=log_delim)
        ax3=fig.add_subplot(2,2,3)
        self.plot2D(field, dir_id=0, val=valx, plot_title='(Y,Z)', translation=translation, superimpose_Serpent=superimpose_Serpent, Serpent_xsize=Serpent_ysize, Serpent_ysize=Serpent_zsize, Serpent_geom_path=Serpent_paths_dirxyz[0], new_fig=False, equal=equal, xlim=xlim, ylim=ylim, clim=clim, field_title=field_title, colormap=colormap, shrink_cbar=shrink_cbar, pad_cbar=pad_cbar, patch_clip=patches_clip_2D[0], verbose=verbose, level=level+1, log_delim=log_delim)
        ax4=fig.add_subplot(2,2,4, projection="3d")
        if not clipping_3D:
            self.plot3D(field, new_fig=False, xlim=xlim, ylim=ylim, zlim=zlim, translation=translation, clim=clim, field_title=field_title, colormap=colormap, view=view, view_dist=view_dist, sample_fraction=sample_fraction, fast=fast, force_slow=False, equal=equal, scatter_size=scatter_size, alpha=alpha, shrink_cbar=shrink_cbar, pad_cbar=pad_cbar, verbose=verbose, level=level+1, log_delim=log_delim)
        else:
            clipped = self.clip(['x', 'y'], direction=[-1,+1], logical_and=False, verbose=verbose, level=level+1, log_delim=log_delim)
            clipped.plot3D(field, new_fig=False, xlim=xlim, ylim=ylim, zlim=zlim, translation=translation, clim=clim, field_title=field_title, colormap=colormap, view=view, sample_fraction=sample_fraction, fast=fast, force_slow=False, equal=equal, scatter_size=scatter_size, alpha=alpha, shrink_cbar=shrink_cbar, pad_cbar=pad_cbar, verbose=verbose, level=level+1, log_delim=log_delim)
        plt.tight_layout()
        if save_fig:
            if isinstance(fig_name, type(None)):
                fig_path = f'{fig_folder}/Summary_plot_{field}{fig_suffix}.png'
            else:
                fig_path = f'{fig_folder}/{fig_name}{fig_suffix}.png'
            log_print(f'Saving figure in {fig_path}', verbose, level+1, log_delim)
            plt.savefig(fig_path, dpi=fig_dpi, bbox_inches='tight')

        log_print(f'Done.', verbose, level, log_delim, end_block=True)
        return (ax1, ax2, ax3, ax4)

    def show_Serpent_plot(self, plot_file_path, title=None, new_fig=True, fig_size=None, save_fig=False, fig_folder='./', fig_name=None, fig_suffix='', fig_dpi=400, xlim=None, ylim=None, patch_clip=None, verbose=True, level=0, log_delim='  '):
        log_print(f'Showing Serpent plot from {plot_file_path}', verbose, level, log_delim)
        if isinstance(title, type(None)):
            title = plot_file_path

        if new_fig:
            if not isinstance(fig_size, type(None)):
                plt.figure(figsize=fig_size)
            else:
                plt.figure()

        ax = plt.gca()
        img = mpimg.imread(plot_file_path)
        if not isinstance(xlim, type(None)) and not isinstance(ylim, type(None)):
            imgplot = plt.imshow(img, extent=[xlim[0], xlim[1], ylim[0], ylim[1]], aspect='auto')
        else:
            imgplot = plt.imshow(img, aspect='auto')
        if not isinstance(patch_clip, type(None)):
            patch_clip.set_transform(plt.gca().transData)
            imgplot.set_clip_path(patch_clip)
        plt.title(title)
        plt.tight_layout()

        if save_fig:
            if isinstance(fig_name, type(None)):
                fig_path = f'{fig_folder}/Serpent_plot_{os.path.splitext(os.path.basename(plot_file_path))[0]}{fig_suffix}.png'
            else:
                fig_path = f'{fig_folder}/{fig_name}{fig_suffix}.png'
            log_print(f'Saving figure in {fig_path}', verbose, level+1, log_delim)
            plt.savefig(fig_path, dpi=fig_dpi, bbox_inches='tight')

        log_print(f'Done.', verbose, level, log_delim, end_block=True)
        return ax

    def construct_3D_Serpent_plot(self, middle_plane_xyz_plot_paths, xlim, ylim, zlim, npix=50, title=None, equal=True, new_fig=True, fig_size=None, save_fig=False, fig_folder='./', fig_name=None, fig_suffix='', fig_dpi=400, verbose=True, level=0, log_delim='  '):
        log_print(f'Contructing 3D geometry from Serpent plot from {middle_plane_xyz_plot_paths}', verbose, level, log_delim)

        if isinstance(title, type(None)):
            title = ''

        if new_fig:
            if not isinstance(fig_size, type(None)):
                plt.figure(figsize=fig_size)
                plt.gcf().add_subplot(111, projection="3d")
            else:
                plt.figure()
                plt.gcf().add_subplot(111, projection="3d")

        ax = plt.gca()
        for i, plot_file_path in enumerate(middle_plane_xyz_plot_paths):
            log_print(f'Picture {i}', verbose, level+1, log_delim)
            array = plt.imread(plot_file_path)
            array = np.swapaxes(array, 0, 1)
            array = np.flip(array, 1)
            rstride=int(array.shape[0]/npix)
            cstride=int(array.shape[1]/npix)
            if i == 0:
                x_1 = np.linspace(ylim[0], ylim[1], array.shape[0])
                y_1 = np.linspace(zlim[0], zlim[1], array.shape[1])
                y_1, x_1 = np.meshgrid(y_1, x_1)
                val = np.min(xlim)
                vals = np.ones((array.shape[0], array.shape[1]))*val
                to_plot = [vals, x_1, y_1]
            elif i == 1:
                x_1 = np.linspace(xlim[0], xlim[1], array.shape[0])
                y_1 = np.linspace(zlim[0], zlim[1], array.shape[1])
                y_1, x_1 = np.meshgrid(y_1, x_1)
                val = np.max(ylim)
                vals = np.ones((array.shape[0], array.shape[1]))*val
                to_plot = [x_1, vals, y_1]
            elif i == 2:
                x_1 = np.linspace(xlim[0], xlim[1], array.shape[0])
                y_1 = np.linspace(ylim[0], ylim[1], array.shape[1])
                y_1, x_1 = np.meshgrid(y_1, x_1)
                val = np.min(zlim)
                vals = np.ones((array.shape[0], array.shape[1]))*val
                to_plot = [x_1, y_1, vals]

            ax.plot_surface(to_plot[0], to_plot[1], to_plot[2], rstride=rstride, cstride=cstride, facecolors=array, zorder=i, shade=False)
        ax.set_xlabel('x [cm]')
        ax.set_ylabel('y [cm]')
        ax.set_zlabel('z [cm]')
        if equal:
            ax.set_box_aspect((np.diff(ax.get_xlim())[0],np.diff(ax.get_ylim())[0],np.diff(ax.get_zlim())[0]))
        plt.title(title)
        plt.tight_layout()

        if save_fig:
            if isinstance(fig_name, type(None)):
                fig_path = f'{fig_folder}/Serpent_plot_{os.path.splitext(os.path.basename(plot_file_path))[0]}{fig_suffix}.png'
            else:
                fig_path = f'{fig_folder}/{fig_name}{fig_suffix}.png'
            log_print(f'Saving figure in {fig_path}', verbose, level+1, log_delim)
            plt.savefig(fig_path, dpi=fig_dpi, bbox_inches='tight')

        log_print(f'Done.', verbose, level, log_delim, end_block=True)
        return ax

    #### WRITING ####

    def write_fields(self, path, fields='all', indices='all', separate_fields=False, separate_indices=False, sep =',', include_headers=True, include_indices=True, prefix='data', file_type='.csv', float_format='%.4E', verbose=True, level=0, log_delim='  '):
        subpbed = self.extract(fields=fields, indices=indices)
        if separate_fields:
            for field in subpbed.data.columns:
                if separate_indices:
                    for index in subpbed.data.index:
                        name = f'{path}/{prefix}_index_{index}_field_{field}{file_type}'
                        subpbed.data.loc[index, field].to_csv(name, header=include_headers, index=include_indices, sep=sep, float_format=float_format)
                else:
                    name = f'{path}/{prefix}_field_{field}{file_type}'
                    subpbed.data[field].to_csv(name, header=include_headers, index=include_indices, sep=sep, float_format=float_format)
        else:
            if separate_indices:
                for index in subpbed.data.index:
                    name = f'{path}/{prefix}_index_{index}{file_type}'
                    subpbed.data.loc[index].to_csv(name, header=include_headers, index=include_indices, sep=sep, float_format=float_format)
            else:
                name = f'{path}/{prefix}{file_type}'
                subpbed.data.to_csv(name, header=include_headers, index=include_indices, sep=sep, float_format=float_format)


    #### DOMAIN DECOMPOSITION ####

    def decompose_in_domains(self, n_domains_list, decomposition_types, filling_domain=None, center=None, shift_sectors=0, keep_subdomains=True, verbose=True, level=0, log_delim='  '):
        to_drop = [i for i in self.data.columns if 'domain_id' in i]
        if len(to_drop)>0:
            self.data.drop(columns=to_drop, inplace=True)

        if isinstance(n_domains_list, int):
            n_domains_list = [n_domains_list]
        N_levels = len(n_domains_list)
        if isinstance(shift_sectors, (int)):
            shift_sectors = [[shift_sectors for j in range(n_domains_list[i])] for i in range(N_levels)]
        elif isinstance(shift_sectors, (tuple, list, np.ndarray)):
            shift_sectors = [shift_sectors for i in range(N_levels)]

        if isinstance(center, (int, type(None))):
            center = [[center for j in range(n_domains_list[i])] for i in range(N_levels)]
        elif isinstance(center, (tuple, list, np.ndarray)):
            center = [center for i in range(N_levels)]

        if N_levels != len(decomposition_types):
            raise Exception('List of #domains must match decomposition types in size (please use 1 letter alias)')

        # First decomposition
        self.decompose_in_domains_simple(n_domains_list[0], decomposition_types[0], filling_domain, center[0][0], shift_sectors[0][0], verbose=verbose, level=level+1, log_delim=log_delim)
        self.data['domain_id_0'] = np.array(self.data.domain_id)

        # Rest of decompositions
        for i_decomp in range(1, N_levels):
            n_domains = n_domains_list[i_decomp]
            level_domains = n_domains_list[i_decomp-1]
            decomposition_type = decomposition_types[i_decomp]
            domain_id = np.ones(len(self.data))*-1
            for i_domain in range(level_domains):
                shift_sector = shift_sectors[i_decomp-1][i_domain]
                ctr = center[i_decomp-1][i_domain]
                sub_pbed = deepcopy(self)
                sub_pbed.data = sub_pbed.data[sub_pbed.data[f'domain_id_{i_decomp-1}'] == i_domain]
                # sub_pbed.data = sub_pbed.data.reset_index(drop=True)
                sub_pbed.decompose_in_domains_simple(n_domains, decomposition_type, filling_domain, ctr, shift_sector, idle=False, verbose=verbose, level=level+1, log_delim=log_delim)
                domain_id[sub_pbed.data.index] = sub_pbed.data.domain_id
            self.data[f'domain_id_{i_decomp}'] = np.array(domain_id)
        
        old_index = list(self.data.index)
        self.data = self.data.sort_values([f'domain_id_{i_decomp}' for i_decomp in range(N_levels)])
        self.data.domain_id = self.data.set_index([f'domain_id_{i_decomp}' for i_decomp in range(N_levels)]).index.factorize()[0]

        self.data = self.data.loc[old_index] #.sort_values('id')
        are_ghosts = self.is_ghost()
        self.data.loc[are_ghosts, 'domain_id'] = np.nan
        self.domains = dict()
        _, self.cnt_domains = np.unique(self.data.domain_id, return_counts=True)
        for i_domain in range(int(self.data.domain_id.max())+1):
            self.domains[i_domain] = self.cnt_domains[i_domain]
        log_print(f'Count in domains: {self.cnt_domains}', verbose, level+1, log_delim)

        if not keep_subdomains:
            self.data.drop(self.data.filter(like='domain_id_'), axis=1, inplace=True)

    def decompose_in_domains_simple(self, n_domains, decomposition_type, filling_domain=None, center=None, shift_sector=0, idle=True, verbose=True, level=0, log_delim='  '):

        if isinstance(filling_domain, type(None)):
            filling_domain = n_domains-1

        if not isinstance(center, type(None)):
            self.center = center
            self.process_data(calculate_center=False, verbose=verbose, level=level+1, log_delim=log_delim)
        else:
            self.process_data(verbose=verbose, level=level+1, log_delim=log_delim)

        decomposition_type = str(decomposition_type).lower()
        target_npebbles_domains = np.floor((len(self.data)- self.N_ghosts)/ n_domains).astype(int)

        if decomposition_type in ['0', 'random', 'n']:

            log_print(f'Decomposing pebbles randomly in {n_domains} domains', verbose, level, log_delim)
            domain_id = []
            i_pebble = 0
            are_ghosts = self.is_ghost()
            for i_domain in range(n_domains):
                cnt_domain = 0
                while cnt_domain < target_npebbles_domains:
                    if not are_ghosts[i_pebble]:
                        domain_id.append(i_domain)
                        cnt_domain += 1
                    else:
                        domain_id.append(np.nan)
                    i_pebble += 1

            # Handles last domain
            cnt_extra = 0
            while len(domain_id) < len(self.data):
                domain_id.append(filling_domain)
                cnt_extra += 1
            if cnt_extra > 0:
                log_print(f'{cnt_extra} extra pebbles in domain {filling_domain}', verbose, level+1, log_delim)

            # Shuffles
            np.random.shuffle(domain_id)

        elif decomposition_type in ['1', 'index', 'i', 'id']:

            log_print(f'Decomposing pebbles by index in {n_domains} domains', verbose, level, log_delim)
            domain_id = []
            i_pebble = 0
            are_ghosts = self.is_ghost()
            for i_domain in range(n_domains):
                cnt_domain = 0
                while cnt_domain < target_npebbles_domains:
                    if not are_ghosts[i_pebble]:
                        domain_id.append(i_domain)
                        cnt_domain += 1
                    else:
                        domain_id.append(np.nan)
                    i_pebble += 1

            # Handles last domain
            cnt_extra = 0
            while len(domain_id) < len(self.data):
                domain_id.append(filling_domain)
                cnt_extra += 1
            if cnt_extra > 0:
                log_print(f'{cnt_extra} extra pebbles in domain {filling_domain}', verbose, level+1, log_delim)

        elif decomposition_type in ['2', 'sectors', 'sector', 's']:

            log_print(f'Decomposing pebbles by azimuthal sectors in {n_domains} domains', verbose, level, log_delim)
            if shift_sector > 0:
                log_print(f'Shifting sectors by {shift_sector} degree', verbose, level+1, log_delim)

            if 'azim_angle' not in self.data.columns:
                self.process_data(verbose, level+1, log_delim)

            # Shift
            angles = self.data.azim_angle - shift_sector
            angles[angles < 0] += 360

            angles = angles.sort_values()
            domain_id = pd.DataFrame(-np.ones(angles.shape), index=angles.index)
            i_pebble = 0
            are_ghosts = self.is_ghost()
            for i_domain in range(n_domains):
                cnt_domain = 0
                while cnt_domain < target_npebbles_domains:
                    if not are_ghosts[angles.index[i_pebble]]:
                        domain_id.loc[angles.index[i_pebble]] = i_domain
                        cnt_domain += 1
                    else:
                        domain_id.loc[angles.index[i_pebble]] = np.nan
                    i_pebble += 1

            # Handles last domain
            domain_id = np.array(domain_id.sort_index())
            negative_id = domain_id < 0
            cnt_extra = np.sum(negative_id)
            domain_id[negative_id] = filling_domain
            if cnt_extra > 0:
                log_print(f'{cnt_extra} extra pebbles in domain {filling_domain}', verbose, level+1, log_delim)

        elif decomposition_type in ['4', 'radial', 'rad', 'r']:

            log_print(f'Decomposing pebbles by radial zones in {n_domains} domains', verbose, level, log_delim)

            if 'r_dist' not in self.data.columns:
                self.process_data(verbose, level+1, log_delim)

            radial_dist = self.data.r_dist.sort_values()
            domain_id = pd.DataFrame(-np.ones(radial_dist.shape), index=radial_dist.index)
            i_pebble = 0
            are_ghosts = self.is_ghost()
            for i_domain in range(n_domains):
                cnt_domain = 0
                while cnt_domain < target_npebbles_domains:
                    if not are_ghosts[radial_dist.index[i_pebble]]:
                        domain_id.loc[radial_dist.index[i_pebble]] = i_domain
                        cnt_domain += 1
                    else:
                        domain_id.loc[radial_dist.index[i_pebble]] = np.nan
                    i_pebble += 1

            # Handles last domain
            domain_id = np.array(domain_id.sort_index())
            negative_id = domain_id < 0
            cnt_extra = np.sum(negative_id)
            domain_id[negative_id] = filling_domain
            if cnt_extra > 0:
                log_print(f'{cnt_extra} extra pebbles in domain {filling_domain}', verbose, level+1, log_delim)

        elif decomposition_type in ['5', 'axial', 'ax', 'a']:

            log_print(f'Decomposing pebbles by axial zones in {n_domains} domains', verbose, level, log_delim)

            axial_dist = self.data.z.sort_values()
            domain_id = pd.DataFrame(-np.ones(axial_dist.shape), index=axial_dist.index)
            i_pebble = 0
            are_ghosts = self.is_ghost()
            for i_domain in range(n_domains):
                cnt_domain = 0
                while cnt_domain < target_npebbles_domains:
                    if not are_ghosts[axial_dist.index[i_pebble]]:
                        domain_id.loc[axial_dist.index[i_pebble]] = i_domain
                        cnt_domain += 1
                    else:
                        domain_id.loc[axial_dist.index[i_pebble]] = np.nan
                    i_pebble += 1

            # Handles last domain
            domain_id = np.array(domain_id.sort_index())
            negative_id = domain_id < 0
            cnt_extra = np.sum(negative_id)
            domain_id[negative_id] = filling_domain
            if cnt_extra > 0:
                log_print(f'{cnt_extra} extra pebbles in domain {filling_domain}', verbose, level+1, log_delim)

        elif decomposition_type in ['6', 'spherical', 'spheric', 'sphere', 'sph', 'o']:

            log_print(f'Decomposing pebbles by spherical zones in {n_domains} domains', verbose, level, log_delim)

            if 'dist' not in self.data.columns:
                self.process_data(verbose, level+1, log_delim)


            dist = self.data.dist.sort_values()
            domain_id = pd.DataFrame(-np.ones(dist.shape), index=dist.index)
            i_pebble = 0
            are_ghosts = self.is_ghost()
            for i_domain in range(n_domains):
                cnt_domain = 0
                while cnt_domain < target_npebbles_domains:
                    if not are_ghosts[dist.index[i_pebble]]:
                        domain_id.loc[dist.index[i_pebble]] = i_domain
                        cnt_domain += 1
                    else:
                        domain_id.loc[dist.index[i_pebble]] = np.nan
                    i_pebble += 1

            # Handles last domain
            domain_id = np.array(domain_id.sort_index())
            negative_id = domain_id < 0
            cnt_extra = np.sum(negative_id)
            domain_id[negative_id] = filling_domain
            if cnt_extra > 0:
                log_print(f'{cnt_extra} extra pebbles in domain {filling_domain}', verbose, level+1, log_delim)

        self.data['domain_id'] = domain_id
        

        if idle:
            self.domains = dict()
            _, cnt_domains = np.unique(self.data.domain_id, return_counts=True)
            for i_domain in range(n_domains):
                self.domains[i_domain] = cnt_domains[i_domain]
            log_print(f'Count in domains: {cnt_domains}', verbose, level+1, log_delim)
        log_print(f'Done.', verbose, level, log_delim, end_block=True)


def log_print(text, printing=True, level=0, log_delim='  ', line_limit=None, end_block=False, end_block_log_delim='\n', returning=True):
    # Line breaks
    if not isinstance(line_limit, type(None)):
        current_word = ''
        current_line = level*log_delim
        s = ''

        # Loop over message
        for c in text:
            if c != ' ':  # If not end of word
                if c == '\'':
                    c = '\"'
                current_word += c # Add to word
            else: # If end of word
                if len(current_line) + len(current_word) <= line_limit: # and if within line limits even when adding word
                    if current_line == level*log_delim:
                        current_line += current_word
                    else:
                        current_line += ' ' + current_word
                else: # if when adding word out of limits
                    s += current_line + '\n' # write line go to next line
                    current_line = level*log_delim + current_word # add word to next line
                current_word = '' # reset word
        s += current_line + ' ' + current_word # add last line
    else:
        s = level*log_delim + text

    if end_block and level==0:
        s += end_block_log_delim
    if printing:
        print(s)
    if returning:
        return s

def find_matching(string, list_strings, infer_if_incorrect=True, lower_ok=True, cutoff=0.4):
    list_strings = list(list_strings)
    if lower_ok:
        keys = [key.lower() for key in list_strings]
        search_parameter = string.lower()
    else:
        keys = [key for key in list_strings]
        search_parameter = str(string)

    if search_parameter not in keys:
        if infer_if_incorrect:
            match_string = difflib.get_close_matches(search_parameter, keys, cutoff=cutoff, n=100000)
            if len(match_string) == 0:
                return None
            match_indices = [keys.index(key) for key in match_string]
            real_matches = [list_strings[i] for i in match_indices]
            print(f'String {string} not found, close matches are: {real_matches}. Picking the first one.')
            string = real_matches[0]
        else:
            return None
    else:
        index = keys.index(string.lower())
        string = list_strings[index]
    return string

def errorbar(x, y, yerr, steps=False, fancy=True, label=None, alpha=0.3, color=None, labelize_error=False, force_continuous=True, linewidth=1, linestyle='-'):
    if isinstance(label, type(None)):
        label='Values'
    if isinstance(labelize_error, str):
        label_error = labelize_error
    elif labelize_error:
        label_error = 'Error'
    else:
        label_error = None

    x = np.array(x)
    y = np.array(y)
    yerr = np.array(yerr)
    if force_continuous:
        nans = np.isnan(y)
        y[nans] = np.interp(x[nans], x[~nans], y[~nans])
        if np.array(yerr).ndim <= 1:
            yerr[nans] = np.interp(x[nans], x[~nans], yerr[~nans])
        else:
            for i in range(yerr.shape[1]):
                yerr[nans, i] = np.interp(x[nans], x[~nans], yerr[~nans, i])

    if not fancy:
        if steps:
            plt.errorbar(x, y, yerr, drawstyle='steps', label=label, color=color, linewidth=linewidth, linestyle=linestyle)
        else:
            plt.errorbar(x, y, yerr, label=label, color=color, linewidth=linewidth, linestyle=linestyle)
    else:
        if np.array(yerr).ndim <= 1:
            yerr_low = yerr
            yerr_high = yerr
        elif np.array(yerr).ndim == 2:
            yerr_low = yerr[:,0]
            yerr_high = yerr[:,1]
        if steps:
            plt.plot(x, y, drawstyle='steps', label=label, zorder=1, color=color, linewidth=linewidth, linestyle=linestyle)
            plt.fill_between(x, y-yerr_low, y+yerr_high, step="pre", label=label_error, alpha=alpha, color=plt.gca().lines[-1].get_color())
        else:
            plt.plot(x, y, label=label, zorder=1, color=color, linewidth=linewidth, linestyle=linestyle)
            plt.fill_between(x, y-yerr_low, y+yerr_high, label=label_error, alpha=alpha, color=plt.gca().lines[-1].get_color())
        plt.legend()



def plot_enveloppe(yvalues, xvalues=None, new_fig=True, show_std=True, show_extremes=True):
    if hasattr(yvalues, 'ngroups'):
        mean = np.ravel(yvalues.mean())
        std = np.ravel(yvalues.std())
        mini = np.ravel(yvalues.min())
        maxi = np.ravel(yvalues.max())
        if isinstance(xvalues, type(None)):
            xvalues = yvalues.groups.keys()
    else:
        mean = np.mean(yvalues)
        std = np.std(yvalues)
        mini = np.min(yvalues)
        maxi = np.max(yvalues)
        if isinstance(xvalues, type(None)):
            xvalues = np.arange(len(mean))
    if new_fig:
        plt.figure()
    if show_std:
        errorbar(xvalues, mean,  std)
    else:
        plt.plot(xvalues, mean, label='Values')
    if show_extremes:
        plt.plot(xvalues, mini, ':' , label='Min', color=plt.gca().lines[-1].get_color())
        plt.plot(xvalues, maxi, ':', label='Max', color=plt.gca().lines[-1].get_color())


    plt.legend()

def swap_axis(ax):
    xlabel, ylabel = ax.get_xlabel(), ax.get_ylabel()
    plt.xlabel(ylabel)
    plt.ylabel(xlabel)

    xlim, ylim = ax.get_xlim(),  ax.get_ylim()
    plt.xlim(ylim)
    plt.ylim(xlim)
    for lines in ax.get_lines():
        try:
            iter(lines)
        except:
            lines = [lines]
        for line in lines:
            xdata, ydata = line.get_xdata(), line.get_ydata()
            line.set_xdata(ydata)
            line.set_ydata(xdata)
            line.axes.autoscale_view()

    return ax

turbo_colormap_data = np.array([[0.18995,0.07176,0.23217], [0.19483,0.08339,0.26149], [0.19956,0.09498,0.29024], [0.20415,0.10652,0.31844], [0.20860,0.11802,0.34607], [0.21291,0.12947,0.37314], [0.21708,0.14087,0.39964], [0.22111,0.15223,0.42558], [0.22500,0.16354,0.45096], [0.22875,0.17481,0.47578], [0.23236,0.18603,0.50004], [0.23582,0.19720,0.52373], [0.23915,0.20833,0.54686], [0.24234,0.21941,0.56942], [0.24539,0.23044,0.59142], [0.24830,0.24143,0.61286], [0.25107,0.25237,0.63374], [0.25369,0.26327,0.65406], [0.25618,0.27412,0.67381], [0.25853,0.28492,0.69300], [0.26074,0.29568,0.71162], [0.26280,0.30639,0.72968], [0.26473,0.31706,0.74718], [0.26652,0.32768,0.76412], [0.26816,0.33825,0.78050], [0.26967,0.34878,0.79631], [0.27103,0.35926,0.81156], [0.27226,0.36970,0.82624], [0.27334,0.38008,0.84037], [0.27429,0.39043,0.85393], [0.27509,0.40072,0.86692], [0.27576,0.41097,0.87936], [0.27628,0.42118,0.89123], [0.27667,0.43134,0.90254], [0.27691,0.44145,0.91328], [0.27701,0.45152,0.92347], [0.27698,0.46153,0.93309], [0.27680,0.47151,0.94214], [0.27648,0.48144,0.95064], [0.27603,0.49132,0.95857], [0.27543,0.50115,0.96594], [0.27469,0.51094,0.97275], [0.27381,0.52069,0.97899], [0.27273,0.53040,0.98461], [0.27106,0.54015,0.98930], [0.26878,0.54995,0.99303], [0.26592,0.55979,0.99583], [0.26252,0.56967,0.99773], [0.25862,0.57958,0.99876], [0.25425,0.58950,0.99896], [0.24946,0.59943,0.99835], [0.24427,0.60937,0.99697], [0.23874,0.61931,0.99485], [0.23288,0.62923,0.99202], [0.22676,0.63913,0.98851], [0.22039,0.64901,0.98436], [0.21382,0.65886,0.97959], [0.20708,0.66866,0.97423], [0.20021,0.67842,0.96833], [0.19326,0.68812,0.96190], [0.18625,0.69775,0.95498], [0.17923,0.70732,0.94761], [0.17223,0.71680,0.93981], [0.16529,0.72620,0.93161], [0.15844,0.73551,0.92305], [0.15173,0.74472,0.91416], [0.14519,0.75381,0.90496], [0.13886,0.76279,0.89550], [0.13278,0.77165,0.88580], [0.12698,0.78037,0.87590], [0.12151,0.78896,0.86581], [0.11639,0.79740,0.85559], [0.11167,0.80569,0.84525], [0.10738,0.81381,0.83484], [0.10357,0.82177,0.82437], [0.10026,0.82955,0.81389], [0.09750,0.83714,0.80342], [0.09532,0.84455,0.79299], [0.09377,0.85175,0.78264], [0.09287,0.85875,0.77240], [0.09267,0.86554,0.76230], [0.09320,0.87211,0.75237], [0.09451,0.87844,0.74265], [0.09662,0.88454,0.73316], [0.09958,0.89040,0.72393], [0.10342,0.89600,0.71500], [0.10815,0.90142,0.70599], [0.11374,0.90673,0.69651], [0.12014,0.91193,0.68660], [0.12733,0.91701,0.67627], [0.13526,0.92197,0.66556], [0.14391,0.92680,0.65448], [0.15323,0.93151,0.64308], [0.16319,0.93609,0.63137], [0.17377,0.94053,0.61938], [0.18491,0.94484,0.60713], [0.19659,0.94901,0.59466], [0.20877,0.95304,0.58199], [0.22142,0.95692,0.56914], [0.23449,0.96065,0.55614], [0.24797,0.96423,0.54303], [0.26180,0.96765,0.52981], [0.27597,0.97092,0.51653], [0.29042,0.97403,0.50321], [0.30513,0.97697,0.48987], [0.32006,0.97974,0.47654], [0.33517,0.98234,0.46325], [0.35043,0.98477,0.45002], [0.36581,0.98702,0.43688], [0.38127,0.98909,0.42386], [0.39678,0.99098,0.41098], [0.41229,0.99268,0.39826], [0.42778,0.99419,0.38575], [0.44321,0.99551,0.37345], [0.45854,0.99663,0.36140], [0.47375,0.99755,0.34963], [0.48879,0.99828,0.33816], [0.50362,0.99879,0.32701], [0.51822,0.99910,0.31622], [0.53255,0.99919,0.30581], [0.54658,0.99907,0.29581], [0.56026,0.99873,0.28623], [0.57357,0.99817,0.27712], [0.58646,0.99739,0.26849], [0.59891,0.99638,0.26038], [0.61088,0.99514,0.25280], [0.62233,0.99366,0.24579], [0.63323,0.99195,0.23937], [0.64362,0.98999,0.23356], [0.65394,0.98775,0.22835], [0.66428,0.98524,0.22370], [0.67462,0.98246,0.21960], [0.68494,0.97941,0.21602], [0.69525,0.97610,0.21294], [0.70553,0.97255,0.21032], [0.71577,0.96875,0.20815], [0.72596,0.96470,0.20640], [0.73610,0.96043,0.20504], [0.74617,0.95593,0.20406], [0.75617,0.95121,0.20343], [0.76608,0.94627,0.20311], [0.77591,0.94113,0.20310], [0.78563,0.93579,0.20336], [0.79524,0.93025,0.20386], [0.80473,0.92452,0.20459], [0.81410,0.91861,0.20552], [0.82333,0.91253,0.20663], [0.83241,0.90627,0.20788], [0.84133,0.89986,0.20926], [0.85010,0.89328,0.21074], [0.85868,0.88655,0.21230], [0.86709,0.87968,0.21391], [0.87530,0.87267,0.21555], [0.88331,0.86553,0.21719], [0.89112,0.85826,0.21880], [0.89870,0.85087,0.22038], [0.90605,0.84337,0.22188], [0.91317,0.83576,0.22328], [0.92004,0.82806,0.22456], [0.92666,0.82025,0.22570], [0.93301,0.81236,0.22667], [0.93909,0.80439,0.22744], [0.94489,0.79634,0.22800], [0.95039,0.78823,0.22831], [0.95560,0.78005,0.22836], [0.96049,0.77181,0.22811], [0.96507,0.76352,0.22754], [0.96931,0.75519,0.22663], [0.97323,0.74682,0.22536], [0.97679,0.73842,0.22369], [0.98000,0.73000,0.22161], [0.98289,0.72140,0.21918], [0.98549,0.71250,0.21650], [0.98781,0.70330,0.21358], [0.98986,0.69382,0.21043], [0.99163,0.68408,0.20706], [0.99314,0.67408,0.20348], [0.99438,0.66386,0.19971], [0.99535,0.65341,0.19577], [0.99607,0.64277,0.19165], [0.99654,0.63193,0.18738], [0.99675,0.62093,0.18297], [0.99672,0.60977,0.17842], [0.99644,0.59846,0.17376], [0.99593,0.58703,0.16899], [0.99517,0.57549,0.16412], [0.99419,0.56386,0.15918], [0.99297,0.55214,0.15417], [0.99153,0.54036,0.14910], [0.98987,0.52854,0.14398], [0.98799,0.51667,0.13883], [0.98590,0.50479,0.13367], [0.98360,0.49291,0.12849], [0.98108,0.48104,0.12332], [0.97837,0.46920,0.11817], [0.97545,0.45740,0.11305], [0.97234,0.44565,0.10797], [0.96904,0.43399,0.10294], [0.96555,0.42241,0.09798], [0.96187,0.41093,0.09310], [0.95801,0.39958,0.08831], [0.95398,0.38836,0.08362], [0.94977,0.37729,0.07905], [0.94538,0.36638,0.07461], [0.94084,0.35566,0.07031], [0.93612,0.34513,0.06616], [0.93125,0.33482,0.06218], [0.92623,0.32473,0.05837], [0.92105,0.31489,0.05475], [0.91572,0.30530,0.05134], [0.91024,0.29599,0.04814], [0.90463,0.28696,0.04516], [0.89888,0.27824,0.04243], [0.89298,0.26981,0.03993], [0.88691,0.26152,0.03753], [0.88066,0.25334,0.03521], [0.87422,0.24526,0.03297], [0.86760,0.23730,0.03082], [0.86079,0.22945,0.02875], [0.85380,0.22170,0.02677], [0.84662,0.21407,0.02487], [0.83926,0.20654,0.02305], [0.83172,0.19912,0.02131], [0.82399,0.19182,0.01966], [0.81608,0.18462,0.01809], [0.80799,0.17753,0.01660], [0.79971,0.17055,0.01520], [0.79125,0.16368,0.01387], [0.78260,0.15693,0.01264], [0.77377,0.15028,0.01148], [0.76476,0.14374,0.01041], [0.75556,0.13731,0.00942], [0.74617,0.13098,0.00851], [0.73661,0.12477,0.00769], [0.72686,0.11867,0.00695], [0.71692,0.11268,0.00629], [0.70680,0.10680,0.00571], [0.69650,0.10102,0.00522], [0.68602,0.09536,0.00481], [0.67535,0.08980,0.00449], [0.66449,0.08436,0.00424], [0.65345,0.07902,0.00408], [0.64223,0.07380,0.00401], [0.63082,0.06868,0.00401], [0.61923,0.06367,0.00410], [0.60746,0.05878,0.00427], [0.59550,0.05399,0.00453], [0.58336,0.04931,0.00486], [0.57103,0.04474,0.00529], [0.55852,0.04028,0.00579], [0.54583,0.03593,0.00638], [0.53295,0.03169,0.00705], [0.51989,0.02756,0.00780], [0.50664,0.02354,0.00863], [0.49321,0.01963,0.00955], [0.47960,0.01583,0.01055]])

def RGBToPyCmap(rgbdata):
    nsteps = rgbdata.shape[0]
    stepaxis = np.linspace(0, 1, nsteps)

    rdata=[]; gdata=[]; bdata=[]
    for istep in range(nsteps):
        r = rgbdata[istep,0]
        g = rgbdata[istep,1]
        b = rgbdata[istep,2]
        rdata.append((stepaxis[istep], r, r))
        gdata.append((stepaxis[istep], g, g))
        bdata.append((stepaxis[istep], b, b))

    mpl_data = {'red':   rdata,
                 'green': gdata,
                 'blue':  bdata}

    return mpl_data

try:
    mpl_data = RGBToPyCmap(turbo_colormap_data)
    plt.register_cmap(name='turbo', data=mpl_data, lut=turbo_colormap_data.shape[0])
    mpl_data_r = RGBToPyCmap(turbo_colormap_data[::-1,:])
    plt.register_cmap(name='turbo_r', data=mpl_data_r, lut=turbo_colormap_data.shape[0])
except:
    pass

if __name__ == '__main__':
    pbed_file = '/home/yryves/serpent_cases/perf_pbed/fpb_pos_500000'
    for step in range(2):
        det_file = f'/home/yryves/serpent_cases/perf_pbed/ndets/ndets3/core500000/input_det{step}.m'
        pbed = Pebble_bed()
        pbed.read_pbed_file(pbed_file)
        pbed.read_detector(det_file)
        pbed.add_field('power_err_pc', pbed.data['power_err']*100)
        clipped = pbed.clip([0, 1], direction=[-1,+1], logical_and=False)
        clipped.plot3D('flux_thermal_E0', fig_size=(6,6), scatter_size=2, shrink_cbar=0.6, pad_cbar=-0.1, save_fig=True, fig_dpi=200, fig_suffix=f'_{step}', no_ax=True, plot_title=f'Step {step}')
        clipped.plot3D('flux_fast_E0', fig_size=(6,6), scatter_size=2, shrink_cbar=0.6, pad_cbar=-0.1, save_fig=True, fig_dpi=200, fig_suffix=f'_{step}', no_ax=True, plot_title=f'Step {step}')
        clipped.plot3D('power', fig_size=(6,6), scatter_size=2, shrink_cbar=0.6, pad_cbar=-0.1, save_fig=True, fig_dpi=200, fig_suffix=f'_{step}', no_ax=True, plot_title=f'Step {step}')
        clipped.plot3D('power_err_pc', fig_size=(6,6), scatter_size=2, shrink_cbar=0.6, pad_cbar=-0.1, save_fig=True, fig_dpi=200, fig_suffix=f'_{step}', no_ax=True, plot_title=f'Step {step}')


#     file = '/home/yryves/serpent_cases/domain_decomposition_dvlpt/test_stl3/Data/data_141.csv'
#     verbose=True
#     pbed = Pebble_bed()
#     pbed.read_datafile(file, verbose=verbose)
#     #pbed.plot_distribution('burnup', normalized=True, save_fig=True)
#     #pbed.plot_res('anakeff', save_fig=True, which_columns=0, which_columns_errors=1, res_file_path='/home/yryves/serpent_cases/domain_decomposition_dvlpt/test_stl3/wrk_Serpent/input_res.m')
#     #pbed.plot1D_mean('burnup', fieldx='z', save_fig=True)
#     #pbed.plot1D_mean('power', error_field='power_err', save_fig=True)
#     #pbed.plot1D_sum('power', error_field='power_err', save_fig=True)

#     Rgeom_lim = [-144.0, 144.0] # cm
#     Zgeom_lim = [0, 487.528]
#     folder = '/home/yryves/serpent_cases/domain_decomposition_dvlpt/test_stl3/'
#     plot_files = [f'{folder}/Plots/input_geom1_bu141.png', f'{folder}/Plots/input_geom2_bu141.png', f'{folder}/Plots/input_geom3_bu141.png']
#     pbed.read_datafile(f'{folder}/Data/data_142.csv')
#     pbed.add_field('powerkW', pbed.data.power/1000)
#     #pbed.plot_summary('powerkW', save_fig=True, fig_size=(10,10), shrink_cbar=0.7, alpha=1, superimpose_Serpent=True, Serpent_xsize=Rgeom_lim, Serpent_ysize=Rgeom_lim, Serpent_zsize=Zgeom_lim, Serpent_paths_dirxyz=plot_files, scatter_size=10, field_title='Pebble power [kW]')

#     #pbed.plot1D_mean('powerkW')
#     # pbed.data['power_contribution'] = pbed.data.power/pbed.data.power.sum()*100
#     # pbed.plot1D_sum('power_contribution', fieldx='z', bins=100, fig_size=(4,6), cumulative=False, save_fig=True, plot_title='', ylabel='Power contribution [%]', xlabel='Elevation [cm]', label='Bin power')
#     # ax1 = plt.gca()
#     # ax1 = swap_axis(ax1)

#     # ax2 = ax1.twiny()
#     # pbed.plot1D_sum('power_contribution', new_fig=False, fieldx='z', color='r', bins=10000, fig_size=(6,3), cumulative=True, save_fig=True, plot_title='', ylabel='Cumulative power contribution [%]', xlabel='Elevation [cm]', label='Cumulative power')
#     # ax2 = swap_axis(ax2)
#     # ax1.set_xlim([0, 2.5])
#     # ax2.set_xlim([0, 100])
#     # ax2.set_axisbelow(True)

#     # plt.axhline(100, linestyle=':', color='k', label='Defueling chute')
#     # h1, l1 = ax1.get_legend_handles_labels()
#     # h2, l2 = ax2.get_legend_handles_labels()
#     # ax2.legend(h1+h2, l1+l2, loc='lower left')
#     # ax1.set_axisbelow(True)

#     # plt.savefig('power_contribution.png', dpi=600)
#     # print(pbed.data.loc[pbed.data.z<100].sort_values('z').power_contribution.cumsum())

#     folder = "/home/yryves/serpent_cases/domain_decomposition_dvlpt/test_stl3/"
#     pbed_path = folder + "fpb_pos"
#     # det0_path = folder + "input_det0.m"
#     dep_path  = folder + "wrk_Serpent/input_dep.m"
#     res_path  = folder + "wrk_Serpent/input_res.m"
#     data_path = folder + "Data/data_142.csv"
#     # plotxy_path  = folder + "input_geom3.png"
#     # plotyz_path  = folder + "input_geom4.png"
#     # R = 80
#     # H = 190
#     dd=False
#     fuel_material='fuel'
#     pbed = Pebble_bed()
#     pbed.read_datafile(data_path, verbose=verbose)
#     # dz = 10
#     # ori_r = np.array(pbed.data.r)
#     # for i in range(10):
#         # pbed.data.z += dz
#         # pbed.data.loc[pbed.data.z + ori_r > 150, "z"] -= (150-40)
#         # pbed.data.loc[pbed.data.z - ori_r < 40, "r"] = 0
#         # pbed.data.loc[pbed.data.z - ori_r >= 40, "r"] = ori_r[pbed.data.z - ori_r >= 40]
#     # plt.figure()
#     # ax= pbed.show_Serpent_plot(plotxy_path, xlim=[-R,R], ylim=[-R,R], verbose=verbose) #, tol=1000)
#     # plt.show()
# #
#     # plt.figure()
#     # ax= pbed.show_Serpent_plot(plotyz_path, xlim=[-R,R], ylim=[0,H], verbose=verbose) #, tol=1000)
#     # plt.show()
#     # ax= pbed.plot2D(field='z', verbose=verbose) #, tol=1000)
#     # plt.show()
#     # pbed.read_detector(det0_path, which_dets='all', verbose=verbose)
#     # ax= pbed.plot2D(field='id', verbose=verbose) #, tol=1000)
#     # plt.show()
#     pbed.read_depletion(dep_path, material_name=fuel_material, dd=dd, fields='all', verbose=verbose)
#     # ax = pbed.plot3D(force_slow=True, sample_fraction=1, field='id', verbose=verbose)
#     # plt.show()
#     # error
#     # pbed.clip('rdist', 150, +1).clip('rdist', 155, -1)
#     # plt.figure()
#     # ax = pbed.plot2D(field='dist', field_title='Distance to center')
#     # plt.show()
#     # plt.figure()
#     # ax= pbed.plot2D(field='r_dist', field_title='Radial distance to center')
#     # plt.show()

# # verbose=True
# # folder = "/home/yryves/serpent_cases/test_larger/"
# # pbed_path = folder + "fpb_pos"
# # pbed = Pebble_bed()
# # pbed.read_pbed_file(pbed_path)
# # pbed.decompose_in_domains_simple(5, 'axial')
# # pbed.add_field('angle', np.arctan2(pbed.data.z-pbed.center[2], pbed.data.r_dist))
# # pbed.plot2D('angle', dir_id='x', val=0, fig_size=(3,5))
# # pbed.plot2D('domain_id', dir_id='x', val=0, fig_size=(3,5))



# # pbed.decompose_in_domains(4, 'r', verbose=False)
# # pbed.plot3D('domain_id', fig_size=(6,8), scatter_size=50, save_fig=True, fig_suffix=f'_o', field_title=f'DD: {t}', verbose=False)
# # plt.show()

# # pbed.decompose_in_domains(6, 'o', verbose=False)
# # pbed.plot3D('domain_id', fig_size=(6,8), scatter_size=50, save_fig=True, fig_suffix=f'_o', field_title=f'DD: {t}', verbose=False)
# # plt.show()

# # pbed.decompose_in_domains([2, 6], 'rs', verbose=False)
# # pbed.plot3D('domain_id', fig_size=(6,8), scatter_size=50, save_fig=True, fig_suffix=f'_rs', field_title=f'DD: {t}', verbose=False)
# # plt.show()

# # pbed.decompose_in_domains([4, 4], 'as', shift_sectors=[90/4*i for i in range(4)], verbose=False)
# # pbed.plot3D('domain_id', fig_size=(6,8), scatter_size=50, save_fig=True, fig_suffix=f'_as', field_title=f'DD: {t}', verbose=False)
# # plt.show()

# # pbed.decompose_in_domains(8, 's', center=[[40,0,0]], verbose=False)
# # pbed.plot3D('domain_id', fig_size=(6,8), scatter_size=50, save_fig=True, fig_suffix=f'_s_shift', field_title=f'DD: {t}', verbose=False)
# # plt.show()

# #pbed.clip(0, direction=-1).plot3D('angle', alpha=0.6, save_fig=True, fig_size=(3,5), xlim=[-40,40])