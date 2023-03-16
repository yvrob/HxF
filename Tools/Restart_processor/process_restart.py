#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__author__ = "Yves Robert"
__credits__ = ["Tatiana Siarafera"]
__version__ = "1.0"
__maintainer__ = "Yves Robert"
__email__ = "yves.robert@berkeley.edu"

""" Description: Library to write and read Serpent restart files, and process them"""

#%% Modules
import struct
import os
import matplotlib.pyplot as plt # can be commented if not installed
import pandas as pd # can be commented if not installed
import numpy as np # can be commented if not installed
import re

#%% Classes definition
class Restart_File:
    """Restart file objects are objects which read/write materials compositions at different points in time. 
       
       The structure is: 
       - Restart file (example: restart): contains snapshots, burnup points, time points. All dictionnaries
            - Snapshot 0: restart[0]. Contains all material objects in snapshot 0, which have the burnup 0 and time 0.
                - Material 0 (example: "mat0"): Material object containing material and nuclides information at given burnup and time.
                - Material 1 (example: "mat1")
                ...
            - Snapshot 1: restart[1]
            ...
        Initialization: only one initialization is needed and possible.
        1) If initialized with a path, it creates an empty object linked to the path, empty of snapshots.
            restart = Restart_File(path_to_file=path_in)
          Then to read the compositions path: 
            restart.read_restart()
          It will fill the object with all the snapshots contained in the restart file.
        2) If initialized with a list of snapshots (which are dictionnaries of Material objects), it creates an object with the corresponding snapshots.
            restart = Restart_File(snapshots=list_snapshots)
           It then needs to be linked to a path to write a binary file or text file.
            restart.path_to_file = path
           Binary: 
            restart.write_binary()
           Text:
            restart.write_text()
        3) If initialized with a dictionary Material objects (aka a snapshot), it creates a snapshot 0 with the materials.
           It then needs to be linked to a path to write a binary file or text file.
        
        
        Useful:
        * To print restart details, just type: <name of Restart_File>
        * To extract a snapshot: <name of Restart_File>.extract_snapshot(<snapshot id>)
        * To extract a specific material of a snapshot: <name of Restart_File>.extract_material("<name of Material>", <snapshot id>)
        * To extract the evolution of a specific material : <name of Restart_File>.follow_material("<name of Material>")
    """   

    def __init__(self, path_to_file=None, snapshots=None, materials=None):
        """Creates object
        Args: (only one to give)
            path_to_file (str, optional): path to restart file. Defaults to None.
            snapshots (str, optional): list of snapshots. Defaults to None.
            materials (str, optional): dictionary materials for a single snapshot. Defaults to None.
        """        
        self.snapshots = dict()
        test = [not isinstance(i, type(None)) for i in [path_to_file, snapshots, materials]]
        if sum(test) != 1:
            raise Exception('Just choose one option: path, snapshots or materials')
        elif test[0]:
            self.path_to_file = path_to_file
        elif test[1]:
            self.snapshots = snapshots
        elif test[2]:
            self.snapshots[0] = materials

    def __repr__(self):
        """Prints information
        Returns:
            str: Restart information
        """        
        s = 'Restart file: {}, snapshots ({} points):\n'.format(self.path_to_file, len(self.snapshots))
        for i in self.snapshots:
            s += '\t{}: BU = {:.2f} MWd/kg, time  = {:.2f} days, {} materials\n'.format(i, self._burnups[i], self._times[i], len(self.snapshots[i]))
        s += '\tWritten: {}'.format(os.path.exists(self.path_to_file))
        return s

    def read_restart(self):
        """Reads the linked file and creates snapshots.
        """ 
        print('Reading snapshots in {}'.format(self.path_to_file))
        
        # Initialization       
        self.snapshots = dict()
        self._burnups = dict()
        self._times = dict()
        current_step = -1

        # Read restart file
        with open(self.path_to_file, mode='rb') as file:  # b is important -> binary
            while True:
                # Create Material object and fill it with information from one material binary block
                mat = Material()
                read_ok = mat.read(file)

                # If error in reading, end of file
                if not read_ok:
                    break

                # Find the right snapshot to put the material in. 
                # If new, add new snapshot and record time and burnup 
                if len(self._burnups) == 0 or mat.bu_global != self._burnups[list(self._burnups.keys())[-1]]:
                    # Iterate step
                    current_step += 1

                    # New snapshot
                    self.snapshots[current_step] = dict()

                    # Store time/bu information
                    self._burnups[current_step] = mat.bu_global
                    self._times[current_step] = mat.bu_days
                # Store material
                self.snapshots[current_step][mat.name] = mat
        print('\tDone reading: found {} snapshots'.format(len(self.snapshots)))


    def write_binary(self, snapshot_ids=None, material_names=None):
        """Writes snapshots in a restart binary file.
           If no snapshot id is given, all time steps are written.
           Otherwise, it writes time steps corresponding to the snapshot_ids list.
           If no material name is given, all materials are written.
           Otherwise, writes the selected materials.
        Args:
            snapshot_ids (int, optional): List of snapshots to write. Defaults to None.
            material_names (str, optional): dictionary materials to write. Defaults to None.
        """        

        # If no id list is given, take all snapshots
        if isinstance(snapshot_ids, type(None)):
            snapshot_ids = list(self.snapshots.keys())

        print('Writing snapshots to binary {} in {}'.format(snapshot_ids, self.path_to_file))

        # Loop over selected snapshots
        contents = []
        for i_step, s in enumerate([self.snapshots[j] for j in snapshot_ids]):
            # If no material name is given, take all materials
            if isinstance(material_names, type(None)):
                material_names = list(s.keys())
            materials = [s[i] for i in material_names]

            # Loop over materials of the snapshot and add to the content to write
            print('\tProcessing snapshot {} with {} materials'.format(snapshot_ids[i_step], len(material_names)))
            
            for mat in materials:
                contents.append(mat.to_binary())
        
        # Write content
        print('\tWriting ...')
        with open(self.path_to_file, 'wb') as f:
            f.write(b''.join(contents))

        print('\tDone writing: wrote {} snapshots'.format(len(snapshot_ids)))

    def write_text(self, snapshot_ids=None, material_names=None, name_out=None):
        # If no output name is given, take the name/path of the linked binary file
        if isinstance(name_out, type(None)):
            prefix = '.'.join(self.path_to_file.split('.')[:-1])
        else:
            prefix = name_out

        # If no id list is given, take all snapshots
        if isinstance(snapshot_ids, type(None)):
            snapshot_ids = list(self.snapshots.keys())

        print('Writing snapshots {} to text with prefix {}'.format(snapshot_ids, prefix))

        # Loop over selected snapshots
        for i_step, s in enumerate([self.snapshots[j] for j in snapshot_ids]):
            # If no material name is given, take all materials
            if isinstance(material_names, type(None)):
                material_names = list(s.keys())

            print('\tWriting snapshot {} with {} materials'.format(snapshot_ids[i_step], len(material_names)))

            # Path where the file will be written, with snapshot id as suffix
            path_out = '{}_{}.txt'.format(prefix, snapshot_ids[i_step]) 

            # Loop over material and write to snapshot text file
            with open(path_out, 'w') as f:
                for j in material_names:
                    mat = s[j]
                    f.write(mat.to_text())
                    f.write('\n\n')

        print('\tDone writing: wrote {} snapshots'.format(len(snapshot_ids)))

    def write_dataframe(self, snapshot_ids=None, cutoff=-np.inf, material_names=None, name_out=None):
        # If no output name is given, take the name/path of the linked binary file
        if isinstance(name_out, type(None)):
            prefix = '.'.join(self.path_to_file.split('.')[:-1])
        else:
            prefix = name_out

        # If no id list is given, take all snapshots
        if isinstance(snapshot_ids, type(None)):
            snapshot_ids = list(self.snapshots.keys())

        print('Writing snapshots {} to dataframe with prefix {}'.format(snapshot_ids, prefix))

        # Loop over selected snapshots
        for i_step, s in enumerate([self.snapshots[j] for j in snapshot_ids]):
            # If no material name is given, take all materials
            if isinstance(material_names, type(None)):
                material_names = list(s.keys())

            print('\tWriting snapshot {} with {} materials'.format(snapshot_ids[i_step], len(material_names)))

            # Create dataframe
            d = dict()
            for j in material_names:
                mat = s[j]
                d[mat.name] = dict()
                d[mat.name]['Atomic Density [at/b.cm]'] = mat.adens
                d[mat.name]['Mass Density [g/cm^3]'] = mat.mdens
                d[mat.name]['Burnup [MWd/kg_HM]'] = mat.bu
                for nuc in mat.nuclides:
                    d[mat.name][nuc] = mat.nuclides[nuc]['adens']
            dataframe = pd.DataFrame.from_dict(d, orient='index')
            dataframe.columns = list(dataframe.columns[:3]) +  list(translate_ZAI_to_name_many(dataframe.columns[3:]))
            dataframe.index.name = 'Material'
            dataframe = dataframe.loc[natural_sort(dataframe.index)]
            dataframe = dataframe.loc[:,(dataframe>cutoff).all(axis=0)]

            # Path where the file will be written, with snapshot id as suffix
            if np.isinf(cutoff):
                path_out = '{}_{}.csv'.format(prefix, snapshot_ids[i_step]) 
            else:
                path_out = '{}_{}_cut{:.1E}.csv'.format(prefix, snapshot_ids[i_step], cutoff) 

            # Write to snapshot csv file
            dataframe.to_csv(path_out)

        print('\tDone writing: wrote {} snapshots'.format(len(snapshot_ids)))
        return dataframe

    def extract_snapshot(self, snapshot_id):
        """Extracts a specific snapshot based on the id.
        If snapshot is -1, takes the latest added snapshot.
        Args:
            snapshot_id (int): step of the snapshot
        Returns:
            dict: specific snapshot (dictionary Material objects)
        """        
        if snapshot_id == -1:
            snapshot_id = list(self.snapshots.keys())[-1]
        return self.snapshots[snapshot_id]

    def extract_material(self, material_name, snapshot_id):
        """Extracts a specific material for a given snapshot.
        Args:
            material_name (str): name of the material to extract
            snapshot_id (id): step of the snapshot
        Returns:
            Material: extracted material
        """        
        if snapshot_id == -1:
            snapshot_id = list(self.snapshots.keys())[-1]
        return self.snapshots[snapshot_id][material_name]

    def follow_material(self, material_name):
        """Extracts all states for specific material.
        Returns:
            dict: dictionary containing all available material states.
        """               
        states = dict()
        for i in self.snapshots:
            try:
                states[i] = self.snapshots[i][material_name]
            except:
                print('No material {} in snapshot {}. Skipped'.format(material_name, i))
        return states

class Material:
    """Material objects correspond to Serpent materials and include all the informations stored/needed in restart files.
    """    

    def __repr__(self):
        """Prints information
        Returns:
            str: Material information
        """        
        # Extract top inventory
        top = sorted(self.nuclides.items(), key=lambda nuc: nuc[1]['adens'], reverse=True)
        s = '{}, adens: {:.2E}, bu: {:.2E}, nnuc: {}, top 5: {}'.format(self.name, self.adens, self.bu, self.nnuc, ' '.join([translate(i[0]) for i in top[:5]]))
        return s

    def read(self, file):
        """Read binary file block
        Args:
            file (<class '_io.BufferedReader'>): open binary restart file 
        Returns:
            bool: True if the material was successfully read, False otherwise (end of file)
        """        
        # Link to binary file
        self.file_name = file.name

        # Populate material fields from linked binary file
        # Binary files are made of blocks of constant size. 
        # We just need to iterate with the right byte size to read fields

        # Read first sub-block, if not readable, the material is not valid, it is the end of file
        s = file.read(8)
        if not s:
            return False

        # Read snapshot/material fields
        n = struct.unpack("q", s)[0]  # length of material name
        self.name = struct.unpack("{}s".format(n), file.read(n))[0].decode('UTF-8') # material name
        self.bu_global = struct.unpack("d", file.read(8))[0] # BU of snapshot
        self.bu_days = struct.unpack("d", file.read(8))[0] # time of snapshot
        self.nnuc = struct.unpack("q", file.read(8))[0] # Number of nuclides in material
        self.adens = struct.unpack("d", file.read(8))[0] # Atomic density of material
        self.mdens = struct.unpack("d", file.read(8))[0] # Mass density of material
        self.bu = struct.unpack("d", file.read(8))[0] # Burnup of material

        # Read nuclides and populate a dictionary
        self.nuclides = dict()
        for i in range(self.nnuc):
            ZAI, adens = struct.unpack("qd", file.read(16))
            self.nuclides[str(ZAI)] = dict()
            self.nuclides[str(ZAI)]['adens'] = adens

        return True

    def to_binary(self):
        """Converts the material information to a binary block, which can be used for writing a binary restart file
        Returns:
            bytes: material information block
        """        

        # Populate block with necessary material information
        content = b''
        content += struct.pack('q', len(self.name))
        content += struct.pack('{}s'.format(len(self.name)), str.encode(self.name))
        content += struct.pack('d', self.bu_global)
        content += struct.pack('d', self.bu_days)
        content += struct.pack('q', self.nnuc)
        content += struct.pack('d', self.adens)
        content += struct.pack('d', self.mdens)
        content += struct.pack('d', self.bu)
        for i in self.nuclides:
            content += struct.pack('q', int(i))
            content += struct.pack('d', self.nuclides[i]['adens'])
        return content

    def to_text(self):  
        """Converts the material information to a text block.
        
        Returns:
            str: string containing material information
        """   
        s = 'Material {}\n'.format(self.name)
        for k in self.__dict__.keys():
            if k != 'nuclides' and k != 'name' and k != 'file_name':
                s += '\t{}\t{}\n'.format(k, getattr(self, k))
        s += '\tnuclides:\n'
        for k in self.nuclides:
            s += '\t\t{}\t{}\n'.format(k, self.nuclides[k]['adens'])
        return s

    def extract_nuclide(self, name_nuclide):
        if not str(name_nuclide).isdigit():
            name_nuclide = translate(name_nuclide)
        return self.nuclides[name_nuclide]['adens']

    def plot_densities(self, nnuc=None, logscale=True, translating=True):
        """Plot histogram of densities for the top nnuc nuclides
        Args:
            nnuc (int, optional): Number of nuclides to plot. If None, plot all. Defaults to None.
            logscale (bool, optional): If need to use logscale and not linear scale. Defaults to True.
            translating (bool, optional): If need to show real nuclides names. Defaults to True.
        """
        fig, ax = plt.subplots()
        ax.set_axisbelow(True)
        plt.grid()
        if isinstance(nnuc, type(None)):
            nnuc = self.nnuc
        top = sorted(self.nuclides.items(), key=lambda nuc: nuc[1]['adens'], reverse=True)[:nnuc]
        if translating:
            ZAI = [translate(i[0]) for i in top]
        else:
            ZAI = [i[0] for i in top]
        adens = [i[1]['adens'] for i in top]
        ax.bar(ZAI, adens)
        plt.xticks(rotation = 90)
        if logscale:
            plt.yscale('log')
        plt.ylabel('Atomic density [at/b.cm]')
        plt.title(f'{self.name} top {nnuc} nuclide densities' )

#%% Functions
def translate(name):
    """Translate Serpent ZAI notation to human-readable notation for nuclides and vice versa
    Args:
        ZAI (str/int): Serpent ZAI/nuclide name
    Returns:
        str: nuclide name/Serpent ZAI
    """    
    elements = {1: 'H', 2: 'He', 3: 'Li', 4: 'Be', 5: 'B', 6: 'C', 7: 'N', 8: 'O', 9: 'F', 10: 'Ne', 11: 'Na', 12: 'Mg', 13: 'Al', 14: 'Si', 15: 'P', 16: 'S', 17: 'Cl', 18: 'Ar', 19: 'K', 20: 'Ca', 21: 'Sc', 22: 'Ti', 23: 'V', 24: 'Cr', 25: 'Mn', 26: 'Fe', 27: 'Co', 28: 'Ni', 29: 'Cu', 30: 'Zn', 31: 'Ga', 32: 'Ge', 33: 'As', 34: 'Se', 35: 'Br', 36: 'Kr', 37: 'Rb', 38: 'Sr', 39: 'Y', 40: 'Zr', 41: 'Nb', 42: 'Mo', 43: 'Tc', 44: 'Ru', 45: 'Rh', 46: 'Pd', 47: 'Ag', 48: 'Cd', 49: 'In', 50: 'Sn', 51: 'Sb', 52: 'Te', 53: 'I', 54: 'Xe', 55: 'Cs', 56: 'Ba', 57: 'La', 58: 'Ce', 59: 'Pr', 60: 'Nd', 61: 'Pm', 62: 'Sm', 63: 'Eu', 64: 'Gd', 65: 'Tb', 66: 'Dy', 67: 'Ho', 68: 'Er', 69: 'Tm', 70: 'Yb', 71: 'Lu', 72: 'Hf', 73: 'Ta', 74: 'W', 75: 'Re', 76: 'Os', 77: 'Ir', 78: 'Pt', 79: 'Au', 80: 'Hg', 81: 'Tl', 82: 'Pb', 83: 'Bi', 84: 'Po', 85: 'At', 86: 'Rn', 87: 'Fr', 88: 'Ra', 89: 'Ac', 90: 'Th', 91: 'Pa', 92: 'U', 93: 'Np', 94: 'Pu', 95: 'Am', 96: 'Cm', 97: 'Bk', 98: 'Cf', 99: 'Es', 100: 'Fm', 101: 'Md', 102: 'No', 103: 'Lr', 104: 'Rf', 105: 'Db', 106: 'Sg', 107: 'Bh', 108: 'Hs', 109: 'Mt', 110: 'Ds', 111: 'Rg', 112: 'Uub'}
    name = str(name)
    
    # Case where ZAI is given -> name
    if name.isdigit():
        if str(name)[-1] == '0':
            ZA = int(name)/10
            suffix = ''
        elif str(name)[-1] == '1':
            ZA = (int(name) - 1)/10
            suffix = 'm'
        elif str(name)[-1] == '2':
            ZA = (int(name) - 2)/10
            suffix = 'm'
        Z = int(ZA/1000)
        element = elements[Z]
        A = int(ZA-Z*1000)
        if A == 0:
            A = 'nat'
        nuclide = '{}{}{}'.format(element, A, suffix)

    # Case where name is given -> ZAI
    else:
        if name[-1] == 'm':
            name = name[:-1]
            I = '1'
        else:
            I = '0'
        
        i = len(name)-1
        while name[i:].isdigit():
            i-=1
        A = name[i+1:]
        element = name[:i+1]
        Z = str(list(elements.keys())[list(elements.values()).index(element)])
        nuclide = Z+A+I

    return nuclide

def translate_ZAI_to_name_many(names):
    elements = {1: 'H', 2: 'He', 3: 'Li', 4: 'Be', 5: 'B', 6: 'C', 7: 'N', 8: 'O', 9: 'F', 10: 'Ne', 11: 'Na', 12: 'Mg', 13: 'Al', 14: 'Si', 15: 'P', 16: 'S', 17: 'Cl', 18: 'Ar', 19: 'K', 20: 'Ca', 21: 'Sc', 22: 'Ti', 23: 'V', 24: 'Cr', 25: 'Mn', 26: 'Fe', 27: 'Co', 28: 'Ni', 29: 'Cu', 30: 'Zn', 31: 'Ga', 32: 'Ge', 33: 'As', 34: 'Se', 35: 'Br', 36: 'Kr', 37: 'Rb', 38: 'Sr', 39: 'Y', 40: 'Zr', 41: 'Nb', 42: 'Mo', 43: 'Tc', 44: 'Ru', 45: 'Rh', 46: 'Pd', 47: 'Ag', 48: 'Cd', 49: 'In', 50: 'Sn', 51: 'Sb', 52: 'Te', 53: 'I', 54: 'Xe', 55: 'Cs', 56: 'Ba', 57: 'La', 58: 'Ce', 59: 'Pr', 60: 'Nd', 61: 'Pm', 62: 'Sm', 63: 'Eu', 64: 'Gd', 65: 'Tb', 66: 'Dy', 67: 'Ho', 68: 'Er', 69: 'Tm', 70: 'Yb', 71: 'Lu', 72: 'Hf', 73: 'Ta', 74: 'W', 75: 'Re', 76: 'Os', 77: 'Ir', 78: 'Pt', 79: 'Au', 80: 'Hg', 81: 'Tl', 82: 'Pb', 83: 'Bi', 84: 'Po', 85: 'At', 86: 'Rn', 87: 'Fr', 88: 'Ra', 89: 'Ac', 90: 'Th', 91: 'Pa', 92: 'U', 93: 'Np', 94: 'Pu', 95: 'Am', 96: 'Cm', 97: 'Bk', 98: 'Cf', 99: 'Es', 100: 'Fm', 101: 'Md', 102: 'No', 103: 'Lr', 104: 'Rf', 105: 'Db', 106: 'Sg', 107: 'Bh', 108: 'Hs', 109: 'Mt', 110: 'Ds', 111: 'Rg', 112: 'Uub'}

    names = np.array(names).astype(str)
    suffixes = np.empty_like(names)
    stable = np.array([name[-1]=='0' for name in names])
    metastable1 = np.array([(name[-1]=='1') & (name!='-1') for name in names])
    metastable2 = np.array([name[-1]=='2' for name in names])
    suffixes[stable] = ''
    suffixes[(metastable1) | (metastable2)] = 'm'
    ZA = np.array([int(name) for name in names], dtype=int)
    ZA[metastable1] -= 1
    ZA[metastable2] -= 2
    ZA = (ZA/10).astype(int)
    Z = (ZA/1000).astype(int)
    element = np.array([elements[z]  if z!=0 else 'Lost' for z in Z])
    A = (ZA-Z*1000).astype(int).astype(str)
    A[A=='0'] = 'nat'
    A[element=='Lost'] = ''
    nuclides = np.core.defchararray.add(element, A)
    nuclides = np.core.defchararray.add(nuclides, suffixes)
    return nuclides

def natural_sort(l): 
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
    return sorted(l, key=alphanum_key)

#%% Main
if __name__ == "__main__":
    ## Input
    file_in1 = './input.inp.wrk_650' 

    ## Read restart
    restart = Restart_File(path_to_file=file_in1)
    restart.read_restart()

    # Remove parent material "fuel"
    materials = [i for i in restart.snapshots[0].keys() if i != 'fuel']

    ## Quick check
    materials1 = restart.extract_snapshot(0)
    print(f'\tCheck: {len(materials1)} materials in restart file'.format())

    # Create and write dataframe
    df = restart.write_dataframe(snapshot_ids=[0], cutoff=1e-16, material_names=materials)

    # Plot stuff for fun
    materials1['fuel'].plot_densities(nnuc=10, logscale=True, translating=True)
    
    plt.figure()
    nbins=5
    field = 'Xe135'
    grouping = 'Burnup [MWd/kg_HM]'
    df.groupby(pd.cut(df[grouping], bins=nbins)).hist(field, bins=75, ax=plt.gca(),alpha=0.4)
    plt.title(f'{grouping} vs {field}')
    _ = [l.set_alpha(1) for l in plt.legend(range(1, nbins+1), title='BU group').legendHandles]
