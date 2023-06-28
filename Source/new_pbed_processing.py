import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import re
from glob import glob

from matplotlib.patches import Circle
from matplotlib.collections import PatchCollection
from matplotlib.colors import Normalize, LogNorm
import matplotlib


def natural_sort(l): 
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
    return sorted(l, key=alphanum_key)

def Z_to_element(Z):
    elements = ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt']
    return elements[Z-1]

def element_to_Z(element):
    elements = ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt']
    return elements.index(element)+1

def ZA_to_name(ZA):
    ZA = str(ZA)
    # Case where ZA (+ temperature card) -> name
    ZA = int(ZA.split('.')[0])
    Z = int(ZA/1000)
    A = int(ZA-Z*1000)
    if A>=300:
        suffix = 'm'
        if Z > 80:
            A -= 100
        else:
            A -= 200
    else:
        suffix = ''
    if A == 0:
        A = 'nat'
    element = Z_to_element(Z)
    name = '{}{}{}'.format(element, str(A).zfill(3), suffix)
    return name

def ZAI_to_name(ZAI):
    # Case where ZAI is given -> name
    if str(ZAI)[-1] == '0':
        ZA = int(ZAI)/10
        suffix = ''
    else:
        ZA = (int(ZAI) - int(str(ZAI)[-1]))/10
        suffix = 'm'

    Z = int(ZA/1000)
    if Z > 120:
        print('nan')
    element = Z_to_element(Z)
    A = int(ZA-Z*1000)
    if A == 0:
        A = 'nat'
    name = '{}{}{}'.format(element, str(A).zfill(3), suffix)
    return name

def name_to_ZAI(name, separated=False):
    # Case where name is given -> ZAI
    if name[-1] == 'm':
        name = name[:-1]
        I = '1'
    else:
        I = '0'

    if name[-3:] == 'nat':
        element = name[:-3]
        Z = str(element_to_Z(element))
        if not separated:
            return int(Z+'000')
        else:
            return int(Z), 0

    i = len(name)-1
    while name[i:].isdigit():
        i-=1
    A = name[i+1:].zfill(3)
    element = name[:i+1]
    Z = str(element_to_Z(element))
    ZAI = int(Z+A+I)
    if not separated:
        return int(ZAI)
    else:
        return int(Z),int(A),int(I)

def name_to_ZA(name, separated=False):
    # Case where name is given -> ZA
    if name[-1] == 'm':
        name = name[:-1]
        isomer = True
    else:
        isomer = False
    if name[-3:] == 'nat':
        element = name[:-3]
        Z = str(element_to_Z(element))
        if not separated:
            return int(Z+'000')
        else:
            return int(Z), 0

    i = len(name)-1
    while name[i:].isdigit():
        i-=1
    A = name[i+1:].zfill(3)
    element = name[:i+1]
    Z = str(element_to_Z(element))
    if isomer:
        if int(A) < 100:
            A = str(int(A)+200)
        else:
            A = str(int(A)+100)
    ZA = Z+A
    if not separated:
        return int(ZA)
    else:
        return int(Z),int(A)


def last_color(last_type='line'):
    if last_type=='line':
        return plt.gca().get_lines()[-1].get_color()

def slice_pbed(data, column='y', val=0.0, tol=None):
    if isinstance(tol, type(None)):
        sliced_data = data[data[column] == val]
    else:
        sliced_data = data[(data[column]-val).abs()<=tol]
    return sliced_data

def projection(data, column, val):
    rel_pos = val - data[column]
    tmp = data['r']**2 - rel_pos**2
    r_projected = np.ones_like(tmp)*np.nan
    r_projected[tmp >= 0] = tmp[tmp >= 0]**0.5
    return r_projected

def mask(data, values, field, direction, vals_mask):
    try:
        _ = len(vals_mask)
    except:
        vals_mask = [vals_mask]

    mask_vals = np.ones_like(values) * np.nan
    for val in vals_mask:
        if direction == '=':
            mask = data[field] == val
        elif direction == '>':
            mask = data[field] >  val
        elif direction == '>=':
            mask = data[field] >= val
        elif direction == '<':
            mask = data[field] <  val
        elif direction == '<=':
            mask = data[field] <= val
        mask_vals[mask] = values[mask]
    return mask_vals

def convert_fima(data, ini_HM, percent):
    hm_iso = [90227, 90228, 90229, 90230, 90232, 90233, 90234, 91231, 91232, 91233, 92232, 92233, 92234, 92235, 92236, 92237, 92238, 92239, 92240, 92241, 93235, 93236, 93237, 93238, 93239, 94236, 94237, 94238, 94240, 94239, 94241, 94242, 94243, 94244, 95241, 95242, 95243, 95244, 96240, 96241, 96242, 96243, 96244, 96245, 96246, 96247, 96248, 96249, 96250, 902270, 902280, 902290, 902300, 902320, 902330, 902340, 912310, 912320, 912330, 922320, 922330, 922340, 922350, 922360, 922370, 922380, 922390, 922400, 922410, 932350, 932360, 932370, 932380, 932390, 942360, 942370, 942380, 942400, 942390, 942410, 942420, 942430, 942440, 952410, 952420, 952430, 952440, 962400, 962410, 962420, 962430, 962440, 962450, 962460, 962470, 962480, 962490, 962500]
    hm_found = []
    for i in hm_iso:
        if i in data.columns:
            hm_found.append(i)
    data['HM'] = data[hm_found].sum(axis=1)
    data['fima'] = (ini_HM-data['HM'])/ini_HM
    if percent:
        data['fima'] *= 100
    return data

def plot_slice(data, field='id', column='y', val=0.0, tol=None, project=True, clim=None, log_color=False, colormap='turbo', field_mask=None, field_dir='=', field_vals=None, nan_col=[0.8,0.8,0.8], clim_rule='global'):
    if column == 'x':
        xdir = 'y'
        ydir = 'z'
        column = 'x'
    elif column == 'y':
        xdir = 'x'
        ydir = 'z'
        column = 'y'
    elif column == 'z':
        xdir = 'x'
        ydir = 'y'
        column = 'z'    
    
    if isinstance(tol, type(None)):
        tol = data['r']
    sliced_data = slice_pbed(data, column, val, tol)

    x = np.array(sliced_data[xdir])
    y = np.array(sliced_data[ydir])
    if project:
        r = projection(sliced_data, column, val)
    else:
        r = np.array(sliced_data['r'])
        
    colors = np.array(sliced_data[field])
    if not isinstance(field_mask, type(None)):
        colors = mask(sliced_data, colors, field_mask, field_dir, field_vals)
    patches = []
    for i in range(len(sliced_data)):
        circle = Circle((x[i], y[i]), r[i])
        patches.append(circle)
        
    if isinstance(clim, type(None)):
        if clim_rule=='mask':
            clim = [np.nanmin(colors), np.nanmax(colors)]
        elif clim_rule=='slice':
            clim = [np.nanmin(sliced_data[field]), np.nanmax(sliced_data[field])]
        elif clim_rule=='global':
            clim = [np.nanmin(data[field]), np.nanmax(data[field])]
            
    cmap = matplotlib.colormaps[colormap]
    cmap.set_bad(nan_col, alpha = 1.) # nans will appear gray
    if log_color:
        p = PatchCollection(patches, cmap=cmap, norm=LogNorm())
    else:
        p = PatchCollection(patches, cmap=cmap)

    p.set_array(colors)
    p.set_clim(clim)   
    ax = plt.gca()
    ax.add_collection(p)
    ax.autoscale_view()
    ax.set_aspect("equal", adjustable="box")
    plt.tight_layout()
    return ax, p

def plot_profile(data, field='id', column='y', plot_type='mean', mean_bins=None, dir_slice=None, val_slice=None, tol=None, field_mask=None, field_dir='=', field_vals=None, label='', linecolor=None, linestyle='-', markerstyle=None, markersize=None, linewidth=None, alpha=0.4):
    if isinstance(dir_slice, type(None)):
        sliced_data = data
    else:
        sliced_data = slice_pbed(data, dir_slice, val_slice, tol)
    if isinstance(linecolor, type(None)):
        color_cycle = plt.gca()._get_lines.prop_cycler
        linecolor = next(color_cycle)["color"]
    if plot_type in ['mean', 'std'] and not isinstance(mean_bins, type(None)):
        groups = sliced_data.groupby(pd.cut(sliced_data[column], bins=mean_bins))
        X = [x.mid for x in groups[field].apply(plot_type).index]
    else:
        groups = sliced_data.groupby(column)
        X = groups.groups.keys()
    if plot_type=='std':
        mean = groups[field].mean()
        std = groups[field].std()
        plt.fill_between(X, mean-std, mean+std, alpha=alpha, label=label, color=linecolor)
    else:
        plt.plot(X, groups[field].apply(plot_type), label=label, color=linecolor, linewidth=linewidth, ls=linestyle, marker=markerstyle, markersize=markersize)

def get_dist(data, field='id', bins=50, norm=True, scaling_factor=1.0):
    counts, bins = np.histogram(data[field], bins=bins)
    if norm:
        counts = counts.astype(float)/np.nanmax(counts)
    counts =  counts.astype(float)*scaling_factor
    return counts, bins

def plot_dist(data, field='id', bins=50, binwidth=1, norm=True, scaling_factor=1.0, dir_slice=None, val_slice=None, tol=None, label='', show_mean=False, mean_linestyle='--', mean_color=None, mean_linewidth=None, color=None, linewidth=None, edgecolor=None, alpha=0.4):
    if isinstance(dir_slice, type(None)):
        sliced_data = data
    else:
        sliced_data = slice_pbed(data, dir_slice, val_slice, tol)
    if isinstance(color, type(None)):
        color_cycle = plt.gca()._get_lines.prop_cycler
        color = next(color_cycle)["color"]
    counts, bins = get_dist(sliced_data, field, bins, norm, scaling_factor)
    plt.bar((bins[1:] + bins[:-1])/2, counts, width=np.diff(bins)*binwidth, label=label, alpha=alpha, color=color, linewidth=linewidth, edgecolor=edgecolor)
    if show_mean:
        if isinstance(mean_color, type(None)):
            mean_color = color
        plt.axvline(sliced_data[field].mean(), linestyle=mean_linestyle, color=mean_color, linewidth=mean_linewidth)
    
def plot_zones(zoned_data, field, pass_number):
    data_pass = zoned_data[zoned_data['passes']==pass_number].sort_values(['radial_ID', 'axial_ID'])
    values = np.empty((Nr, Nz))
    for r in range(Nr):
        for z in range(Nz):
            values[r, z] = data_pass[(data_pass['radial_ID']==r+1) & (data_pass['axial_ID']==z+1)][field]
    fig, ax = plt.subplots()
    ax.set_xlabel('R')
    ax.set_ylabel('Z')
    ax.set_xticks(range(1, Nr+1))
    ax.set_yticks(range(1, Nz+1))
    pcolor = ax.pcolor(range(1, Nr+1), range(1, Nz+1), values.T, cmap='turbo')
    plt.colorbar(pcolor)

def determine_zones(core_data, Nr, Nz, Rmax, Zmin, Zmax):
    core_data['radial_ID'] = pd.cut(core_data['r_dist'], bins=np.linspace(0, Rmax, Nr+1), include_lowest=True, labels=False).astype(int) + 1
    core_data['axial_ID'] = pd.cut(core_data['z'], bins=np.linspace(Zmin, Zmax, Nz+1), include_lowest=True, labels=False).astype(int) + 1
    return core_data

def zone_pbed(core_data, Nr, Nz, Rmax, Zmin, Zmax):
    core_data = determine_zones(core_data, Nr, Nz, Rmax, Zmin, Zmax)
    zoned_pbed = core_data.groupby(['axial_ID', 'radial_ID', 'passes']).mean().reset_index()
    return zoned_pbed

def plot_evolution(cycle_data, field, field_vs='passes', field_err=None, rel_err=True, multiplier_err=1.0, alpha=0.4, savefig=None, dpi=300):
    x = cycle_data[field_vs]
    y = cycle_data[field]
    plt.plot(x, y)
    if not isinstance(field_err, type(None)):
        error = cycle_data[field_err]
        if rel_err:
            error *= y
        plt.fill_between(x, y-error*multiplier_err, y+error*multiplier_err, color=last_color(), alpha=alpha)
    plt.tight_layout()
    if not isinstance(savefig, type(None)):
        if savefig[-4:] != '.png':
            savefig = savefig + '.png'
        plt.savefig(savefig, bbox_inches='tight', dpi=dpi)

def plot_core_fields(core_data, fields, num_cols=4, savefig=None, dpi=200):
    num_fields = len(fields)
    num_rows = max(1, (num_fields + num_cols - 1) // num_cols)  # Calculate the number of rows
    if num_fields > 30:
        return
    fig, axes = plt.subplots(num_rows, num_cols, figsize=(3*num_cols, 6*num_rows))
    for i, field in enumerate(list(fields)):
        plt.sca(axes.flatten()[i])
        ax, p = plot_slice(core_data, field)
        cbar = plt.colorbar(p, shrink=0.8, orientation='horizontal')
        cbar.formatter.set_powerlimits((0, 0))
        plt.axis('off')
        title = '_'.join(field.split('_'))  # Splitting the title string by underscores
        title_words = title.split('_')
        title_with_newlines = '\n'.join(' '.join(title_words[i:i+3]) for i in range(0, len(title_words), 3))
        plt.title(title_with_newlines)

    if len(fields) < num_cols*num_rows:
        [ax.set_visible(False) for ax in axes.flatten()[len(fields):num_cols*num_rows]]
    plt.tight_layout()
    if not isinstance(savefig, type(None)):
        if savefig[-4:] != '.png':
            savefig = savefig + '.png'
        plt.savefig(savefig, bbox_inches='tight', dpi=dpi)