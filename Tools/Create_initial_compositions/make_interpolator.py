import serpentTools as st
from scipy.interpolate import interp1d
import pickle

depletion_output_file = '/home/yryves/serpent_cases/single_pebble/input_dep.m'
fuel_material = 'fuel'

dep = st.read(depletion_output_file)
mat = dep.materials[fuel_material]
list_iso = mat.zai
interpolator = interp1d(mat.burnup, mat.adens, kind="linear")
with open("{}.interpolator".format(depletion_output_file.split(".m")[0]), "wb") as file_handle:
    pickle.dump((interpolator, list_iso), file_handle)
