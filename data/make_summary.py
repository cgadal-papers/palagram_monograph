import numpy as np
import glob
import os
from netCDF4 import Dataset

# %%
path_data = '/output_data'
list_runs = sorted(glob.glob(os.path.join(path_data, '*.nc')))
