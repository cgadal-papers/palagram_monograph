import numpy as np
import glob
import os
from netCDF4 import Dataset

# %%
path_data = './output_data'
list_runs = sorted(glob.glob(os.path.join(path_data, '*.nc')))

for i, run in enumerate(list_runs):
    d = Dataset(run)
    attributes = sorted(d.__dict__.keys())
    variables = [var for var in d.variables.keys() if var not in [
        'x_front', 't', 'v_front']]
    if i == 0:  # writing headline
        column_heads = (attributes +
                        ['{} ({})'.format(par, d.variables[par].unit)
                            for par in variables]
                        )
        line = ','.join(column_heads) + '\n'
        with open('dataset_summary.csv', 'w') as the_file:
            the_file.write(line)
            for i in range(2):
                the_file.write('\n')
    # writing data to csv
    line_values = ([d.__dict__[attr] for attr in attributes] +
                   ['{:.3e}'.format(d.variables[var][:].data)
                    for var in variables]
                   )
    line = ','.join(line_values) + '\n'
    with open('dataset_summary.csv', 'a') as the_file:
        the_file.write(line)

[d.variables[var][:].data for var in variables if var not in [
    'x_front', 't']]
