import numpy as np
import glob
import os
from netCDF4 import Dataset

# %%
path_data = './output_data'
list_runs = sorted(glob.glob(os.path.join(path_data, '*.nc')))

file = '../dataset_summary.csv'

for i, run in enumerate(list_runs):
    d = Dataset(run)
    attributes = sorted(d.__dict__.keys())
    variables = [var for var in d.variables.keys() if var not in [
        'x_front', 't', 'v_front']]
    if i == 0:  # writing headline
        var_save = np.copy(variables)
        column_heads = (['file_ID'] + attributes +
                        ['{} ({})'.format(par, d.variables[par].unit)
                            for par in variables]
                        )
        line = ','.join(column_heads) + '\n'
        with open(file, 'w') as the_file:
            the_file.write(line)
            for i in range(2):
                the_file.write('\n')
    # writing data to csv
    line_values = ([run.split(os.sep)[-1]] + [d.__dict__[attr] for attr in attributes] +
                   ['{:.3e}'.format(d.variables[var][:].data) if (var in variables) else '-'
                    for var in var_save]
                   )
    line = ','.join(line_values) + '\n'
    with open(file, 'a') as the_file:
        the_file.write(line)
