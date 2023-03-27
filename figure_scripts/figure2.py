import glob
import os
import sys

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import template as tp
from matplotlib.colors import to_rgba
from matplotlib.path import Path
from netCDF4 import Dataset


# %% Load data
path_data = '../data/output_data'
list_runs = glob.glob(os.path.join(path_data, '*.nc'))
datasets = np.array([Dataset(run) for run in list_runs])
color_setups = {'Cyril': -9, 'Rastello': -10,
                'Jean': -6, 'Julien': -7, 'Cyril/Marie': -9}

# %% mask data
phi = np.array([d.variables['phi'][:].data for d in datasets])
alpha = np.array([d.variables['alpha'][:].data for d in datasets])
rho_p = np.array([d.variables['rho_p'][:].data for d in datasets])

# mask = (phi > 0.001) & (phi < 0.6) & (alpha < 70) & (rho_p > 1000)
mask = np.ones_like(alpha).astype('bool')

# %% figure
fig, ax = plt.subplots(1, 1, constrained_layout=True,
                       figsize=tp.large_figure_size)
for d in datasets[mask]:
    ax.scatter(d.variables['t'][:].data, d.variables['x_front'][:].data,
               s=2 if d.author == 'Julien' else 0.3, color=tp.color_setups[d.author],
               zorder=color_setups[d.author], rasterized=True)

ax.set_ylabel('Front position, $x_{f}$ [m]')
ax.set_xlabel('Time, $t$ [s]')
ax.set_xlim(0, 100)
ax.set_ylim(bottom=0)

legend_elements = [Line2D([0], [0], marker='.', ls='none', color=c, label=key)
                   for key, c in sorted(tp.color_setups.items())]
ax.legend(handles=legend_elements)

fig.savefig(
    '../paper/figures/{}.pdf'.format(sys.argv[0].split(os.sep)[-1].replace('.py', '')), dpi=600)
