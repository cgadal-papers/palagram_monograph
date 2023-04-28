import glob
import os
import sys

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import template as tp
from netCDF4 import Dataset
from lmfit.model import load_modelresult


# %% Load data
path_data = '../data/output_data'
list_runs = sorted(glob.glob(os.path.join(path_data, '*.nc')))
list_fitresults = sorted(glob.glob(os.path.join(path_data, 'fitresult*')))

datasets = np.array([Dataset(run) for run in list_runs])

# %% mask data
alpha, phi, rho_p = np.array(
    [[d.variables['alpha'][:].data, d.variables['phi'][:].data, d.variables['rho_p'][:].data]
     for d in datasets]).T

# mask = (phi > 0.001) & (phi < 0.6) & (alpha < 70) & (rho_p > 1000)
mask = np.ones_like(alpha).astype('bool')

zorder_setups = {'Cyril': -9, 'Rastello': -10,
                 'Jean': -6, 'Julien': -7, 'Cyril/Marie': -9}

# %% figure
layout = [['legend', 'legend'], [('(a)'), '(b)']]

fig, axarr = plt.subplot_mosaic(layout, figsize=tp.large_figure_size, constrained_layout=True,
                                gridspec_kw={'height_ratios': [0.001, 1]})
# All runs
ax = axarr['(a)']
for d in datasets[mask]:
    ax.scatter(d.variables['t'][:].data, d.variables['x_front'][:].data,
               s=2 if d.author == 'Julien' else 0.3, color=tp.color_setups[d.author],
               zorder=zorder_setups[d.author], rasterized=True)

ax.set_ylabel(r'Front position, $x_{\rm f}$ [m]')
ax.set_xlabel(r'Time, $t$ [s]')
ax.set_xlim(0, 100)
ax.set_ylim(bottom=0)

# selected run, non-dimensional
ax = axarr['(b)']
runs = [50, 100, 150, 200, 250]

for run in runs:
    d = datasets[list_runs.index(os.path.join(
        path_data, 'run_{:03d}.nc'.format(run)))]
    print(d.author)
    #
    t_ad = d.variables['L0'][:].data/d.variables['u0'][:].data
    x_ad = d.variables['L0'][:].data
    #
    x_axis = d.variables['t'][:].data/t_ad
    y_axis = d.variables['x_front'][:].data/x_ad
    ax.scatter(x_axis, y_axis, color=tp.color_setups[d.author], s=2)
    # plotting fit
    fitresult = load_modelresult(os.path.join(
        path_data, 'fitresult_run_{:03d}.save'.format(run)))
    xplot = np.linspace(
        fitresult.userkws['t'].min(), fitresult.userkws['t'].max(), 500)
    ax.plot(xplot, fitresult.eval(t=xplot),
            color='k', lw=1, ls='--')

ax.set_ylabel(r'Front position, $x_{\rm f}/L_{0}$')
ax.set_xlabel(r'Time, $t/t_{0}$')

# axarr['legend'].axis('off')
# handles, labels = axarr['(a)'].get_legend_handles_labels()
# leg = axarr['legend'].legend(handles, labels, loc="upper center", ncol=5,
#                              borderaxespad=0, title='Datasets')

fig.savefig(
    '../paper/figures/{}.pdf'.format(sys.argv[0].split(os.sep)[-1].replace('.py', '')), dpi=600)
