import glob
import os
import sys

import matplotlib.pyplot as plt
import matplotlib.transforms as mtransforms
import numpy as np
import template as tp
from netCDF4 import Dataset

# plt.rcParams['figure.constrained_layout.hspace'] = 0
# plt.rcParams['figure.constrained_layout.h_pad'] = 0.0005


# %% Load data
path_data = '../data/output_data'
list_runs = glob.glob(os.path.join(path_data, '*.nc'))
datasets = np.array([Dataset(run) for run in list_runs])

# %% create data vectors
alpha, phi, St, Re, = np.array([[d.variables['alpha'][:].data, d.variables['phi'][:].data,
                                d.variables['St'][:].data, d.variables['Re'][:].data]
                               for d in datasets]).T

authors = np.array([d.author for d in datasets])

# %% graphic specifications
# changing order depending on authors
author_zorder = ['Cyril', 'Cyril/Marie', 'Jean', 'Rastello', 'Julien']

figsize = (tp.half_figure_width, 2.4*tp.half_figure_width)
fig, axarr = plt.subplots(4, 1, figsize=figsize, constrained_layout=True,
                          gridspec_kw={'height_ratios': [0.3, 1, 1, 1]})

for ax, var in zip(axarr[1:].flatten(), [phi, alpha, St]):
    for author in author_zorder:
        mask = (authors == author)
        marker = 's' if author == 'Julien' else None
        ax.scatter(Re[mask], var[mask], c=tp.color_setups[author],
                   label=author, marker=marker)
        ax.set_xscale('log')
        ax.set_yscale('log')

axarr[0].axis('off')
handles, labels = axarr[1].get_legend_handles_labels()
leg = axarr[0].legend(handles, labels, loc="upper center", ncol=3,
                      borderaxespad=0, title='Datasets', mode='expand')
#
axarr[1].set_xticks([])
axarr[2].set_xticks([])
axarr[3].set_xlabel(r'Reynolds number, $\mathcal{R}e$')
#
axarr[1].set_ylabel(r'Volume fraction, $\phi$')
axarr[2].set_ylabel(r'angle, $\alpha~[^\circ]$')
axarr[3].set_ylabel(r'Stokes number, $\mathcal{S}t$')
#
axarr[2].set_yscale('linear')
axarr[2].set_ylim(bottom=-1)

for ax, l in zip(axarr[1:].flatten(), 'abcdefgh'):
    trans = mtransforms.ScaledTranslation(
        5/72, -5/72, fig.dpi_scale_trans)
    label = '({})'.format(l)
    ax.text(0.0, 1.0, label, transform=ax.transAxes + trans, color='k',
            va='top', ha='left')

fig.savefig(
    '../paper/figures/{}.pdf'.format(sys.argv[0].split(os.sep)[-1].replace('.py', '')), dpi=600)
