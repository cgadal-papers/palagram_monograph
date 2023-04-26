import glob
import os
import sys

import matplotlib.pyplot as plt
import matplotlib.transforms as mtransforms
import numpy as np
import template as tp
from netCDF4 import Dataset

# %% Load data
path_data = '../data/output_data'
list_runs = glob.glob(os.path.join(path_data, '*.nc'))
datasets = np.array([Dataset(run) for run in list_runs])

# %% create data vectors
alpha, phi, Fr, St, L = np.array([[d.variables['alpha'][:].data, d.variables['phi'][:].data,
                                   d.variables['Fr'][:].data, d.variables['St'][:].data,
                                   d.variables['L'][:].data]
                                  for d in datasets]).T

H0 = np.array([d.variables['H0'][:].data for d in datasets])
L0 = np.array([d.variables['L0'][:].data for d in datasets])
a = H0/L0

authors = np.array([d.author for d in datasets])
particle_type = np.array([d.particle_type for d in datasets])

# %% graphic specifications
# changing order depending on authors
author_zorder = ['Cyril', 'Cyril/Marie', 'Julien', 'Rastello', 'Jean']

# %% masks for plot
alpha0 = 7
alpha_pad = 1.5
mask_alpha = (alpha > alpha0 - alpha_pad) & (alpha < alpha0 + alpha_pad)

fig, axarr = plt.subplots(2, 1, constrained_layout=True,
                          figsize=tp.half_figure_size_inv, sharex=True)
for var, ax in zip([Fr, L], axarr.flatten()):
    for author in author_zorder:
        mask = (authors == author) & mask_alpha
        marker = 's' if author == 'Julien' else None
        #
        ax.scatter(St[mask]/a[mask], var[mask], c=tp.color_setups[author],
                   label=author, marker=marker)
    # saline current
    mask_saline = (particle_type == 'saline water') & mask_alpha
    moy, std = np.nanmean(var[mask_saline]), np.nanstd(var[mask_saline])
    ax.axhline(moy, color='k', ls='--', zorder=-10)
    ax.axhspan(moy - std, moy + std, color='k', alpha=0.2, zorder=-10)

axarr[1].set_xlabel(r'Stokes number, $\mathcal{S}t/a$')
axarr[1].set_xscale('log')
# axarr[1].set_yscale('asinh', linear_width=0.006)

axarr[0].set_ylabel(r'Froude number, $\mathcal{F}r$')
axarr[1].set_ylabel(r'Non dim. dissipation, $\tilde{\lambda}$')
axarr[0].set_ylim(0, 1.4)
# axarr[1].set_ylim(0, 1.4)

for ax, l in zip(axarr.flatten(), ['a', 'b']):
    trans = mtransforms.ScaledTranslation(
        5/72, -5/72, fig.dpi_scale_trans)
    label = '({})'.format(l)
    ax.text(0.0, 1.0, label, transform=ax.transAxes + trans, color='k',
            va='top', ha='left')

fig.savefig(
    '../paper/figures/{}.pdf'.format(sys.argv[0].split(os.sep)[-1].replace('.py', '')), dpi=600)
