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
author_zorder = ['Cyril', 'Cyril/Marie', 'Jean', 'Rastello', 'Julien']

# %% masks for plot
alphas = [0, 45]
alpha_pad = 1

fig, axarr = plt.subplots(2, 1, constrained_layout=True,
                          figsize=tp.half_figure_size_inv, sharex=True)

for alpha0, ax in zip(alphas, axarr.flatten()):
    mask_alpha = (alpha > alpha0 - alpha_pad) & (alpha < alpha0 + alpha_pad)
    if alpha0 == 0:
        axins = ax.inset_axes([0.43, 0.53, 0.55, 0.43])
        axins.set_ylim(0, 0.8)
        axins.set_xlim(right=0.8)
        axins.set_xlabel('$\phi$', labelpad=0)
        axins.set_ylabel(r'$\mathcal{F}r$', labelpad=2)
    for author in author_zorder:
        mask = (authors == author) & mask_alpha
        marker = 's' if author == 'Julien' else None
        #
        ax.scatter(phi[mask], Fr[mask], c=tp.color_setups[author],
                   label=author, marker=marker)
        if alpha0 == 0:
            axins.scatter(phi[mask], Fr[mask], c=tp.color_setups[author],
                          label=author, marker=marker)
    # saline current
    mask_saline = (particle_type == 'saline water') & mask_alpha
    moy, std = np.nanmean(Fr[mask_saline]), np.nanstd(Fr[mask_saline])
    ax.axhline(moy, color='k', ls='--', zorder=-10)
    ax.axhspan(moy - std, moy + std, color='k', alpha=0.2, zorder=-10)
    #
    ax.set_ylabel(r'Froude number, $\mathcal{F}r$')
    ax.set_ylim([0, 1.4])

ax.set_xlabel(r'Volume fraction, $\phi$')
ax.set_xscale('log')

for ax, l in zip(axarr.flatten(), ['a', 'b']):
    trans = mtransforms.ScaledTranslation(
        5/72, -5/72, fig.dpi_scale_trans)
    label = '({})'.format(l)
    ax.text(0.0, 1.0, label, transform=ax.transAxes + trans, color='k',
            va='top', ha='left')

fig.savefig(
    '../paper/figures/{}.pdf'.format(sys.argv[0].split(os.sep)[-1].replace('.py', '')), dpi=600)
