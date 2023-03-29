import glob
import os

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
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
        ax.scatter(St[mask], var[mask], c=tp.color_setups[author],
                   label=author, marker=marker)
        #
# saline current
# ATTENTION Saline all NAn -> check fits ?
mask_saline = (particle_type == 'saline water') & mask_alpha
moy, std = np.nanmean(L[mask_saline]), np.nanstd(L[mask_saline])
axarr[1].axhline(moy, color='k', ls='--', zorder=-10)
axarr[1].axhspan(moy - std, moy + std, color='k', alpha=0.2, zorder=-10)

axarr[1].set_xlabel(r'Stokes number, $\mathcal{S}t$')
axarr[1].set_xscale('log')

axarr[0].set_ylabel(r'Froude number, $\mathcal{F}r$')
axarr[1].set_ylabel(r'Non dim. dissipation, $\tilde{\lambda}$')
axarr[0].set_ylim(0, 1.4)
# axarr[1].set_ylim(0, 1.4)

plt.show()
