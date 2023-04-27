import glob
import os
import sys

import matplotlib.pyplot as plt
import matplotlib.transforms as mtransforms
import numpy as np
import template as tp
from netCDF4 import Dataset

plt.rcParams['figure.constrained_layout.hspace'] = 0
plt.rcParams['figure.constrained_layout.h_pad'] = 0

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
alpha0 = [3, 7, 15, 45]
alpha_pad = 1.5

figsize = (tp.large_figure_width, tp.golden*tp.large_figure_width/2)
fig, axarr = plt.subplots(4, 2, constrained_layout=True,
                          figsize=figsize, sharex=True)
for a0, axarr_sub in zip(alpha0, axarr):
    mask_alpha = (alpha > a0 - alpha_pad) & (alpha < a0 + alpha_pad)
    for i, (var, ax) in enumerate(zip([Fr, L], axarr_sub.flatten())):
        for author in author_zorder:
            mask = (authors == author) & mask_alpha
            marker = 's' if author == 'Julien' else None
            #
            ax.scatter(St[mask]/a[mask], var[mask], c=tp.color_setups[author],
                       label=author, marker=marker)
        #
        if a0 == 7:
            # saline current
            mask_saline = (particle_type == 'saline water') & mask_alpha
            moy, std = np.nanmean(
                var[mask_saline]), np.nanstd(var[mask_saline])
            ax.axhline(moy, color='k', ls='--', zorder=-10)
            ax.axhspan(moy - std, moy + std, color='k', alpha=0.2, zorder=-10)
        else:
            if i == 0:
                moy, std = np.nanmean(
                    var[mask_alpha]), np.nanstd(var[mask_alpha])
                ax.axhline(moy, color='k', ls='--', zorder=-10)
                ax.axhspan(moy - std, moy + std, color='k',
                           alpha=0.2, zorder=-10)
            else:
                ax.axhline(0, color='k', ls=':', zorder=-10, lw=1)

    axarr_sub[0].set_ylabel(r'$\sqrt{a}\mathcal{F}r$')
    axarr_sub[1].set_ylabel(r'Dissipation, $\tilde{\lambda}$')

    axarr_sub[1].set_xscale('log')

    # axarr_sub[0].set_ylim(0, 1.85)
    axarr_sub[0].set_ylim(0, 1.4)
    axarr_sub[1].set_ylim(-0.02, 0.07)
    axarr_sub[1].set_xlim(1.5e-3, 1)
    axarr_sub[0].set_xlim(1.5e-3, 1)
    #
    if a0 != 45:
        axarr_sub[0].set_xticklabels([])
        axarr_sub[1].set_xticklabels([])
    else:
        axarr_sub[1].set_xlabel(r'Stokes number, $\mathcal{S}t/a$')
        axarr_sub[0].set_xlabel(r'Stokes number, $\mathcal{S}t/a$')


for ax, l in zip(axarr.flatten(), 'abcdefgh'):
    trans = mtransforms.ScaledTranslation(
        5/72, -5/72, fig.dpi_scale_trans)
    label = '({})'.format(l)
    ax.text(0.0, 1.0, label, transform=ax.transAxes + trans, color='k',
            va='top', ha='left')

fig.savefig(
    '../paper/figures/{}.pdf'.format(sys.argv[0].split(os.sep)[-1].replace('.py', '')), dpi=600)
