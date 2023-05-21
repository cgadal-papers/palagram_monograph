import glob
import os
import sys

import matplotlib.pyplot as plt
import matplotlib.transforms as mtransforms
import numpy as np
import template as tp
from netCDF4 import Dataset

plt.rcParams['figure.constrained_layout.hspace'] = 0
plt.rcParams['figure.constrained_layout.h_pad'] = 0.0005

# %% Load data
path_data = '../data/output_data'
list_runs = glob.glob(os.path.join(path_data, '*.nc'))
datasets = np.array([Dataset(run) for run in list_runs])

# %% create data vectors
alpha, phi, St, H0, L0, Fr, L = np.array([[d.variables['alpha'][:].data, d.variables['phi'][:].data,
                                           d.variables['St'][:].data, d.variables['H0'][:].data,
                                           d.variables['L0'][:].data, d.variables['Fr'][:].data,
                                           d.variables['L'][:].data
                                           ] for d in datasets]).T

Ha = np.array([d.variables['H_a'][:].data if 'H_a' in d.variables.keys()
              else d.variables['H0'][:].data for d in datasets])
authors, particles = np.array([[d.author, d.particle_type,
                                ] for d in datasets]).T
a = H0/L0

# %% graphic specifications
dataset_idx = np.vectorize(tp.datasets.get)(authors)

markers = np.vectorize(tp.marker_style.get)(particles)
markers[dataset_idx == 'SedFoam'] = 's'

facecolors = np.vectorize(tp.color_datasets.get)(dataset_idx)
edgecolors = np.full_like(facecolors, 'k')

mask_nosuspended = (authors == 'Rastello') & (H0/Ha < 1)
edgecolors[mask_nosuspended] = 'tab:red'

zorders = np.vectorize(lambda dataset: tp.datset_zorder[dataset])(dataset_idx)
random_order = np.arange(zorders.size)
rng = np.random.default_rng(1994)
rng.shuffle(random_order)
plot_idxs = np.lexsort((random_order, zorders))

# %% masks for plot
alpha0 = [0, 7, 15, 45]
alpha_pad = 1.5

figsize = (tp.large_figure_width, tp.golden*tp.large_figure_width/2)
fig, axarr = plt.subplots(4, 3, constrained_layout=True,
                          figsize=figsize, gridspec_kw={'width_ratios': [0.2, 1, 1]})
for a0, axarr_sub in zip(alpha0, axarr[:, 1:]):
    mask_alpha = (alpha > a0 - alpha_pad) & (alpha < a0 + alpha_pad)
    mask = (mask_alpha & (phi < 0.449))[plot_idxs]
    for i, (var, ax) in enumerate(zip([Fr, L], axarr_sub.flatten())):
        tp.mscatter((St/a)[plot_idxs][mask], var[plot_idxs][mask], ax=ax, m=markers[plot_idxs][mask],
                    facecolors=facecolors[plot_idxs][mask], edgecolors=edgecolors[plot_idxs][mask], lw=0.5)
        #
        if i == 0:
            moy, std = np.nanmean(
                var[plot_idxs][mask]), np.nanstd(var[plot_idxs][mask])
            ax.axhline(moy, color='k', ls=':', zorder=-10, lw=1)
            # ax.axhspan(moy - std, moy + std, color='k', alpha=0.2, zorder=-10)
        else:
            ax.axhline(0, color='k', ls=':', zorder=-10, lw=1)
        #
        if a0 == 7:
            mask_saline = ((particles == 'saline water')
                           & mask_alpha)
            moy, std = np.nanmean(
                var[mask_saline]), np.nanstd(var[mask_saline])
            ax.axhline(moy, color='k', ls='--', zorder=-10, lw=1)
            # ax.axhspan(moy - std, moy + std, color='k', alpha=0.2, zorder=-10)

    axarr_sub[0].set_ylabel(r'Froude number, $\mathcal{F}r$')
    axarr_sub[1].set_ylabel(r'Attenuation, $\tilde{\lambda}$')

    axarr_sub[1].set_xscale('log')
    axarr_sub[0].set_xscale('log')

    axarr_sub[0].set_ylim(0, 1.59)
    axarr_sub[1].set_ylim(-0.02, 0.07)
    axarr_sub[1].set_xlim(1.5e-3, 1)
    axarr_sub[0].set_xlim(1.5e-3, 1)
    #
    axarr_sub[1].set_yticks([0, 0.03, 0.06])
    #
    if a0 != 45:
        axarr_sub[0].set_xticklabels([])
        axarr_sub[1].set_xticklabels([])
    else:
        axarr_sub[1].set_xlabel(r'Stokes number, $\mathcal{S}t/a$')
        axarr_sub[0].set_xlabel(r'Stokes number, $\mathcal{S}t/a$')

xline = 0.8
for ax, a0 in zip(axarr[:, 0].flatten(), alpha0):
    ax.set_axis_off()
    ax.set_xlim([0, 1])
    ax.set_ylim([0, 1])
    ax.axvline(x=xline, color='k')
    ax.text(xline/2, 0.5, r'${:.0f}^\circ$'.format(
        a0), ha='center', va='center')

axarr[0, 0].axhline(y=0.8, xmax=xline, color='k')
axarr[0, 0].text(xline/2, 0.83, r'$\alpha~[^\circ]$', ha='center', va='bottom')

for ax, l in zip(axarr[:, 1:].flatten(), 'abcdefgh'):
    trans = mtransforms.ScaledTranslation(
        5/72, -5/72, fig.dpi_scale_trans)
    label = '({})'.format(l)
    ax.text(0.0, 1.0, label, transform=ax.transAxes + trans, color='k',
            va='top', ha='left')

fig.savefig(
    '../paper/figures/{}.pdf'.format(sys.argv[0].split(os.sep)[-1].replace('.py', '')), dpi=600)
