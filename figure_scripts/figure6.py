import glob
import os
import sys

import matplotlib.pyplot as plt
import matplotlib.transforms as mtransforms
import numpy as np
import template as tp
from netCDF4 import Dataset

plt.rcParams['figure.constrained_layout.hspace'] = 0
plt.rcParams['figure.constrained_layout.h_pad'] = 0.005


# %% Load data
path_data = '../data/output_data'
list_runs = glob.glob(os.path.join(path_data, '*.nc'))
datasets = np.array([Dataset(run) for run in list_runs])

# %% create data vectors
alpha, phi, St, H0, L0, Fr = np.array([[d.variables['alpha'][:].data, d.variables['phi'][:].data,
                                        d.variables['St'][:].data, d.variables['H0'][:].data,
                                        d.variables['L0'][:].data, d.variables['Fr'][:].data,
                                        ] for d in datasets]).T

Ha = np.array([d.variables['H_a'][:].data if 'H_a' in d.variables.keys()
              else d.variables['H0'][:].data for d in datasets])
authors, particles = np.array([[d.author, d.particle_type,
                                ] for d in datasets]).T

# %% graphic specifications
dataset_idx = np.vectorize(tp.datasets.get)(authors)

markers = np.vectorize(tp.marker_style.get)(particles)
markers[dataset_idx == 'SedFoam'] = 's'

facecolors = np.vectorize(tp.color_datasets.get)(dataset_idx)
edgecolors = np.full_like(facecolors, 'k')
edgecolors[H0/Ha < 0.2] = 'tab:red'

zorders = np.vectorize(lambda dataset: tp.datset_zorder[dataset])(dataset_idx)
random_order = np.arange(zorders.size)
rng = np.random.default_rng(1994)
rng.shuffle(random_order)
plot_idxs = np.lexsort((random_order, zorders))

# %% masks for plot
alphas = [0, 45]
alpha_pad = 1.5

fig, axarr = plt.subplots(2, 1, constrained_layout=True,
                          figsize=tp.half_figure_size_inv, sharex=True)

for alpha0, ax in zip(alphas, axarr.flatten()):
    mask = ((alpha > alpha0 - alpha_pad) &
            (alpha < alpha0 + alpha_pad))[plot_idxs]
    if alpha0 == 0:
        axins = ax.inset_axes([0.43, 0.53, 0.55, 0.45])
        axins.set_ylim(0, 0.8)
        axins.set_xlim(right=0.8)
        axins.set_xlabel('$\phi$', labelpad=0)
        axins.set_ylabel(r'$\mathcal{F}r$', labelpad=2)
        tp.mscatter(phi[plot_idxs][mask], Fr[plot_idxs][mask], ax=axins, m=markers[plot_idxs][mask],
                    facecolors=facecolors[plot_idxs][mask], edgecolors=edgecolors[plot_idxs][mask], lw=0.5)
    #
    tp.mscatter(phi[plot_idxs][mask], Fr[plot_idxs][mask], ax=ax, m=markers[plot_idxs][mask],
                facecolors=facecolors[plot_idxs][mask], edgecolors=edgecolors[plot_idxs][mask], lw=0.5)
    #
    ax.set_ylabel(r'Froude number, $\mathcal{F}r$')
    ax.set_ylim([0, 1.6])

ax.set_xlabel(r'Volume fraction, $\phi$')
ax.set_xscale('log')

for ax, l in zip(axarr.flatten(), 'ab'):
    trans = mtransforms.ScaledTranslation(
        5/72, -5/72, fig.dpi_scale_trans)
    label = '({})'.format(l)
    ax.text(0.0, 1.0, label, transform=ax.transAxes + trans, color='k',
            va='top', ha='left')

fig.savefig(
    '../paper/figures/{}.pdf'.format(sys.argv[0].split(os.sep)[-1].replace('.py', '')), dpi=600)
