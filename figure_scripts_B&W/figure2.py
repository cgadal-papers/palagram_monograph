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
datasets = np.array([Dataset(run) for run in list_runs if Dataset(
    run).particle_type != 'saline water'])

# %% create data vectors
alpha, phi, Re, St, H0 = np.array([[d.variables['alpha'][:].data, d.variables['phi'][:].data,
                                    d.variables['Re'][:].data, d.variables['St'][:].data,
                                    d.variables['H0'][:].data,
                                    ] for d in datasets]).T

authors, particles = np.array([[d.author, d.particle_type,
                                ] for d in datasets]).T

Ha = np.array([d.variables['H_a'][:].data if 'H_a' in d.variables.keys()
              else d.variables['H0'][:].data for d in datasets])

# %% graphic vector for plots
dataset_idx = np.vectorize(tp.datasets.get)(authors)

markers = np.vectorize(tp.marker_style.get)(particles)
markers[dataset_idx == '4'] = 's'

# facecolors = np.vectorize(tp.color_datasets_BW.get)(dataset_idx)
# edgecolors = np.full_like(facecolors, 'k')
facecolors = np.array([tp.color_datasets_BW[d] for d in dataset_idx])
edgecolors = np.array(['k' for d in dataset_idx])

mask_nosuspended = (authors == 'Rastello') & (H0/Ha < 1)
edgecolors[mask_nosuspended] = tp.to_grayscale('tab:red')[0]

zorders = np.vectorize(lambda dataset: tp.dataset_zorder[dataset])(dataset_idx)
random_order = np.arange(zorders.size)
rng = np.random.default_rng(1994)
rng.shuffle(random_order)
plot_idxs = np.lexsort((random_order, zorders))

# %% figure

figsize = (tp.half_figure_width, 2.4*tp.half_figure_width)
fig, axarr = plt.subplots(5, 1, figsize=figsize, constrained_layout=True,
                          gridspec_kw={'height_ratios': [0.2, 0.2, 1, 1, 1]})

for ax, var in zip(axarr[2:].flatten(), [phi, alpha, St]):
    tp.mscatter(Re[plot_idxs], var[plot_idxs], ax=ax, m=markers[plot_idxs],
                facecolors=facecolors[plot_idxs], edgecolors=edgecolors[plot_idxs], lw=0.5)
    ax.set_xscale('log')
    ax.set_yscale('log')

# legends
axarr[0].axis('off')
axarr[1].axis('off')

leg = axarr[0].legend(handles=tp.legend_datasets, loc="upper center",
                      ncol=4, borderaxespad=0, title='Datasets', mode='expand')
leg = axarr[1].legend(handles=tp.legend_particles, loc="upper center",
                      ncol=2, borderaxespad=0, title='Particles', mode='expand')
#
axarr[2].set_xticks([])
axarr[3].set_xticks([])
axarr[4].set_xlabel(r'Reynolds number, $\mathcal{R}e$')
#
axarr[2].set_ylabel(r'Volume fraction, $\phi$')
axarr[3].set_ylabel(r'angle, $\alpha~[^\circ]$')
axarr[4].set_ylabel(r'Stokes number, $\mathcal{S}t$')
#
axarr[3].set_yscale('linear')
axarr[3].set_ylim(bottom=-1)

for ax, l in zip(axarr[2:].flatten(), 'abcdefgh'):
    trans = mtransforms.ScaledTranslation(
        5/72, -5/72, fig.dpi_scale_trans)
    label = '({})'.format(l)
    ax.text(0.0, 1.0, label, transform=ax.transAxes + trans, color='k',
            va='top', ha='left')

fig.align_labels()

fig.savefig(
    '../paper/figures_B&W/{}.pdf'.format(sys.argv[0].split(os.sep)[-1].replace('.py', '')), dpi=600)
