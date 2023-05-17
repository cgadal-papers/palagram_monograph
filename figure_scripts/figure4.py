import glob
import os
import sys

import matplotlib.pyplot as plt
import matplotlib.transforms as mtransforms
from matplotlib.lines import Line2D
import numpy as np
import template as tp
from netCDF4 import Dataset
from matplotlib.colors import to_rgba


plt.rcParams['xtick.top'] = False


def ang2sin(ang):
    return np.sin(np.radians(ang))


def sin2ang(sin):
    return np.degrees(np.arcsin(sin))


def Birman(ang):
    return -0.1824*np.radians(ang)**2 + 0.2781*np.radians(ang) + 0.4871


# %% Load data
path_data = '../data/output_data'
list_runs = glob.glob(os.path.join(path_data, '*.nc'))
datasets = np.array([Dataset(run) for run in list_runs])

# %% create data vectors
alpha, phi, Fr, St = np.array([[d.variables['alpha'][:].data, d.variables['phi'][:].data,
                                d.variables['Fr'][:].data, d.variables['St'][:].data]
                               for d in datasets]).T
authors = np.array([d.author for d in datasets])
particles = np.array([d.particle_type for d in datasets])

H0 = np.array([d.variables['H0'][:].data for d in datasets])
L0 = np.array([d.variables['L0'][:].data for d in datasets])
Ha = np.array([d.variables['H_a'][:].data if 'H_a' in d.variables.keys()
              else d.variables['H0'][:].data for d in datasets])
aspect = H0/L0

# %% graphic specifications
# %% masks for plot
# mask_Stokes = (St > 4e-2) & (St < 6e-2)
# mask_Stokes = (St < 6e-2)
mask_phi = (phi < 0.45)

# %% graphic vector for plots
alphas = np.ones_like(Fr)
# alphas[(~mask_phi | ~mask_Stokes)] = 0.2
alphas[~mask_phi] = 0.2

dataset_idx = np.vectorize(tp.datasets.get)(authors)

markers = np.vectorize(tp.marker_style.get)(particles)
markers[dataset_idx == 'SedFoam'] = 's'

facecolors = np.vectorize(tp.color_datasets.get)(dataset_idx)
facecolors = np.array([to_rgba(c, a) for c, a in zip(facecolors, alphas)])

# edgecolors = np.copy(facecolors)
edgecolors = np.array([to_rgba('k', a) for a in alphas])
edgecolors[H0/Ha < 0.2] = np.array([to_rgba('tab:red', 1)
                                   for a in alphas[H0/Ha < 0.2]])

zorders = np.vectorize(lambda dataset: tp.datset_zorder[dataset])(dataset_idx)
random_order = np.arange(zorders.size)
rng = np.random.default_rng(1994)
rng.shuffle(random_order)
plot_idxs = np.lexsort((random_order, zorders))


# %% figure
fig, ax = plt.subplots(1, 1, figsize=tp.large_figure_size,
                       constrained_layout=True)

alpha_plot = np.linspace(0, 50, 500)

tp.mscatter(ang2sin(alpha)[plot_idxs], Fr[plot_idxs], ax=ax, m=markers[plot_idxs],
            facecolors=facecolors[plot_idxs], edgecolors=edgecolors[plot_idxs], lw=0.5)

ax.plot(ang2sin(alpha_plot), Birman(alpha_plot),
        ls='--', color='k', zorder=-8)

ax.set_xlabel(r'$\sin \alpha$')
ax.set_xlim(left=-0.012)

secax = ax.secondary_xaxis('top', functions=(sin2ang, ang2sin))
secax.set_xlabel(r'angle, $\alpha$ [deg.]')

ax.set_ylabel(r'Froude number, $\mathcal{F}r$')
ax.set_ylim(0, 1.6)
ax.set_xlim(right=ang2sin(48))

other_elements = [Line2D([0], [0], color='k', ls='--',
                         label='Birman et al. 2007')]
leg1 = ax.legend(handles=tp.legend_datasets +
                 other_elements, ncol=2, title='Datasets')
leg2 = ax.legend(handles=tp.legend_particles, ncol=2,
                 title='Particles', bbox_to_anchor=(0.365, 0.725), loc='lower left')
# leg2 = ax.legend(handles=tp.legend_particles, ncol=2,
#                  title='Particles', bbox_to_anchor=(0, 0.475), loc='lower left')
ax.add_artist(leg1)

fig.align_labels()

fig.savefig(
    '../paper/figures/{}.pdf'.format(sys.argv[0].split(os.sep)[-1].replace('.py', '')), dpi=600)
