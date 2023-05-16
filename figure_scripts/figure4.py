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
# changing order depending on authors
author_zorder = ['Cyril', 'Rastello', 'Cyril/Marie', 'Jean', 'Julien']

# %% masks for plot
# mask_Stokes = (St > 4e-2) & (St < 6e-2)
# mask_Stokes = (St < 6e-2)
mask_phi = (phi < 0.45)

# %% graphic vector for plots
alphas = np.ones_like(Fr)
# alphas[(~mask_phi | ~mask_Stokes)] = 0.2
alphas[~mask_phi] = 0.2

markers = np.full_like(Fr, 'o', dtype=str)
markers[particles == 'Hydrogels'] = '*'
markers[authors == 'Julien'] = 's'
markers[H0/Ha < 0.2] = 'd'


# %% figure
layout = [['legend', 'legend'], ['(a)', '(b)']]
fig, axarr = plt.subplot_mosaic(layout, figsize=tp.large_figure_size, constrained_layout=True,
                                gridspec_kw={'height_ratios': [0.001, 1]})

alpha_plot = np.linspace(0, 50, 500)
for label in ['(a)', '(b)']:
    ax = axarr[label]
    for author in author_zorder:
        mask = (authors == author)
        #
        Y = Fr if label == '(a)' else Fr*np.sqrt(aspect)
        #
        edgecolors = [to_rgba(tp.color_setups[author], a)
                      for a in alphas[mask]]
        # edgecolors = [(0.4, 0.4, 0.4, a)
        #               for a in alphas[mask]]
        facecolors = [to_rgba(tp.color_setups[author], a)
                      for a in alphas[mask]]
        tp.mscatter(ang2sin(alpha[mask]), Y[mask], ax=ax,
                    edgecolors=edgecolors, label=author, m=markers[mask], facecolors=facecolors)
    #
    fact = 1 if label == '(a)' else np.sqrt(0.1)
    ax.plot(ang2sin(alpha_plot), fact*Birman(alpha_plot),
            ls='--', color='k', zorder=-8)

    ax.set_xlabel(r'$\sin \alpha$')
    ax.set_xlim(left=-0.012)

    secax = ax.secondary_xaxis('top', functions=(sin2ang, ang2sin))
    secax.set_xlabel(r'angle, $\alpha$ [deg.]')

axarr['(a)'].set_ylabel(r'Froude number, $\mathcal{F}r$')
axarr['(b)'].set_ylabel(r'$\sqrt{a}\mathcal{F}r$')
axarr['(a)'].set_ylim(0, 1.6)
axarr['(b)'].set_ylim(0, 1.85)
axarr['(a)'].set_xlim(right=ang2sin(48))
axarr['(b)'].set_xlim(right=ang2sin(48))

axarr['legend'].axis('off')
other_elements = [Line2D([0], [0], color='k', ls='--',
                         label='Birman et al. 2007')]
leg = axarr['legend'].legend(handles=tp.legend_elements + other_elements,
                             ncol=5, title='Datasets', loc="upper center", borderaxespad=0)

for label, ax in axarr.items():
    if label not in ['legend']:
        trans = mtransforms.ScaledTranslation(
            5/72, -5/72, fig.dpi_scale_trans)
        ax.text(0.0, 1.0, label, transform=ax.transAxes + trans, color='k',
                va='top', ha='left')

fig.align_labels()

fig.savefig(
    '../paper/figures/{}.pdf'.format(sys.argv[0].split(os.sep)[-1].replace('.py', '')), dpi=600)
