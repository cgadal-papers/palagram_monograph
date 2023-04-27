import glob
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import template as tp
from netCDF4 import Dataset

plt.rcParams['xtick.top'] = False


def ang2sin(ang):
    return np.sin(np.radians(ang))


def sin2ang(sin):
    return np.degrees(np.arcsin(sin))


def mscatter(x, y, ax=None, m=None, **kw):
    import matplotlib.markers as mmarkers
    if not ax:
        ax = plt.gca()
    sc = ax.scatter(x, y, **kw)
    if (m is not None) and (len(m) == len(x)):
        paths = []
        for marker in m:
            if isinstance(marker, mmarkers.MarkerStyle):
                marker_obj = marker
            else:
                marker_obj = mmarkers.MarkerStyle(marker)
            path = marker_obj.get_path().transformed(
                marker_obj.get_transform())
            paths.append(path)
        sc.set_paths(paths)
    return sc


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
author_zorder = ['Rastello', 'Cyril', 'Cyril/Marie', 'Jean', 'Julien']

# %% masks for plot
mask_Stokes = (St > 4e-2) & (St < 6e-2)
# mask_phi = (phi > 0.12) & (phi < 0.18)
# mask_Stokes = (St < 5.5e-2)
mask_phi = (phi < 0.40)

# %% graphic vector for plots
alphas = np.ones_like(Fr)
alphas[(~mask_phi) & (authors == 'Jean')] = 0.2
alphas[(~mask_Stokes) & (authors != 'Jean')] = 0.2

markers = np.full_like(Fr, 'o', dtype=str)
markers[particles == 'Hydrogels'] = '*'
markers[authors == 'Julien'] = 's'
markers[H0/Ha < 0.2] = 'd'


# %% figure
layout = [['legend', 'legend'], [('(a)'), '(b)']]
fig, axarr = plt.subplot_mosaic(layout, figsize=tp.large_figure_size, constrained_layout=True,
                                gridspec_kw={'height_ratios': [0.001, 1]})

for label in ['(a)', '(b)']:
    ax = axarr[label]
    for author in author_zorder:
        mask = (authors == author)
        #
        Y = Fr if label == '(a)' else Fr*np.sqrt(aspect)
        #
        mscatter(np.sin(np.radians(alpha[mask])), Y[mask], ax=ax,
                 c=tp.color_setups[author], alpha=alphas[mask], label=author, m=markers[mask])
    ax.set_xlabel(r'$\sin \alpha$')
    ax.set_xlim(left=-0.012)

    secax = ax.secondary_xaxis('top', functions=(sin2ang, ang2sin))
    secax.set_xlabel(r'angle, $\alpha$ [deg.]')

axarr['(a)'].set_ylabel(r'Froude number, $\mathcal{F}r$')
axarr['(b)'].set_ylabel(r'$\sqrt{a}\mathcal{F}r$')
axarr['(a)'].set_ylim(0, 1.35)
axarr['(b)'].set_ylim(0, 1.85)

axarr['legend'].axis('off')
handles, labels = axarr['(a)'].get_legend_handles_labels()
leg = axarr['legend'].legend(handles, labels, loc="upper center", ncol=5,
                             borderaxespad=0, title='Datasets')
for lh in leg.legend_handles:
    lh.set_alpha(1)

fig.savefig(
    '../paper/figures/{}.pdf'.format(sys.argv[0].split(os.sep)[-1].replace('.py', '')), dpi=600)
