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

# %% graphic specifications
# changing order depending on authors
author_zorder = ['Rastello', 'Cyril', 'Cyril/Marie', 'Julien', 'Jean']

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
fig, ax = plt.subplots(1, 1, constrained_layout=True,
                       figsize=tp.large_figure_size)
for author in author_zorder:
    mask = (authors == author)
    #
    mscatter(np.sin(np.radians(alpha[mask])), Fr[mask], ax=ax,
             c=tp.color_setups[author], alpha=alphas[mask], label=author, m=markers[mask])
ax.set_xlabel(r'$\sin \alpha$')
ax.set_ylabel(r'Froude number, $\mathcal{F}r$')
ax.set_xlim(left=-0.012)
ax.set_ylim(0, 1.35)

secax = ax.secondary_xaxis('top', functions=(sin2ang, ang2sin))
secax.set_xlabel('angle [deg.]')

leg = ax.legend(ncols=3)
for lh in leg.legend_handles:
    lh.set_alpha(1)

fig.savefig(
    '../paper/figures/{}.pdf'.format(sys.argv[0].split(os.sep)[-1].replace('.py', '')), dpi=600)
