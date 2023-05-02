import glob
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import Dataset


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
path_data = '../../data/output_data'
list_runs = glob.glob(os.path.join(path_data, '*.nc'))
datasets = np.array([Dataset(run) for run in list_runs])

# %% create data vectors
alpha, phi, Fr, St, L0, H0 = np.array([[d.variables['alpha'][:].data, d.variables['phi'][:].data,
                                        d.variables['Fr'][:].data, d.variables['St'][:].data,
                                        d.variables['L0'][:].data, d.variables['H0'][:].data,
                                        ] for d in datasets]).T
authors = np.array([d.author for d in datasets])
particles = np.array([d.particle_type for d in datasets])
Ha = np.array([d.variables['H_a'][:].data if 'H_a' in d.variables.keys()
              else d.variables['H0'][:].data for d in datasets])

a = H0/L0
r = H0/Ha

# %% graphic specifications
# changing order depending on authors
author_zorder = ['Rastello', 'Cyril', 'Cyril/Marie', 'Julien', 'Jean']

# %% masks for plot
# mask_alpha = (alpha < 0.5)
mask_alpha = (alpha > 5) & (alpha < 11)
mask_Stokes = (St > 4e-2) & (St < 6e-2)

# %% graphic vector for plots
alphas = np.ones_like(Fr)
# alphas[(~mask_phi) & (authors == 'Jean')] = 0.2
# alphas[((~mask_alpha) | (~mask_Stokes)) & (authors != 'Jean')] = 0.2
alphas[(~mask_Stokes) & (authors != 'Jean')] = 0.2

markers = np.full_like(Fr, 'o', dtype=str)
markers[particles == 'Hydrogels'] = '*'
markers[authors == 'Julien'] = 's'
markers[H0/Ha < 0.2] = 'd'


# %% figure
fig, ax = plt.subplots(1, 1, constrained_layout=True)
for author in author_zorder:
    mask = (authors == author) & mask_alpha
    # mask = (authors == author)
    #
    if mask.any():
        mscatter(a[mask], Fr[mask], ax=ax,
                 label=author, m=markers[mask], alpha=alphas[mask])
        print(Fr[mask].mean())
ax.set_xlabel(r'$a$')
ax.set_ylabel(r'Froude number, $\mathcal{F}r$')
# ax.set_yscale('log')
# ax.set_xscale('log')
# ax.set_xlim(0, 5)
# ax.set_ylim(0, 1.35)
xtest = np.linspace(0, 5, 100)
ax.plot(xtest, 0.5/np.sqrt(xtest), color='k')

# secax = ax.secondary_xaxis('top', functions=(sin2ang, ang2sin))
# secax.set_xlabel('angle [deg.]')

ax.legend()
# for lh in leg.legend_handles:
#     lh.set_alpha(1)

# fig.savefig(
#     '../paper/figures/{}.pdf'.format(sys.argv[0].split(os.sep)[-1].replace('.py', '')), dpi=600)

plt.show()
