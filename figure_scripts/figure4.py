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
alpha, phi, Fr, St = np.array([[d.variables['alpha'][:].data, d.variables['phi'][:].data,
                                d.variables['Fr'][:].data, d.variables['St'][:].data]
                               for d in datasets]).T
authors = np.array([d.author for d in datasets])

# %% graphic specifications
# changing order depending on authors
author_zorder = ['Rastello', 'Cyril', 'Cyril/Marie', 'Julien', 'Jean']

# %% masks for plot
mask_Stokes = (St > 4e-2) & (St < 6e-2)
# mask_phi = (phi > 0.12) & (phi < 0.18)
# mask_Stokes = (St < 5.5e-2)
mask_phi = (phi < 0.20) & (phi < 0.40)

# %% graphic vector for plots
alphas = np.ones_like(Fr)
alphas[(~mask_phi) & (authors == 'Jean')] = 0.2
alphas[(~mask_Stokes) & (authors != 'Jean')] = 0.2

# %% figure
fig, ax = plt.subplots(1, 1, constrained_layout=True,
                       figsize=tp.large_figure_size)
for author in author_zorder:
    mask = (authors == author)
    marker = 's' if author == 'Julien' else None
    #
    ax.scatter(np.sin(np.radians(alpha[mask])), Fr[mask],
               c=tp.color_setups[author], alpha=alphas[mask], label=author, marker=marker)
ax.set_xlabel(r'$\sin \alpha$')
ax.set_ylabel(r'Froude number, $\mathcal{F}r$')
ax.set_xlim(left=-0.012)
ax.set_ylim(0, 1.35)

leg = ax.legend(ncols=3)
for lh in leg.legend_handles:
    lh.set_alpha(1)

plt.show()
