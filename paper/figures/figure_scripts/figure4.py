import glob
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import template as tp
from matplotlib.colors import to_rgba
from models import Birman, Froude
from netCDF4 import Dataset

plt.rcParams["xtick.top"] = False


def ang2sin(ang):
    return np.sin(np.radians(ang))


def sin2ang(sin):
    return np.degrees(np.arcsin(sin))


# %% Load data
path_data = "../../../data/output_data"
list_runs = glob.glob(os.path.join(path_data, "*.nc"))
datasets = np.array([Dataset(run) for run in list_runs])

# %% create data vectors
alpha, phi, Fr, St, H0 = np.array(
    [
        [
            d.variables["alpha"][:].data,
            d.variables["phi"][:].data,
            d.variables["Fr"][:].data,
            d.variables["St"][:].data,
            d.variables["H0"][:].data,
        ]
        for d in datasets
    ]
).T

Ha = np.array(
    [d.variables["H_a"][:].data if "H_a" in d.variables.keys() else d.variables["H0"][:].data for d in datasets]
)
authors, particles = np.array(
    [
        [
            d.author,
            d.particle_type,
        ]
        for d in datasets
    ]
).T

# %% graphic specifications
# %% masks for plot
mask_phi = phi < tp.phi_c

# %% graphic vector for plots
alphas = np.ones_like(Fr)
alphas[~mask_phi] = 0.4

dataset_idx = np.vectorize(tp.datasets.get)(authors)

markers = np.vectorize(tp.marker_style.get)(particles)
markers[dataset_idx == "4"] = "s"

facecolors = np.vectorize(tp.color_datasets.get)(dataset_idx)
facecolors = np.array([to_rgba(c, a) for c, a in zip(facecolors, alphas)])

mask_nosuspended = (authors == "Rastello") & (H0 / Ha < 1)
edgecolors = np.array([to_rgba("k", a) for a in alphas])
edgecolors[mask_nosuspended] = np.array([to_rgba("tab:red", 0.4) for a in alphas[mask_nosuspended]])

zorders = np.vectorize(lambda dataset: tp.dataset_zorder2[dataset])(dataset_idx)
random_order = np.arange(zorders.size)
rng = np.random.default_rng(1994)
rng.shuffle(random_order)
plot_idxs = np.lexsort((random_order, zorders))

# %% theory
phi_m = 0.585
alpha_plot = np.linspace(0, 50, 500)
eta = 0
Fr_th = Froude(alpha_plot, eta, 7e4)

# %% figure
fig, ax = plt.subplots(1, 1, figsize=tp.large_figure_size, constrained_layout=True)


tp.mscatter(
    ang2sin(alpha)[plot_idxs],
    Fr[plot_idxs],
    ax=ax,
    m=markers[plot_idxs],
    facecolors=facecolors[plot_idxs],
    edgecolors=edgecolors[plot_idxs],
    lw=0.5,
)

(line_birman,) = ax.plot(
    ang2sin(alpha_plot),
    Birman(alpha_plot),
    ls="--",
    color="k",
    zorder=-8,
    lw=1,
    label="simulations of\nBirman et al. 2007",
)

(line_theory,) = ax.plot(ang2sin(alpha_plot), Fr_th, ls="-", color="k", zorder=-8, lw=1, label="eq. (15)")

ax.set_xlabel(r"$\sin \alpha$")
ax.set_xlim(left=-0.012)

secax = ax.secondary_xaxis("top", functions=(sin2ang, ang2sin))
secax.set_xlabel(r"angle, $\alpha$ [deg.]")

ax.set_ylabel(r"Froude number, $\mathcal{F}r$")
ax.set_ylim(0, 1.59)
ax.set_xlim(right=ang2sin(48))

other_elements = [line_birman, line_theory]
leg1 = ax.legend(handles=tp.legend_datasets + other_elements, ncol=2, title="Datasets")
leg2 = ax.legend(
    handles=tp.legend_particles, ncol=2, title="Particles", bbox_to_anchor=(0.365, 0.725), loc="lower left"
)

ax.add_artist(leg1)

fig.align_labels()

fig.savefig("../{}.pdf".format(sys.argv[0].split(os.sep)[-1].replace(".py", "")), dpi=600)
