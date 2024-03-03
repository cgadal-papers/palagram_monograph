import glob
import os
import sys

import matplotlib.pyplot as plt
import matplotlib.transforms as mtransforms
import numpy as np
import template as tp
from lmfit.model import load_modelresult
from matplotlib.lines import Line2D
from netCDF4 import Dataset


def polyFULL(t, Fr, lamb, xi, c, d):
    return xi + Fr * (t - lamb * t**2 + c * t**3 + d * t**4)


# %% Load data
path_data = "../../../data/output_data"
list_runs = sorted(glob.glob(os.path.join(path_data, "*.nc")))
list_fitresults = sorted(glob.glob(os.path.join(path_data, "fitresult*")))

datasets = np.array([Dataset(run) for run in list_runs])

# %% mask data
particles = np.array([d.particle_type for d in datasets]).T
mask = particles != "saline water"

zorder_setups = {
    "Cyril Gadal": -9,
    "Marie Rastello": -10,
    "Jean Schneider": -6,
    "Julien Chauchat": -7,
    "Marie Rastello - Cyril Gadal": -10,
}

# %% figure
layout = [["legend", "legend"], ["(a)", "(b)"]]
fig, axarr = plt.subplot_mosaic(
    layout,
    figsize=(tp.large_figure_width, 0.55 * tp.large_figure_width),
    constrained_layout=True,
    gridspec_kw={"height_ratios": [0.001, 1]},
)

# All runs
ax = axarr["(a)"]
for d in datasets[mask]:
    t_ad = d.variables["t0"][:].data
    x_ad = d.variables["L0"][:].data
    #
    x_axis = d.variables["t"][:].data / t_ad
    y_axis = d.variables["x_front"][:].data / x_ad
    if d.author == "Marie Rastello":
        y_axis[x_axis < 0.1] = np.nan
    ax.scatter(
        x_axis,
        y_axis,
        color=tp.color_datasets_BW[tp.datasets[d.author]],
        s=3,
        marker=tp.marker_style[d.particle_type] if d.author != "Julien" else "s",
        rasterized=True,
        zorder=zorder_setups[d.author],
    )

ax.set_ylabel(r"Front position, $x_{\rm f}/L_{0}$")
ax.set_xlabel(r"Time, $t/t_{0}$")
ax.set_xlim(-1, 150)
# ax.set_ylim(-0.5, 17)

# selected run, non-dimensional
ax = axarr["(b)"]
i_runs = [0, 158, 94, 197, 210, 202]

for i in i_runs:
    d = datasets[list_runs.index(os.path.join(path_data, f"run_{i:03d}.nc"))]
    # print(d.author)
    #
    t_ad = d.variables["t0"][:].data
    x_ad = d.variables["L0"][:].data
    #
    x_axis = d.variables["t"][:].data / t_ad
    y_axis = d.variables["x_front"][:].data / x_ad
    ax.scatter(
        x_axis,
        y_axis,
        color=tp.color_datasets_BW[tp.datasets[d.author]],
        s=3,
        marker=tp.marker_style[d.particle_type] if d.author != "Julien" else "s",
        rasterized=True,
        zorder=zorder_setups[d.author],
    )
    # # plotting fit
    fitresult = load_modelresult(
        os.path.join(path_data, f"fitresult_run_{i:03d}.save"), funcdefs={"polyFULL": polyFULL}
    )
    xplot = np.linspace(fitresult.userkws["t"].min(), fitresult.userkws["t"].max(), 500)
    ax.plot(xplot, fitresult.eval(t=xplot), color="k", lw=1, ls="--")

# annotations
ax.annotate(
    r"",
    xy=(19.3, 10.3),
    xycoords="data",
    xytext=(9, 12),
    textcoords="data",
    arrowprops=dict(arrowstyle="<-", shrinkA=0, shrinkB=0, connectionstyle="angle3,angleA=0,angleB=120"),
)
ax.text(15.5, 12, r"$\alpha$")

ax.annotate(
    r"",
    xy=(29.9, 13.45),
    xycoords="data",
    xytext=(36.6, 8.16),
    textcoords="data",
    arrowprops=dict(arrowstyle="<-", shrinkA=0, shrinkB=0, connectionstyle="angle3,angleA=-90,angleB=-20"),
)
ax.text(33, 10.4, r"$\mathcal{S}t$")


ax.set_ylabel(r"Front position, $x_{\rm f}/L_{0}$")
ax.set_xlabel(r"Time, $t/t_{0}$")
ax.set_xlim(-1, 40)
ax.set_ylim(-0.5, 17)

axarr["legend"].axis("off")
other_elements = [Line2D([0], [0], color="k", ls="--", label="fit")]
leg = axarr["legend"].legend(
    handles=tp.legend_datasets + other_elements, ncol=5, title="Datasets", loc="upper center", borderaxespad=0
)

for label, ax in axarr.items():
    if label not in ["legend"]:
        trans = mtransforms.ScaledTranslation(5 / 72, -5 / 72, fig.dpi_scale_trans)
        ax.text(0.0, 1.0, label, transform=ax.transAxes + trans, color="k", va="top", ha="left")

fig.align_labels()

fig.savefig("../{}.pdf".format(sys.argv[0].split(os.sep)[-1].replace(".py", "")), dpi=600)
