from matplotlib.lines import Line2D
from matplotlib.legend_handler import HandlerTuple
import os

import cmasher as cmr
import cmocean as cmo
import matplotlib.pyplot as plt
import matplotlib.style as style
import numpy as np
from matplotlib import rcParams
from scipy.constants import golden

# %% Default parameters (matplotlib.rcParams)
style.use('tableau-colorblind10')

# %% Loading custom style
plt.style.use(os.path.join(os.path.dirname(__file__), 'src/style.mplstyle'))

# %%
# Defining some constants
inches_per_cm = 0.3937
regular_aspect_ratio = 1/golden


# %% Figure sizing

full_text_width = 492.982  # page width, in pt
column_width = 240.089  # one column width, in pt

large_figure_width = full_text_width*0.35136*0.1*inches_per_cm  # in inches
half_figure_width = column_width*0.35136*0.1*inches_per_cm  # in inches

half_figure_size = np.array([1, regular_aspect_ratio])*half_figure_width
large_figure_size = np.array([1, regular_aspect_ratio])*large_figure_width

half_figure_size_inv = np.array([1, 1/regular_aspect_ratio])*half_figure_width
large_figure_size_inv = np.array(
    [1, 1/regular_aspect_ratio])*large_figure_width

rcParams['figure.figsize'] = large_figure_size

# %% colors

# cmap_images = cmo.cm.ice
cmap_images = cmo.cm.gray
# cmap_images = cmo.cm.dense_r
# cmap_images = cmr.arctic

color_cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']

# color_setups = {'Julien': '#823329', 'Jean': '#FE7F2D',
#                 'Cyril': '#FCCA46', 'Cyril/Marie': '#619B8A', 'Rastello': '#A1C181'}

datasets = {'Julien': 'SedFoam', 'Jean': '3',
            'Cyril': '1', 'Cyril/Marie': '2', 'Rastello': '2'}

color_datasets = {'SedFoam': '#823329', '3': '#FE7F2D',
                  '1': '#FCCA46', '2': '#619B8A'}

marker_style = {'glass beads': 'o', 'silica sand': 'h', 'Hydrogels': '*',
                'PMMA': 'D', 'polystyren beads': 'X', 'SedFoam': 's'}

datset_zorder = {'1': 0, '2': 1, '3': 2, 'SedFoam': 2}

# %% corresponding legend

# leg_labels = sorted(color_setups.keys())
# legend_elements = [
#     Line2D([0], [0], marker='s' if author == 'Julien' else 'o', color=color_setups[author], ls='none') for author in sorted(color_setups.keys())
# ]

# legend_elements[leg_labels.index('Rastello')] = (
#     Line2D([0], [0], marker='o', color=color_setups['Rastello'], ls='none'),
#     Line2D([0], [0], marker='d', color=color_setups['Rastello'], ls='none'),
#     Line2D([0], [0], marker='*', color=color_setups['Rastello'], ls='none'),
# )

legend_datasets = [
    Line2D([0], [0], marker='s' if dataset == 'SedFoam' else 'o', color=color_datasets[dataset], ls='none', label=dataset) for dataset in sorted(color_datasets.keys())
]

legend_particles = [
    Line2D([0], [0], marker=marker_style[particle], markerfacecolor='none', markeredgecolor='k', ls='none', label=particle) for particle in sorted(marker_style.keys())
]

# %% plot functions


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
