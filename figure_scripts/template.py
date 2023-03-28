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
# color_setups = {'Cyril': color_cycle[0], 'Rastello': color_cycle[2],
#                 'Jean': color_cycle[1], 'Julien': color_cycle[5],
#                 'Cyril/Marie': color_cycle[4]}
color_setups = {'Cyril': color_cycle[0], 'Rastello': color_cycle[2],
                'Jean': color_cycle[1], 'Julien': 'tab:green',
                'Cyril/Marie': 'tab:purple'}
# color_setups = {'IMFT': color_cycle[0], 'LEGI': color_cycle[1],
#                 'LEMTA': color_cycle[2], 'NUM': color_cycle[3]}
