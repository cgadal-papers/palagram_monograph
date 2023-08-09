import glob
import os
import sys

import matplotlib.pyplot as plt
from PIL import Image
import numpy as np
import template as tp

img = np.array(Image.open(
    'src/image_fig7/run07_2D_3DCyclic_t=15s_V2.png'))[15:290, 60:]

figsize = (tp.large_figure_width, tp.large_figure_width *
           img.shape[0]/img.shape[1])

fig, ax = plt.subplots(1, 1, constrained_layout=True,
                       figsize=figsize, sharex=True)

ax.imshow(img)

ax.set_xticks([])
ax.set_yticks([])

fig.savefig(
    '../{}.pdf'.format(sys.argv[0].split(os.sep)[-1].replace('.py', '')), dpi=600)
