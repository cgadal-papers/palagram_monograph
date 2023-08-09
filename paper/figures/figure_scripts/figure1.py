import glob
import os
import sys

import matplotlib.pyplot as plt
import matplotlib.transforms as mtransforms
import numpy as np
import template as tp
from matplotlib.colors import to_rgba
from matplotlib.path import Path

plt.rcParams['figure.constrained_layout.hspace'] = 0
plt.rcParams['figure.constrained_layout.wspace'] = 0
plt.rcParams['figure.constrained_layout.h_pad'] = 0.005
plt.rcParams['figure.constrained_layout.w_pad'] = 0.005


def sind(x):
    return np.sin(x*np.pi/180.)


def cosd(x):
    return np.cos(x*np.pi/180.)


def Rotation_matrix(theta):
    return np.array([[cosd(theta), -sind(theta)], [sind(theta), cosd(theta)]])


# %%
# ## Sketches parameters
# colors
# color_water = 'aliceblue'
color_water = 'lightcyan'
# color_water_salt = 'lightgreen'
color_water_salt = '#F3F7D4'
color_sed = 'peru'
color_walls = 'k'
alpha_water = 1
color_mixing = 'grey'

# dimension parameters
tank_height = 45  # cm
tank_length = 165  # cm
door_pos = 15  # cm
door_height = 1.1*tank_height
door_pad = 0*tank_height
y_bottom = 0  # cm, at the end of the canal
x_bottom = 0  # cm, at the end of the canal
mixing_height = tank_height
mixing_pad = 0.2*tank_height
mixing_width = 0.15*mixing_height


SLOPES = np.array([-7, -40, 0])
water_height_facts = [0.89, 1, 0.89]  # proportion of tank height
width_ratios = tank_length*cosd(SLOPES) + np.abs(tank_height*sind(SLOPES))


figsize = (tp.large_figure_width, 0.84*tp.large_figure_width)
fig = plt.figure(figsize=figsize, layout='compressed')

subfigs = fig.subfigures(2, 1, height_ratios=(0.364, 0.636), hspace=0.005)

# %%%%%%% Sketches

axs_sketches = subfigs[0].subplots(
    1, 3, width_ratios=width_ratios, gridspec_kw={'wspace': 0.02})

ax_labels = ['a', 'b', 'c']
for i, (ax, ax_label, slope, water_height_fact) in enumerate(zip(axs_sketches.flatten(),
                                                                 ax_labels, SLOPES, water_height_facts)):
    ax.set_aspect('equal')
    ax.set_axis_off()
    # ax.set_xticks([])
    # ax.set_yticks([])

    # parameters
    water_height = water_height_fact*tank_height  # cm, at end of the canal

    # ## reference points
    slope_vec = np.array([cosd(slope), sind(slope)])
    slope_vec_up = np.array([-sind(slope), cosd(slope)])
    down_vec = np.array([0, -1])
    # tank
    bottom_right = np.array([x_bottom, y_bottom])
    bottom_left = bottom_right - tank_length*slope_vec
    top_right = bottom_right + tank_height*slope_vec_up
    top_left = bottom_left + tank_height*slope_vec_up
    # door
    bottom_door = bottom_left + door_pos*slope_vec + door_pad*slope_vec_up
    top_door = bottom_door + door_height*slope_vec_up
    # water (ambient)
    bottom_door_water = bottom_door
    bottom_right_water = bottom_right
    top_right_water = bottom_right_water + water_height*slope_vec_up
    if ax_label == 'b':
        top_door_water = top_right_water - \
            (tank_length - door_pos)*slope_vec
    else:
        top_door_water = top_right_water - \
            (tank_length-door_pos)*np.array([1, 0])/cosd(slope)
    #
    # water (reservoir)
    bottom_right_water_reservoir = bottom_door_water
    top_right_water_reservoir = top_door_water
    bottom_left_water_resevoir = bottom_left
    top_left_water_reservoir = top_door_water - \
        door_pos*np.array([1, 0])/cosd(slope)

    xy_water = [bottom_door_water, top_door_water,
                top_right_water, bottom_right_water]
    xy_water_reservoir = [bottom_left_water_resevoir, top_left_water_reservoir,
                          top_right_water_reservoir, bottom_right_water_reservoir]

    # mixing
    bottom_mixing = bottom_left + 0.5*door_pos*slope_vec + mixing_pad*slope_vec_up
    top_mixing = bottom_mixing + mixing_height*slope_vec_up
    bottom_left_mixing = bottom_mixing - 0.5*mixing_width*slope_vec
    bottom_right_mixing = bottom_mixing + 0.5*mixing_width*slope_vec

    # sediment position generation
    ngrains = 70
    # np.random.seed(220212021)
    np.random.seed(999999)
    xsed = door_pos*np.random.random((ngrains, ))
    ysed = 1.2*water_height*np.random.random((ngrains, ))
    # if i == 0:
    #     xsed, ysed = door_pos * \
    #         np.random.random((ngrains, )), water_height * \
    #         (1+sind(slope))*np.random.random((ngrains, ))
    # else:
    #     xsed, ysed = door_pos*np.random.random((ngrains, )), 0.5*water_height*(
    #         1+sind(slope))*np.random.random((ngrains, ))
    xsed, ysed = (np.dot(Rotation_matrix(slope), np.array(
        [xsed, ysed])).T - tank_length*slope_vec).T

    # ## ploting sketch
    hpad = 0.007*tank_length
    ax.set_xlim(bottom_left[0] - 0.2*door_height*np.abs(sind(slope)) - hpad,
                top_right[0] + hpad)
    ax.set_ylim(bottom_right[1] - 0.2*door_height,
                top_door[1]+0.4*door_height)
    # ax.set_xlim(bottom_left[0] - 0.05*tank_length,
    #             top_right[0] + 0.05*tank_length)
    # ax.set_ylim(bottom_right[1] - 0.2*door_height, top_door[1]+0.4*door_height)

    #
    # ## tank walls
    ax.plot([bottom_left[0], bottom_right[0]], [
            bottom_left[1], bottom_right[1]], color=color_walls)
    ax.plot([bottom_left[0], top_left[0]], [
            bottom_left[1], top_left[1]], color=color_walls)
    ax.plot([bottom_right[0], top_right[0]], [
            bottom_right[1], top_right[1]], color=color_walls)
    ax.plot([bottom_door[0], top_door[0]], [
            bottom_door[1], top_door[1]], color=color_walls)
    if i == 1:
        ax.plot([top_door_water[0], top_right[0]], [
                top_door_water[1], top_right[1]], color=color_walls)
    # ## water
    poly_water = plt.Polygon(
        xy_water, facecolor=color_water, alpha=alpha_water, edgecolor=None)
    ax.add_patch(poly_water)
    poly_water_reservoir = plt.Polygon(xy_water_reservoir,
                                       facecolor=color_water_salt if ax_label == 'c' else color_water,
                                       alpha=alpha_water, edgecolor=None)
    ax.add_patch(poly_water_reservoir)

    # ## sediments
    mask_in = Path(xy_water_reservoir).contains_points(
        np.array([xsed, ysed]).T, radius=3)
    ax.scatter(xsed[mask_in], ysed[mask_in], color=color_sed, s=0.7)
    # ax.scatter(xsed[ysed < water_height], ysed[ysed <
    #            water_height], color=color_sed, s=0.7)

    # # ## mixing
    # ax.plot([bottom_mixing[0], top_mixing[0]], [
    #         bottom_mixing[1], top_mixing[1]], color=color_mixing)
    # ax.plot([bottom_left_mixing[0], bottom_right_mixing[0]],
    #         [bottom_left_mixing[1], bottom_right_mixing[1]], color=color_mixing, lw=2)

    # ## annotations
    ax.plot([bottom_right[0], bottom_right[0]-0.4*tank_length],
            [bottom_right[1], bottom_right[1]], ls='--', color='k')
    theta = np.linspace(180, 180+slope, 100)
    x, y = 0.35*tank_length * \
        np.array([cosd(theta), sind(theta)]) + bottom_right[:, None]
    ax.plot(x, y, color='k')
    #
    if ax_label in ['a', 'b']:
        ax.text(bottom_right[0]-0.4*tank_length, 0.75*y.mean(),
                r'$\alpha$', ha='right', va='center')

    ax.annotate("", xytext=top_door, xy=top_door+0.4*door_height*slope_vec_up,
                arrowprops=dict(arrowstyle="-|>", shrinkA=4, shrinkB=0, color='k'))

    xy = bottom_right - 0.195*door_height*slope_vec_up
    xytext = bottom_door - 0.195*door_height*slope_vec_up
    ax.annotate("", xytext=xytext, xy=xy, arrowprops=dict(arrowstyle="<->",
                                                          shrinkA=0, shrinkB=0,
                                                          color='k'))
    xytext = (xy + xytext)/2 - slope_vec_up*0.05*door_height
    ax.text(xytext[0], xytext[1],
            r'$L_{1}$', ha='center', va='top')
    #
    xy = bottom_left - 0.195*door_height*slope_vec_up
    xytext = bottom_door - 0.195*door_height*slope_vec_up
    ax.annotate("", xytext=xytext, xy=xy, arrowprops=dict(arrowstyle="<->",
                                                          shrinkA=0, shrinkB=0,
                                                          color='k'))
    xytext = (xy + xytext)/2
    ax.text(xytext[0], xytext[1] - 0.06*door_height,
            r'$L_{0}$', ha='center', va='top')

# %%%%%%% experimental images
subfigures_images = subfigs[1].subfigures(
    2, 1, height_ratios=(0.77, 0.20), hspace=0.022)

subfigures_images_top = subfigures_images[0].subfigures(
    1, 3, width_ratios=width_ratios)

# ## IMFT (sand80m_H19/run03)
path_imgs = 'src/images_figure1/IMFT'
images = [plt.imread(img) for img in sorted(
    glob.glob(os.path.join(path_imgs, '*.tif')))]

axs_IMFT = subfigures_images_top[0].subplots(
    len(images), 1, gridspec_kw={'hspace': 0.03})

for ax, img in zip(axs_IMFT, images):
    ax.imshow(img[100:, :150:-1], cmap=tp.cmap_images, vmax=2e4)
    ax.set_xticks([])
    ax.set_yticks([])

# ## LEGI (PMMA, à 45° (run7 dans les fichiers netcdf))
path_imgs = 'src/images_figure1/LEGI'
images = [plt.imread(img) for img in sorted(
    glob.glob(os.path.join(path_imgs, '*.tif')))]

img_sizes = [225, 490, 950]
axs_LEGI = subfigures_images_top[1].subplots(
    len(images), 1, height_ratios=img_sizes, gridspec_kw={'hspace': 0.02})

for ax, img, size in zip(axs_LEGI, images, img_sizes):
    ax.imshow(img[:size, :], cmap=tp.cmap_images)
    ax.set_xticks([])
    ax.set_yticks([])

# ## LEMTA (? manips à 55% )
path_imgs = 'src/images_figure1/LEMTA'
images = [plt.imread(img) for img in sorted(
    glob.glob(os.path.join(path_imgs, '*.tiff')))]

axs_LEMTA = subfigures_images_top[2].subplots(
    len(images), 1, gridspec_kw={'hspace': 0.03})

for ax, img in zip(axs_LEMTA, images):
    ax.imshow(img, cmap=tp.cmap_images, vmax=5.5e4)
    ax.set_xticks([])
    ax.set_yticks([])

# ## SEDFoam
path_imgs = 'src/images_figure1/SEDFOAM'
data = np.load(os.path.join(
    path_imgs, 'run07_snpachots_extracts.npy'), allow_pickle=True).item()

X, Y = np.meshgrid(data['x'], data['y'])

axs_SEDFOAM = subfigures_images[1].subplots(
    2, 5, width_ratios=[0.0001, 1, 1, 1, 0.0001], gridspec_kw={'hspace': 0.01, 'wspace': 0.03})
for ax, arr_phi in zip(axs_SEDFOAM[:, 1:-1].flatten(), data['phi']):
    ax.pcolormesh(X.T, Y.T, arr_phi,
                  cmap=tp.cmap_images2, vmax=0.025, rasterized=True)
    ax.set_aspect('equal')
    ax.set_xticks([])
    ax.set_yticks([])

for ax in axs_SEDFOAM[:, 0].flatten():
    ax.set_axis_off()
for ax in axs_SEDFOAM[:, -1].flatten():
    ax.set_axis_off()

letters = 'abcdefghijklmnopqrstuvwxyz'
axs = np.concatenate(
    [axs_sketches.flatten(), axs_IMFT.flatten(), axs_LEGI.flatten(), axs_LEMTA.flatten(), axs_SEDFOAM[:, 1:-1].flatten()])

for l, ax in zip(letters[:3], axs[:3]):
    trans_base = mtransforms.blended_transform_factory(
        ax.transAxes, subfigs[0].transFigure)
    trans = mtransforms.ScaledTranslation(
        1/72, -5/72, fig.dpi_scale_trans)
    label = '({})'.format(l)
    subfigs[0].text(0, 1, label, color='k', transform=trans_base + trans,
                    va='top', ha='left')

for l, ax in zip(letters[3:], axs[3:]):
    trans = mtransforms.ScaledTranslation(
        -3/72, -3/72, fig.dpi_scale_trans)
    label = '({})'.format(l)
    ax.text(1.0, 1.0, label, transform=ax.transAxes + trans, color='k',
            va='top', ha='right', bbox=dict(facecolor='w', edgecolor='none', alpha=0.8, pad=2.2))


# plt.show()
fig.savefig(
    '../{}.pdf'.format(sys.argv[0].split(os.sep)[-1].replace('.py', '')), dpi=400)
