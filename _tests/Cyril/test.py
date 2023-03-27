import glob
import os

import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import Dataset

data_path = '../data'

# list_runs = glob.glob(os.path.join(data_path, 'runs_MARIE/*.nc'))
list_runs = glob.glob(os.path.join(data_path, 'runs_CYRIL/*.nc'))


# list_runs = glob.glob(os.path.join(data_path, '*/*.nc'))

datasets = [Dataset(run) for run in sorted(list_runs)]
authors = [d.author for d in datasets]

RHO_F = []

fig, axarr = plt.subplots(
    1, 2, constrained_layout=True, sharex=True, sharey=True)
#
for i, d in enumerate(datasets):
    # print(d.author if hasattr(d, 'author') else 'Marie')
    t = d.variables['t'][:].data
    x_front = d.variables['x_front'][:].data
    #
    H = d.variables['H0'][:].data
    alpha = d.variables['alpha'][:].data
    rho_f = d.variables['rho_f'][:].data
    rho_a = d.variables['rho_a'][:].data
    #
    rho_c = rho_f + \
        d.variables['phi'][:].data * \
        (d.variables['rho_p'][:].data - rho_f)
    #
    gprime = 9.81*(rho_c - rho_a)/rho_a
    #
    slope_corr = np.cos(
        alpha*np.pi/180) + 0*np.sin(alpha*np.pi/180)
    u0 = np.sqrt(gprime*H*slope_corr)
    t_ad = H/u0
    #
    u0 = np.sqrt(gprime*H)
    t_ad2 = H/u0
    #
    # if (x_front > 100).any():
    # x_front = x_front/100
    # print(d.author)
    # if d.variables['rho_p'][:].data in [1070, 1003, 1005, 1020]:
    if (d.variables['rho_p'][:].data > 1000) & (d.variables['rho_p'][:].data < 1025):
        ls = '--'
    else:
        ls = '-'
    # axarr[0].plot(t/t_ad, x_front/d.variables['L0'][:].data, lw=1, ls=ls)
    axarr[0].plot(t, x_front, lw=1, ls=ls)
    axarr[1].plot(t/t_ad2, x_front/d.variables['L0'][:].data, lw=1, ls=ls)

    #
    RHO_F.append(rho_f)

plt.show()

PHI = [d.variables['phi'][:].data for d in datasets]


# RHO_P = [d.variables['rho_p'][:].data for d in datasets]
# ALPHA = [d.variables['alpha'][:].data for d in datasets]
# RHO_A = [d.variables['rho_a']
#          [:].data for d in datasets if 'rho_f' in d.variables.keys()]

# RHO_F = [d.variables['rho_f']
#          [:].data for d in datasets if 'rho_f' in d.variables.keys()]

plt.figure()
plt.plot(RHO_F, '.')
plt.show()


# for d in datasets:
#     if d.author == 'Cyril/Marie':
#         print('{} : {}'.format(d.run_oldID, np.isnan(
#             d.variables['x_front'][:].data).all()))
