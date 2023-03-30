import glob
import os

import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import Dataset
import numpy as np
from scipy.optimize import curve_fit


# def logistiqueNUM(x, a, b, c, d, e):
#     y = a*(x-0) - b*(x-0)**(2)/2 + c + d*x**3 + e*x**4
#     return y

def logistiqueNUM(x, a, b, c, d, e):
    y = a*(x-0) - b*(x-0)**(2)/2 + c + d*x**3 + e*x**4
    return y


def logistique(x, a, b, c):
    return logistiqueNUM(x, a, c, b=0, d=0, e=0)


def compute_rsquared(ydata, yfit):
    residuals = ydata - yfit
    ss_res = np.sum(residuals**2)
    ss_tot = np.sum((ydata-np.nanmean(ydata))**2)
    r2 = 1 - (ss_res / ss_tot)
    return r2


def Stokes_Velocity(d, mu, rho_p, rho_f, g):
    return d**2*(rho_p - rho_f) * g/mu/18


# %% Variable definition
# paths
input_path = '../../data/input_data'

# physical parameters
g = 9.81  # [m/s2]
mu = 1e-3  # [kg/m/s]

# fit specifications
BOUNDS_FIT = {'Cyril': (3, 40), 'Cyril/Marie': (8, 100), 'Jean': (0.5, 20),
              'Julien': (3, 80), 'Rastello': (0, 60)}

FUNC_FIT = {'Cyril': logistiqueNUM, 'Cyril/Marie': logistiqueNUM, 'Jean': logistique,
            'Julien': logistiqueNUM, 'Rastello': logistiqueNUM}

# other
SETUPS = {'Cyril': 'IMFT', 'Cyril/Marie': 'LEGI', 'Jean': 'LEMTA',
          'Julien': 'NUM', 'Rastello': 'LEGI'}

# %% Loading data
# list_runs = np.array(glob.glob(os.path.join(input_path, 'runs_JULIEN*/*.nc')))
list_runs = np.array(glob.glob(os.path.join(input_path, 'runs_JEAN/*.nc')))
# list_runs = np.array(glob.glob(os.path.join(input_path, 'runs_MARIE/*.nc')))
# list_runs = np.array(glob.glob(os.path.join(input_path, 'runs_CYRIL/*.nc')))
datasets = np.array([Dataset(run) for run in list_runs])

alphas = np.array([np.degrees(d.variables['alpha'][:].data) if d.author ==
                  'Julien' else d.variables['alpha'][:].data for d in datasets])
grains = np.array([d.particle_type for d in datasets])
authors = np.array([d.author for d in datasets])
# mask_runs = alphas < 0.5
mask_runs = ((alphas > 6) & (alphas < 10)) & (
    grains == 'glass beads') & (authors == 'Cyril')
# mask_runs = ((alphas > 6) & (alphas < 10)) & (
#     grains == 'silica sand') & (authors == 'Cyril')
# mask_runs = (alphas > 30)
# mask_runs = (authors == 'Cyril')
mask_runs = np.ones_like(alphas).astype('bool')

# %% Loop over data file and analysis

fig_fit, ax = plt.subplots(1, 1, constrained_layout=True)

for i, d in enumerate(datasets[mask_runs]):
    run = list_runs[mask_runs][i].split(os.sep)[-1]
    print(run)
    #
    t = d.variables['t'][:].data
    x_front = d.variables['x_front'][:].data
    #
    # Loading some variables
    H0 = d.variables['H0'][:].data        # lock characteristic height, [m]
    L0 = d.variables['L0'][:].data        # lock characteristic width, [m]
    rho_p = d.variables['rho_p'][:].data  # particle velocity, [kg/m3]
    rho_f = d.variables['rho_f'][:].data  # lock fluid density, [kg/m3]
    rho_a = d.variables['rho_a'][:].data  # ambiant fluid density, [kg/m3]
    alpha = d.variables['alpha'][:].data  # bottom slope, [deg.]
    diam = d.variables['d'][:].data  # grain size, [m]
    phi = d.variables['phi'][:].data
    if d.author == 'Julien':
        H0 = H0/100
        L0 = L0/100
        rho_f, rho_p, rho_a = rho_f*1e3, rho_p*1e3, rho_a*1e3
        alpha = alpha*180/np.pi
        diam = diam*1e-6  # grain size, [m]
    #
    # Computing variables for adi time
    rho_c = rho_f + phi * (rho_p - rho_f)  # average lock density, [kg/m3]
    gprime = g*(rho_c - rho_a)/rho_a  # specific gravity
    u0 = np.sqrt(gprime*H0)
    t_ad = L0/u0
    #
    #
    # #### Fitting front position curves
    # defining fitting masks
    mask_ok = ~np.isnan(x_front)
    t_ok, x_ok = t[mask_ok], x_front[mask_ok]
    bounds_fit = BOUNDS_FIT[d.author]
    if (d.author == 'Cyril'):
        if run == 'run_012.nc':
            bounds_fit = (bounds_fit[0], 30)
        elif run == 'run_006.nc':
            bounds_fit = (bounds_fit[0], 42)
        elif run == 'run_020.nc':
            bounds_fit = (bounds_fit[0], 16)
        elif run == 'run_116.nc':
            bounds_fit = (bounds_fit[0], 20)
        elif run == 'run_114.nc':
            bounds_fit = (bounds_fit[0], 17)
    #
    mask = (t_ok > bounds_fit[0]*t_ad) & (t_ok < bounds_fit[1]*t_ad)
    #
    # Make fit
    func_fit = FUNC_FIT[d.author]
    if mask.any():
        p, pcov = curve_fit(func_fit, t_ok[mask], x_ok[mask])
        perr = np.sqrt(np.diag(pcov))
        r_squared = compute_rsquared(
            x_ok[mask], func_fit(t_ok[mask], *p))
        #
        print('author: {}, r2: {:.3f}'.format(d.author, r_squared))
        # if (d.author == 'Rastello') & ~((perr[0] < 1e-5) or (t[-1]/t_ad > 30)):
        #     mask = (t_ok < 13*t_ad) & (t_ok > 3*t_ad)
        #     p, pcov = curve_fit(func_fit, t_ok[mask], x_ok[mask])
        #     perr = np.sqrt(np.diag(pcov))
        #     r_squared = compute_rsquared(
        #         x_ok[mask], func_fit(t_ok[mask], *p), p.size)
        #     print('Bad Fit, new r2: {:.3f}'.format(r_squared))
        #
        #
        ax.plot(t_ok[mask]/t_ad, x_ok[mask]/L0,
                lw=5, alpha=0.5, color='tab:orange')
        ax.plot(t/t_ad, x_front/L0, '.-', color='tab:blue', lw=1)
        ax.plot(t_ok[mask]/t_ad, func_fit(t_ok[mask], *p) /
                L0, color='k', ls='--')
        if func_fit == logistique:
            ax.text(t_ok[mask][-1]/t_ad, x_ok[mask][-1]/L0,
                    '{}, {:.0f}, {:.0f}:{} \n {:.1e}, {:.1e}'.format(phi, alpha, diam*1e6, d.particle_type, p[0]/u0, p[1]/gprime), ha='left')
        else:
            ax.text(t_ok[mask][-1]/t_ad, x_ok[mask][-1]/L0,
                    '{}: {}, {:.0f}, {:.0f}:{} \n {:.0e}, {:.0e}, {:.0e}, {:.0e}'.format(run, phi, alpha, diam*1e6, d.particle_type, p[0]/u0, p[1]/gprime, p[3]/gprime/t_ad, p[4]/gprime/t_ad**2), ha='left')
    else:
        p = np.array([np.nan, np.nan, np.nan])
        perr = np.copy(p)
        r_squared = np.nan
        print('All nans')

    # ax.plot(t, x_front, color='tab:blue')
    # ax.plot(t_ok[mask], x_ok[mask], lw=3, alpha=0.5, color='tab:orange')
    # ax.plot(t_ok[mask], func_fit(t_ok[mask], *p), color='k', ls='--')
    # ax.text(t_ok[-1], x_ok[-1],
    #         '{}, {:.0f}, {:.0f}:{} \n {:.0e}, {:.0e}'.format(phi, alpha, diam*1e6, d.particle_type, p[0]/u0, p[1]/gprime), ha='left')

ax.set_xlim(left=0)
ax.set_ylim(bottom=0)
plt.show()
