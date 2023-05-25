import glob
import os

import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import Dataset
import numpy as np
from lmfit import Model


def polyFULL(t, Fr, L, xi, c, d):
    return xi + Fr*t - L*t**2/2 + c*t**3 + d*t**4


def Stokes_Velocity(d, mu, rho_p, rho_f, g):
    return d**2*(rho_p - rho_f) * g/mu/18


def determine_fit_props(author, alpha, run, params):
    params['xi'].vary = False
    params['Fr'].vary = True
    params['L'].vary = True
    params['c'].vary = False
    params['d'].vary = False
    t_bounds = [5, 25]
    #
    if author == 'Jean':
        t_bounds = [0, 40]
        params['L'].vary = False
        params['xi'].vary = True
    if author == 'Julien':
        params['xi'].vary = True
        if run == 'run42b_front.nc':
            t_bounds = [0, 8]
        if run == 'run37_front.nc':
            t_bounds = [0, 20]
        if run == 'run42_front.nc':
            t_bounds = [0, 20]
    if author == 'Cyril':
        if run == 'run_012.nc':
            t_bounds[-1] = min(t_bounds[-1], 30)
        elif run == 'run_006.nc':
            t_bounds[-1] = min(t_bounds[-1], 42)
        elif run == 'run_020.nc':
            t_bounds[-1] = min(t_bounds[-1], 16)
        elif run == 'run_116.nc':
            t_bounds[-1] = min(t_bounds[-1], 20)
        elif run == 'run_114.nc':
            t_bounds[-1] = min(t_bounds[-1], 17)
    return t_bounds


# %% Variable definition
# paths
input_path = '../../data/input_data'

# physical parameters
g = 9.81  # [m/s2]
mu = 1e-3  # [kg/m/s]

# %% fit objects definition

# %% fit objects definition
# model object creation
model = Model(polyFULL)
params = model.make_params()

# parameter properties (Non dim.)
p0 = {'Fr': 0.4, 'xi': 0, 'L': 0, 'c': 0, 'd': 0}

lower_bounds = {'Fr': 0, 'xi': -np.inf,
                'L': -0.05, 'c': -1e-3, 'd': -1e-4}
higher_bounds = {'Fr': 2, 'xi': np.inf, 'L': 0.1, 'c': 1e-3, 'd': 1e-4}

# set parameter bounds
for par in params.keys():
    params[par].set(value=p0[par], min=lower_bounds[par],
                    max=higher_bounds[par])

# other
SETUPS = {'Cyril': 'IMFT', 'Cyril/Marie': 'LEGI', 'Jean': 'LEMTA',
          'Julien': 'NUM', 'Rastello': 'LEGI'}

# %% Loading data
list_runs = np.array(glob.glob(os.path.join(input_path, 'runs_JULIEN2/*.nc')))
# list_runs = np.array(glob.glob(os.path.join(input_path, 'runs_JEAN/*.nc')))
# list_runs = np.array(glob.glob(os.path.join(input_path, 'runs_MARIE/*.nc')))
# list_runs = np.array(glob.glob(os.path.join(input_path, 'runs_CYRIL/*.nc')))
datasets = np.array([Dataset(run) for run in list_runs])

alphas = np.array([np.degrees(d.variables['alpha'][:].data) if d.author ==
                  'Julien' else d.variables['alpha'][:].data for d in datasets])
grains = np.array([d.particle_type for d in datasets])
authors = np.array([d.author for d in datasets])
diams = np.array([d.variables['d'][:].data for d in datasets])
phi = np.array([d.variables['phi'][:].data for d in datasets])
H0 = np.array([d.variables['H0'][:].data for d in datasets])
# Ha = np.array([d.variables['H_a'][:].data for d in datasets])


authors = np.array([d.author for d in datasets])
# mask_runs = alphas < 0.5
# mask_runs = ((alphas > 6) & (alphas < 10)) & (
#     grains == 'glass beads') & (authors == 'Cyril') & (np.abs(diams - 6.4*1e-5) < 1e-5) & (phi > 3/100)
# mask_runs = (grains == 'saline water') & (authors == 'Cyril')
# mask_runs = ((alphas > 6) & (alphas < 10)) & (
#     grains == 'silica sand') & (authors == 'Cyril')
# mask_runs = (authors == 'Cyril/Marie')
# mask_runs = np.ones_like(alphas).astype('bool')
# mask_runs = (alphas > 6) & (alphas < 8)
mask_runs = (alphas > 40)
# mask_runs = (grains == 'PMMA')
# mask_runs = (H0/Ha < 0.3)
# %% Loop over data file and analysis

fig_fit, ax = plt.subplots(1, 1, constrained_layout=True)

for i, d in enumerate(datasets[mask_runs]):
    run = list_runs[mask_runs][i].split(os.sep)[-1]
    # run = list_runs[mask_runs][i].split('input_data')[-1]
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
        alpha = alpha*180/np.pi
    if rho_f < 500:
        H0 = H0/100
        L0 = L0/100
        rho_f, rho_p, rho_a = rho_f*1e3, rho_p*1e3, rho_a*1e3
        diam = diam*1e-6  # grain size, [m]
    #
    # Computing variables for adi time
    rho_c = rho_f + phi * (rho_p - rho_f)  # average lock density, [kg/m3]
    gprime = g*(rho_c - rho_a)/rho_a  # specific gravity
    # u0 = np.sqrt(gprime*H0*np.cos(np.radians(alpha)))
    u0 = np.sqrt(gprime*H0)
    t_ad = L0/u0
    # if diam > 400e-6:
    #     print(run)
    #
    #
    # #### Fitting front position curves
    t_bounds = determine_fit_props(d.author, alpha, run, params)
    # defining fitting masks
    mask_ok = ~np.isnan(x_front)
    t_ok, x_ok = t[mask_ok]/t_ad, x_front[mask_ok]/L0
    mask = (t_ok > t_bounds[0]) & (t_ok < t_bounds[1])
    #
    # Make fit
    ax.plot(t/t_ad, x_front/L0, '.-', color='tab:blue', lw=1)
    method = 'trust-constr'
    if mask.sum() > 3:
        result = model.fit(x_ok[mask], params, t=t_ok[mask], method=method)
        L, L_err = result.params['L'].value, result.params['L'].stderr
        if L > 2e-2:
            print('Quartic')
            params['d'].vary = True
            result = model.fit(x_ok[mask], params, t=t_ok[mask], method=method)
        r_squared = result.rsquared
        print('author: {}, r2: {:.3f}'.format(d.author, r_squared))
        #
        ax.plot(t_ok[mask], x_ok[mask],
                lw=10, alpha=0.5, color='tab:orange')
        ax.plot(t_ok[mask], result.best_fit, color='k', ls='--')
        # print(result.fit_report())

        par_str = ', '.join(['{:.1e}'.format(result.best_values[key])
                            for key in result.best_values.keys()])
        vs = Stokes_Velocity(diam, mu, rho_p, rho_f, g)  # [m/s]
        St = vs/u0
        ax.text(t_ok[mask][-1], x_ok[mask][-1],
                '{}: {}, {:.0f}, {:.0f}:{}, {:.1e} \n'.format(run,
                                                              phi, alpha, diam*1e6, d.particle_type, St) + par_str)
        # if func_fit == logistique:
        #     ax.text(t_ok[mask][-1]/t_ad, x_ok[mask][-1]/L0,
        #             '{}, {:.0f}, {:.0f}:{} \n {:.1e}, {:.1e}'.format(phi, alpha, diam*1e6, d.particle_type, p[0]/u0, p[1]/gprime), ha='left')
        # else:
        #     ax.text(t_ok[mask][-1]/t_ad, x_ok[mask][-1]/L0,
        #             '{}: {}, {:.0f}, {:.0f}:{} \n {:.0e}, {:.0e}, {:.0e}, {:.0e}'.format(run, phi, alpha, diam*1e6, d.particle_type, p[0]/u0, p[1]/(u0/t_ad), p[3]/(u0/t_ad**2), p[4]/(u0/t_ad**3)), ha='left')
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
