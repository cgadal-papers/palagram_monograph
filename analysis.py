import glob
import os
import shutil

import numpy as np
from lmfit import Model
from lmfit.model import save_modelresult
from netCDF4 import Dataset


def create_variable(netcdf_group, name, data, dimensions=None, std=None,
                    unit=None, comments=None, type='float64'):
    if dimensions is not None:
        var = netcdf_group.createVariable(name, type, (dimensions))
    else:
        var = netcdf_group.createVariable(name, type)
    var[:] = data
    if std is not None:
        var.std = std
    if unit is not None:
        var.unit = unit
    if comments is not None:
        var.comments = comments


def polyFULL(t, Fr, lamb, xi, c, d):
    return xi + Fr*(t - lamb*t**2) + c*t**3 + d*t**4


def Stokes_Velocity(d, mu, rho_p, rho_f, g):
    return d**2*(rho_p - rho_f) * g/mu/18


def determine_fit_props(author, alpha, run, params):
    params['xi'].vary = False
    params['Fr'].vary = True
    params['lamb'].vary = True
    params['c'].vary = False
    params['d'].vary = False
    t_bounds = [5, 30]
    #
    if author == 'Jean':
        t_bounds = [0, 40]
        params['lamb'].vary = False
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
input_path = 'data/input_data'
output_path = 'data/output_data'
shutil.rmtree(output_path)
os.makedirs(output_path)

# physical parameters
g = 9.81  # [m/s2]
mu = 1e-3  # [kg/m/s]

# other
SETUPS = {'Cyril': 'IMFT', 'Cyril/Marie': 'LEGI', 'Jean': 'LEMTA',
          'Julien': 'NUM', 'Rastello': 'LEGI'}

# %% fit objects definition
# model object creation
model = Model(polyFULL)
params = model.make_params()

# parameter properties (Non dim.)
p0 = {'Fr': 0.4, 'xi': 0, 'lamb': 0, 'c': 0, 'd': 0}

# lower_bounds = {'Fr': 0, 'xi': -np.inf,
#                 'L': -0.05, 'c': -1e-3, 'd': -1e-4}
# higher_bounds = {'Fr': 2, 'xi': np.inf, 'L': 0.1, 'c': 1e-3, 'd': 1e-4}

lower_bounds = {'Fr': 0, 'xi': -np.inf,
                'lamb': -0.1, 'c': -1e-3, 'd': -1e-4}
higher_bounds = {'Fr': 2, 'xi': np.inf, 'lamb': 0.2, 'c': 1e-3, 'd': 1e-4}

# set parameter bounds
for par in params.keys():
    params[par].set(value=p0[par], min=lower_bounds[par],
                    max=higher_bounds[par])

# %% Loading data
list_runs = sorted(glob.glob(os.path.join(input_path, 'runs*/*.nc')))
datasets = [Dataset(run) for run in list_runs]

# %% Loop over data file and analysis

for i, d in enumerate(datasets):
    run = list_runs[i].split(os.sep)[-1]
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
        if d.variables['d'].unit == 'mum':
            print(run, '--d')
            diam = diam*1e-6  # grain size, [m]
        if d.variables['rho_f'].unit == 'g/cm3':
            print(run, '-rho')
            rho_f, rho_p, rho_a = rho_f*1e3, rho_p*1e3, rho_a*1e3
        if d.variables['H0'].unit == 'cm':
            print(run, '-L0, H0')
            H0 = H0/100
            L0 = L0/100
    #
    # Computing variables for adi time
    rho_c = rho_f + phi * (rho_p - rho_f)  # average lock density, [kg/m3]
    gprime = g*(rho_c - rho_a)/rho_a  # specific gravity
    u0 = np.sqrt(gprime*H0*np.cos(np.radians(alpha)))
    t_ad = L0/u0
    #
    # #### Fitting front position curves
    t_bounds = determine_fit_props(d.author, alpha, run, params)
    # defining fitting masks
    mask_ok = ~np.isnan(x_front)
    t_ok, x_ok = t[mask_ok]/t_ad, x_front[mask_ok]/L0
    mask = (t_ok > t_bounds[0]) & (t_ok < t_bounds[1])
    #
    # Make fit
    if mask.sum() > 5:
        result = model.fit(x_ok[mask], params, t=t_ok[mask])
        lamb, lamb_err = result.params['lamb'].value, result.params['lamb'].stderr
        if lamb > 1.2e-2:
            params['d'].vary = True
            result = model.fit(x_ok[mask], params, t=t_ok[mask])
        r_squared = result.rsquared
        print('author: {}, r2: {:.3f}'.format(d.author, r_squared))
        Fr, Fr_err = result.params['Fr'].value, result.params['Fr'].stderr
        lamb, lamb_err = result.params['lamb'].value, result.params['lamb'].stderr
        if not params['lamb'].vary:
            lamb, lamb_err = np.nan, np.nan
    else:
        Fr, Fr_err = np.nan, np.nan
        lamb, lamb_err = np.nan, np.nan
        print('All NaNs')
    #
    # ### Writting corresponding netcdf file
    path_dataset = os.path.join(output_path, 'run_{:03d}.nc'.format(i))
    # creating netcdf file and groups
    newfile = Dataset(path_dataset, "w", format="NETCDF4")
    newfile.setup = SETUPS[d.author]
    if not hasattr(d, 'run_oldID'):
        newfile.run_oldID = '{} -- {}'.format(d.author, run)
    # copy input data
    newfile.setncatts(d.__dict__)  # global attributes
    for name, dimension in d.dimensions.items():  # dimensions
        newfile.createDimension(
            name, (len(dimension) if not dimension.isunlimited() else None))
    newfile.label = 'run_{:03d}.nc'.format(i)
    for name, variable in d.variables.items():  # variables
        create_variable(
            newfile, name, d.variables[name][:], dimensions=variable.dimensions, type=variable.datatype)
        newfile[name].setncatts(d[name].__dict__)  # copy variable attributes
    # correct Julien stuff
    if newfile.author == 'Julien':
        newfile.variables['alpha'][:] = newfile.variables['alpha'][:].data*180/np.pi
        if newfile.variables['d'].unit == 'mum':
            newfile.variables['d'][:] = newfile.variables['d'][:].data*1e-6
            newfile.variables['d'].unit = 'm'
        if d.variables['rho_f'].unit == 'g/cm3':
            for var in ['rho_p', 'rho_f', 'rho_a']:
                newfile.variables[var][:] = newfile.variables[var][:].data*1e3
                newfile.variables[var].unit = 'kg/m3'
        if d.variables['H0'].unit == 'cm':
            for var in ['H0', 'L0', 'W0']:
                newfile.variables[var][:] = newfile.variables[var][:].data/100
                newfile.variables[var].unit = 'm'

    # other variables
    vs = Stokes_Velocity(diam, mu, rho_p, rho_f, g)  # [m/s]
    create_variable(newfile, 'vs', vs, unit='m/s',
                    comments='particle Stokes velocity')
    create_variable(newfile, 'rho_c', rho_c, unit='kg/m3',
                    comments='lock average density')
    create_variable(newfile, 'gprime', gprime, unit='m/s2',
                    comments='specific gravity')
    create_variable(newfile, 'u0', u0, unit='m/s',
                    comments='characteristic velocity')
    create_variable(newfile, 't0', t_ad, unit='s',
                    comments='characteristic timescale')
    # Non-dimensional numbers
    create_variable(newfile, 'a', H0/d.variables['L0'][:].data, unit='-',
                    comments='lock aspect ratio')
    create_variable(newfile, 'Re', u0*H0*rho_c/mu, unit='-',
                    comments='Reynolds number')
    create_variable(newfile, 'At', (rho_c - rho_a)/rho_a, unit='-',
                    comments='Atwood number')
    create_variable(newfile, 'St', vs/u0, unit='-',
                    comments='Stokes number')
    # Non-dimensional variables and fit results
    create_variable(newfile, 'Fr', Fr, unit='-', std=Fr_err,
                    comments='Froude number (adi. initial current velocity)')
    create_variable(newfile, 'lamb', lamb, unit='-', std=lamb_err,
                    comments='adi. attenuation parameter')

    # saving and closing netcdf file
    newfile.close()
    # fit results
    if mask.sum() > 5:
        save_modelresult(result, os.path.join(output_path,
                                              'fitresult_run_{:03d}.save'.format(i)))
