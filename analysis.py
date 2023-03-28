import numpy as np
import glob
import os

import numpy as np
from netCDF4 import Dataset
from scipy.optimize import curve_fit


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


def logistiqueNUM(x, a, b, c):
    y = a*(x-0) - b*(x-0)**(2)/2 + c
    return y


def logistique(x, a, b):
    return logistiqueNUM(x, a, b, c=0)


def compute_rsquared(ydata, yfit, ndeg):
    residuals = ydata - yfit
    ss_res = np.sum(residuals**2)
    ss_tot = np.sum((ydata-np.nanmean(ydata))**2)
    r2 = 1 - (ss_res / ss_tot)
    return r2


def Stokes_Velocity(d, mu, rho_p, rho_f, g):
    return d**2*(rho_p - rho_f) * g/mu/18


# %% Variable definition
# paths
input_path = 'data/input_data'
output_path = 'data/output_data'

# physical parameters
g = 9.81  # [m/s2]
mu = 1e-3  # [kg/m/s]

# fit specifications
BOUNDS_FIT = {'Cyril': (3, 15), 'Cyril/Marie': (3, 15), 'Jean': (0.5, 20),
              'Julien': (2, 15), 'Rastello': (4, 60)}

FUNC_FIT = {'Cyril': logistique, 'Cyril/Marie': logistique, 'Jean': logistique,
            'Julien': logistiqueNUM, 'Rastello': logistique}

# other
SETUPS = {'Cyril': 'IMFT', 'Cyril/Marie': 'LEGI', 'Jean': 'LEMTA',
          'Julien': 'NUM', 'Rastello': 'LEGI'}

# %% Loading data
list_runs = glob.glob(os.path.join(input_path, 'runs*/*.nc'))
datasets = [Dataset(run) for run in list_runs]

# %% Loop over data file and analysis

for i, d in enumerate(datasets):
    t = d.variables['t'][:].data
    x_front = d.variables['x_front'][:].data
    #
    # Loading some variables
    H0 = d.variables['H0'][:].data        # lock characteristic height, [cm]
    rho_p = d.variables['rho_p'][:].data  # particle velocity, [kg/m3]
    rho_f = d.variables['rho_f'][:].data  # lock fluid density, [kg/m3]
    rho_a = d.variables['rho_a'][:].data  # ambiant fluid density, [kg/m3]
    alpha = d.variables['alpha'][:].data  # bottom slope, [deg.]
    diam = d.variables['d'][:].data*1e-6  # grain size, [mum]
    phi = d.variables['phi'][:].data
    if d.author == 'Julien':
        H0 = H0/100
        rho_f, rho_p, rho_a = rho_f*1e3, rho_p*1e3, rho_a*1e3
        alpha = alpha*180/np.pi
    #
    # Computing other variables
    a = H0/d.variables['L0'][:].data  # lock aspect ratio
    vs = Stokes_Velocity(diam, mu, rho_p, rho_f, g)  # [m/s]
    rho_c = rho_f + phi * (rho_p - rho_f)  # average lock density, [kg/m3]
    gprime = g*(rho_c - rho_a)/rho_a  # specific gravity
    u0 = np.sqrt(gprime*H0)
    t_ad = H0/u0
    #
    #
    # #### Fitting front position curves
    # defining fitting masks
    mask_ok = ~np.isnan(x_front)
    t_ok, x_ok = t[mask_ok][:-2], x_front[mask_ok][:-2]
    bounds_fit = BOUNDS_FIT[d.author]
    mask = (t_ok > bounds_fit[0]*t_ad) & (t_ok < bounds_fit[1]*t_ad)
    #
    # Make fit
    func_fit = FUNC_FIT[d.author]
    if mask.any():
        p, pcov = curve_fit(func_fit, t_ok[mask], x_ok[mask])
        perr = np.sqrt(np.diag(pcov))
        r_squared = compute_rsquared(
            x_ok[mask], func_fit(t_ok[mask], *p), p.size)
        #
        print('author: {}, r2: {:.3f}'.format(d.author, r_squared))
        if (d.author == 'Rastello') & ~((perr[0] < 1e-5) or (t[-1]/t_ad > 30)):
            mask = (t_ok < 13*t_ad) & (t_ok > 3*t_ad)
            p, pcov = curve_fit(func_fit, t_ok[mask], x_ok[mask])
            perr = np.sqrt(np.diag(pcov))
            r_squared = compute_rsquared(
                x_ok[mask], func_fit(t_ok[mask], *p), p.size)
            print('Bad Fit, new r2: {:.3f}'.format(r_squared))
    else:
        p = np.array([np.nan, np.nan, np.nan])
        perr = np.copy(p)
        r_squared = np.nan
    #
    #
    # ### Writting corresponding netcdf file
    path_dataset = os.path.join(output_path, 'run_{:03d}.nc'.format(i))
    # creating netcdf file and groups
    newfile = Dataset(path_dataset, "w", format="NETCDF4")
    newfile.setup = SETUPS[d.author]
    # copy input data
    newfile.setncatts(d.__dict__)  # global attributes
    for name, dimension in d.dimensions.items():  # dimensions
        newfile.createDimension(
            name, (len(dimension) if not dimension.isunlimited() else None))
    for name, variable in d.variables.items():  # variables
        create_variable(
            newfile, name, d.variables[name][:], dimensions=variable.dimensions, type=variable.datatype)
        newfile[name].setncatts(d[name].__dict__)  # copy variable attributes
    # correct Julien stuff
    if newfile.author == 'Julien':
        newfile.variables['H0'][:] = newfile.variables['H0'][:].data/100
        newfile.variables['alpha'][:] = newfile.variables['alpha'][:].data*180/np.pi
        for var in ['rho_p', 'rho_f', 'rho_a']:
            newfile.variables[var][:] = newfile.variables[var][:].data*1e3
    # fit results
    fitdim = newfile.createDimension("fitdim", 3)
    if p.size == 2:
        p = np.append(p, 0)
        perr = np.append(perr, 0)
    create_variable(newfile, 'p', p, std=perr, unit=['cm/s', 'cm/s2', 'cm'],
                    comments='fit results of logistic curve', dimensions=(fitdim))
    # other variables
    create_variable(newfile, 'a', a, unit='None', comments='lock aspect ratio')
    create_variable(newfile, 'vs', vs, unit='m/s',
                    comments='particle Stokes velocity')
    create_variable(newfile, 'rho_c', rho_c, unit='kg/m3',
                    comments='lock average density')
    create_variable(newfile, 'gprime', gprime, unit='m/s2',
                    comments='specific gravity')
    create_variable(newfile, 'u0', u0, unit='cm/s',
                    comments='characteristic velocity')
    # saving and closing netcdf file
    newfile.close()
