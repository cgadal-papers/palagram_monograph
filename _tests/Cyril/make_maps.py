import glob
import os

import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import Dataset


def compute_drag_coeffs(P, csf):
    """
    Drag coefficients from
    Simple and General Formula for the Settling Velocity of Particles,
    Camenen B., 2007
    """
    a1, a2 = 24, 100
    b1, b2 = 0.39 + 0.22*(6-P), 20
    m1, m2 = 1.2 + 0.12*P, 0.47
    a3, b3 = 2.1 + 0.06*P, 1.75 + 0.35*P
    A = a1 + a2*(1 - np.sin(np.pi*csf/2))**a3
    B = b1 + b2*(1 - np.sin(np.pi*csf/2))**b3
    m = m1*np.sin(np.pi*csf/2)**m2
    return A, B, m


def Generalized_settling_velocity(phi, d, rho_s, m=1, A=24.4, B=0.4,
                                  nu_f=1e-6, rho_f=0.998, g=9.8, alpha=3.1):
    """
    Simple and General Formula for the Settling Velocity of Particles,
    Camenen B., 2007
    -> balancing buoyancy with drag, Cd = ((A/Rp)**(1/m) + B**(1/m))**m
    default coeffs A, B, m are from Richardson and Zacki
    """
    s = rho_s/rho_f
    d_star = (1 - phi)**(alpha/3)*d*((s - 1)*g/nu_f**2)**(1/3)  # -> alpha/3 ?
    R_star = (np.sqrt((1/4)*(A/B)**(2/m)
                      + ((4/3)*(d_star**3/B))**(1/m))
              - (1/2)*(A/B)**(1/m))**(m)
    return nu_f*R_star/d


# ## general parameters
data_path = '../data'

g = 9.81  # gravity, [m/s2]
nu_f_0 = 1e-6  # kinmatic viscosity, [m2/s]

# ## Load datas
list_runs = glob.glob(os.path.join(data_path, '*/*.nc'))
datasets = [Dataset(run) for run in sorted(list_runs)]
authors = np.array([d.author if hasattr(d, 'author')
                   else 'Marie' for d in datasets])


# ## Get variables

# lock geometry
H0 = np.array([d.variables['H0'][:].data for d in datasets])
L0 = np.array([d.variables['L0'][:].data for d in datasets])
W0 = np.array([d.variables['W0'][:].data for d in datasets])

# lock content
rho_f = np.array([d.variables['rho_f'][:].data for d in datasets])
rho_p = np.array([d.variables['rho_p'][:].data for d in datasets])
d = np.array([d.variables['d'][:].data for d in datasets])
phi = np.array([d.variables['phi'][:].data for d in datasets])
vs = np.array([d.variables['v_s'][:].data for d in datasets])
nu_f = np.array([d.variables['nu_f'][:].data for d in datasets])
T_f = np.array([d.variables['T_f'][:].data for d in datasets])

# ambient
alpha = np.array([d.variables['alpha'][:].data for d in datasets])
H_a = np.array([d.variables['alpha'][:].data for d in datasets])
rho_a = np.array([d.variables['rho_a'][:].data for d in datasets])
nu_a = np.array([d.variables['nu_a'][:].data for d in datasets])
T_a = np.array([d.variables['T_a'][:].data for d in datasets])


# #### correcting quantities
# nu_f --> Cyril = Nans, Jean = 1e-3, Marie = 1e-6
nu_f[:] = nu_f_0
# V_s
A, B, m = compute_drag_coeffs(6, 1)  # spherical particles
mask = np.isnan(vs)
vs[mask] = Generalized_settling_velocity(0, d[mask], rho_p[mask],
                                         m=m, A=A, B=B, nu_f=nu_f[mask], rho_f=rho_f[mask])
# phi = 0 for Saline currents
phi[np.isnan(phi) & (authors == 'Cyril')] = 0


# #### New quantities

rho_m = rho_f + phi*(rho_p - rho_f)
delta_rho = (rho_m - rho_a)/rho_a
u0 = np.sqrt(H0*delta_rho*g)
Reynolds = u0*H0/nu_f
Stokes = vs/u0
