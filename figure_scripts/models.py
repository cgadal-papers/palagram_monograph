import numpy as np

# %% Birman et al 2007


def Birman(ang):
    return -0.1824*np.radians(ang)**2 + 0.2781*np.radians(ang) + 0.4871


# %% energetic bilan model (extension of Gadal et al. 2023)

def Froude(theta, eta, Re, Fr0=0.5, a=1, r=1, Cd=0.4, E=500, A=3):
    """
    Calculate $Fr = U_c/u_0$, with $u_0 = \sqrt{\cos\theta\delta_\rho g h_0/rho_c}$ using an energetic bilan similar to that of Gadal et al. 2023. 

    Parameters
    ----------
    Fr0 : scalar, array_like
        Limit for $r=1$, $theta \to 0$, $Cd \to 0$, $Re \to \infty$.
    theta : scalar, array_like
        Bottom slope, in degree.
    eta : scalar, array_like
        Non-dimensinal viscosity $\eta/\eta_f$.
    Re : scalar, array_like
        Reynlods number $Re = \rho_0 u_0 h_{0}/eta_{f}$
    a : scalar, array_like
        Lock aspect ratio, $a = h_0/L_0$
    r : scalar, array_like
        Density ratio, $r = \rho_0/rho_a$
    Cd : scalar, array_like
        Drag coefficient
    E : scalar, array_like
        Viscous dissipation coefficient
    A : scalar, array_like
        Slope effect coefficient

    Returns
    -------
    scalar, array_like
        Non-dimensional velocity $Fr = U_c/u_0$, with $u_0 = \sqrt{\cos\theta\delta_\rho g h_0/rho_c}$.
    """
    drag = (1 + Cd*a/r)
    viscous = E*eta/Re
    Fr = (np.sqrt(viscous**2 + Fr0**2*drag*(1 +
          A*np.tan(np.radians(theta)))) - viscous)/drag
    return Fr

# %% Viscosity models


def Krieger_viscosity(phi, phi_m=0.585, eta_cr=5/2):
    # (Krieger and Dougherty (1959) model, from Stickel and Powell (2005))
    return (1 - phi/phi_m)**(-eta_cr*phi_m)


def Boyer_viscosity(phi, phi_m=0.585, mu1=0.32, mu2=0.7, I0=0.005):
    # Unifying Suspension and Granular Rheology, Boyer et al. 2010
    mu_c = mu1 + (mu2 - mu1)/(1 + I0*phi**2*(phi_m - phi)**(-2))
    return (1 + (5/2)*phi*(1 - phi/phi_m) + mu_c*(phi/(phi_m - phi))**2)


def Ferrini_viscosity(phi, phi_m=0.585, eta_cr=5/2):
    # (Ferrini et al. (1979) model, from Stickel and Powell (2005))
    return (1 + 0.5*eta_cr*phi/(1 - phi/phi_m))**2


if __name__ == '__main__':
    phi = np.linspace(0, 0.58, 200)
    import matplotlib.pyplot as plt
    plt.figure()
    plt.plot(phi, Krieger_viscosity(phi), label='Krieger')
    plt.plot(phi, Boyer_viscosity(phi), label='Boyer')
    plt.plot(phi, Ferrini_viscosity(phi), label='Ferrini')
    plt.gca().set_yscale('log')
    plt.legend()
    plt.show()
