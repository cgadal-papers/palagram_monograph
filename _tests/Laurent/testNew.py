import os
import glob
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy.polynomial.polynomial as nppol


from scipy.optimize import curve_fit


def logistique(x, a, b, c):

    #y = a*x+c*np.log(b*x+1)
    y = a*(x-0) - b*(x-0)**(2)/2 + 0*x**3 + 0
    return y


def logistiqueNum(x, a, b, c):

    #y = a*x+c*np.log(b*x+1)
    y = a*(x-0) - b*(x-0)**(2)/2 + 0*x**3 + c
    return y


#plt.rcParams['text.usetex'] = True
plt.rcParams["font.serif"] = 'StixGeneral'
plt.rcParams["font.family"] = "serif"
plt.rcParams["mathtext.fontset"] = 'cm'

plt.rcParams['font.size'] = '22'

font = {'size': 28}


data_path = '../../data/input_data'

#list_runs = glob.glob(os.path.join(data_path, '0runs_MARIE/*.nc'))

#list_runs = glob.glob(os.path.join(data_path, 'SimusLEGIPALAGRAM/*.nc'))

list_runs = glob.glob(os.path.join(data_path, '*/*.nc'))

datasets = [Dataset(run) for run in list_runs]

RHO_F = []
vit = []
vitdim = []
atedim = []
atead = []
theta = []
jumanip = []
jmanip = []
mmanip = []
cmanip = []
cmmanip = []
ph = []
vsed = []
aspect = []
vstokes = []
At = []
Re = []
Stokes = []

fig, axarr = plt.subplots(
    1, 2, constrained_layout=True, sharex=True, sharey=True)
fig10, ax10 = plt.subplots(1, 1, figsize=(10, 6))
plt.subplots_adjust(left=.14, bottom=.15, right=0.9,
                    top=0.9, wspace=0, hspace=0)
#
for i, d in enumerate(datasets):
    # print(d.author if hasattr(d, 'author') else 'Marie')
    t = d.variables['t'][:].data
    x_front = d.variables['x_front'][:].data
    #
    if hasattr(d, 'author'):
        if d.author == ('Julien'):
            g = 9.81
            H = d.variables['H0'][:].data/100
            a = H/d.variables['L0'][:].data
            rho_f = d.variables['rho_f'][:].data*1000
            rho_p = d.variables['rho_p'][:].data*1000
            rho_a = d.variables['rho_a'][:].data*1000
            alpha = (d.variables['alpha'][:].data)*180/np.pi
            vs = (d.variables['d'][:].data*1e-6)**2*(rho_p-rho_f)*9.81/18/1.e-3
        else:
            g = 9.81
            H = d.variables['H0'][:].data
            a = H/d.variables['L0'][:].data
            rho_f = d.variables['rho_f'][:].data
            rho_p = d.variables['rho_p'][:].data
            rho_a = d.variables['rho_a'][:].data
            alpha = d.variables['alpha'][:].data
            vs = d.variables['d'][:].data**2 * \
                (d.variables['rho_p'][:].data-rho_f)*9.81/18/1.e-3
    else:
        g = 9.81
        H = d.variables['H0'][:].data
        a = H/d.variables['L0'][:].data
        rho_f = d.variables['rho_f'][:].data
        rho_p = d.variables['rho_p'][:].data
        rho_a = d.variables['rho_a'][:].data
        alpha = d.variables['alpha'][:].data
        vs = d.variables['d'][:].data**2 * \
            (d.variables['rho_p'][:].data-rho_f)*9.81/18/1.e-3

    rho_c = rho_f + \
        d.variables['phi'][:].data * (rho_p - rho_f)
    #
    gprime = g*(rho_c - rho_a)/rho_a

    #

    slope_corr = 0.25*np.cos(alpha*np.pi/180)+0.25*6*np.sin(alpha*np.pi/180)
    u0 = np.sqrt(gprime*H)
    t_ad = H/u0
    #
    u0 = np.sqrt(gprime*H)
    t_ad2 = H/u0

    # regression lineair

    tn = t[~np.isnan(x_front)]
    tn = tn[0:-2:1]
    xn = x_front[~np.isnan(x_front)]
    xn = xn[0:-2:1]

    if (any(tn)):
        if (d.variables['phi'][:].data > 0) & (tn[0] >= 0):
            if hasattr(d, 'author'):

                if d.author == ('Cyril/Marie'):
                    masque = (tn < 15*t_ad) & (tn > 3*t_ad)
                    resultat, covar = curve_fit(
                        logistique, tn[masque], xn[masque])
                    perr = np.sqrt(np.diag(covar))
                    aopt = resultat[0]
                    bopt = resultat[1]
                    copt = 0
                    ttest = tn[masque]
                    xtest = xn[masque]

                elif d.author == ('Cyril'):
                    masque = (tn < 15*t_ad) & (tn > 3*t_ad)
                #resultat = np.linalg.lstsq(tn[masque], xn[masque])
                #resultat = nppol.polyfit(tn[masque], xn[masque], 1)
                #aopt = resultat[1]
                    init_vals = [H*100, np.sqrt(gprime/H)/100]
                # plt.figure()
                #plt.plot(tn[masque], xn[masque])
                # plt.show()
                    resultat, covar = curve_fit(
                        logistique, tn[masque], xn[masque])
                    perr = np.sqrt(np.diag(covar))
                    aopt = resultat[0]
                    bopt = resultat[1]
                    copt = 0
                    ttest = tn[masque]
                    xtest = xn[masque]

                elif d.author == ('Jean'):
                    masque = (tn < 20*t_ad) & (tn > 0.5*t_ad)
                #resultat = np.linalg.lstsq(tn[masque], xn[masque])
                #resultat = nppol.polyfit(tn[masque], xn[masque], 1)
                #aopt = resultat[1]
                    init_vals = [H*100, np.sqrt(gprime/H)/100]
                # plt.figure()
                #plt.plot(tn[masque], xn[masque])
                # plt.show()
                    resultat, covar = curve_fit(
                        logistique, tn[masque], xn[masque])
                    perr = np.sqrt(np.diag(covar))
                    aopt = resultat[0]
                    bopt = resultat[1]
                    copt = 0
                    ttest = tn[masque]
                    xtest = xn[masque]

                elif d.author == ('Julien'):
                    masque = (tn < 15*t_ad) & (tn > 2*t_ad)
                #resultat = np.linalg.lstsq(tn[masque], xn[masque])
                #resultat = nppol.polyfit(tn[masque], xn[masque], 1)
                #aopt = resultat[1]
                    init_vals = [H*100, np.sqrt(gprime/H)/100]
                # plt.figure()
                #plt.plot(tn[masque], xn[masque])
                # plt.show()
                    resultat, covar = curve_fit(
                        logistiqueNum, tn[masque], xn[masque])
                    perr = np.sqrt(np.diag(covar))
                    aopt = resultat[0]
                    bopt = resultat[1]
                    copt = resultat[2]
                    ttest = tn[masque]
                    xtest = xn[masque]

                elif d.author == ('Rastello'):
                    masque = (tn < 60*t_ad) & (tn > 4*t_ad)
                #resultat = np.linalg.lstsq(tn[masque], xn[masque])
                #resultat = nppol.polyfit(tn[masque], xn[masque], 1)
                #aopt = resultat[1]
                    init_vals = [H*100, np.sqrt(gprime/H)/100]
                # plt.figure()
                #plt.plot(tn[masque], xn[masque])
                # plt.show()
                    resultat, covar = curve_fit(
                        logistique, tn[masque], xn[masque])
                    perr = np.sqrt(np.diag(covar))
                    if (perr[0] < 1e-5) or (t[-1]/t_ad > 30):
                        aopt = resultat[0]
                        bopt = resultat[1]
                        copt = 0
                        ttest = tn[masque]
                        xtest = xn[masque]
                    else:
                        masque = (tn < 13*t_ad) & (tn > 3*t_ad)
                        if (any(masque)):
                            resultat, covar = curve_fit(
                                logistique, tn[masque], xn[masque])
                            perr = np.sqrt(np.diag(covar))
                        # if (perr[0]<1):
                            aopt = resultat[0]
                            bopt = resultat[1]
                            copt = 0
                            ttest = tn[masque]
                            xtest = xn[masque]
                        # else:
                        #    aopt = 0
                        #    bopt = 0
                        #    copt = 0
                        #    ttest = []
                        #    xtest = []
                        else:
                            aopt = 0
                            bopt = 0
                            copt = 0
                            ttest = []
                            xtest = []
        else:
            aopt = 0
            bopt = 0
            copt = 0
    else:
        aopt = 0
        bopt = 0
        copt = 0

    # print(perr)

    #
    # if (x_front > 100).any():
    # x_front = x_front/100
    # print(d.author)
    # if d.variables['rho_p'][:].data in [1070, 1003, 1005, 1020]:
    if (d.variables['phi'][:].data > 0.001) & (d.variables['phi'][:].data < 0.6) & (alpha < 70) & (rho_p > 1000):
        # if (aopt>0):
        # vit.append(aopt/np.sqrt(gprime*H*np.cos(alpha*np.pi/180)))
        vit.append(aopt/np.sqrt(gprime*H))
        vitdim.append(aopt)
        aspect.append(H/d.variables['L0'][:].data)
        atedim.append(bopt)
        atead.append(bopt/gprime)
        theta.append(alpha)
        ph.append(d.variables['phi'][:].data)
        vsed.append(d.variables['v_s'][:].data)
        vstokes.append(vs)
        At.append((rho_c-rho_a)/rho_a)
        Re.append(np.sqrt(gprime*H)*H*rho_c/1.e-3)
        Stokes.append(vs/np.sqrt(gprime*H))
        if (alpha < 45):
            lw = 2
            ls = '-'
        else:
            lw = 1
            ls = '--'
        if hasattr(d, 'author'):
            if d.author == 'Jean':
                jumanip.append(0)
                jmanip.append(1)
                cmanip.append(0)
                cmmanip.append(0)
                mmanip.append(0)
                col = 'tab:blue'
            elif d.author == ('Cyril/Marie'):
                jumanip.append(0)
                jmanip.append(0)
                cmanip.append(0)
                cmmanip.append(1)
                mmanip.append(0)
                col = 'tab:gray'
            elif d.author == ('Cyril'):
                jumanip.append(0)
                jmanip.append(0)
                cmanip.append(1)
                cmmanip.append(0)
                mmanip.append(0)
                col = 'tab:cyan'
            elif d.author == ('Julien'):
                jumanip.append(1)
                jmanip.append(0)
                cmanip.append(0)
                cmmanip.append(0)
                mmanip.append(0)
                col = 'black'
            elif d.author == ('Rastello'):
                jumanip.append(0)
                jmanip.append(0)
                cmanip.append(0)
                cmmanip.append(0)
                mmanip.append(1)
                col = 'tab:green'
    else:
        ls = 'none'
        lw = 0
        col = 'r'
    # axarr[0].plot(t/t_ad, x_front/d.variables['L0'][:].data, lw=1, ls=ls)
    if not np.isnan(t/t_ad).all():
        axarr[0].plot(t/t_ad, x_front/H, lw=lw, ls=ls, color=col)
        axarr[1].plot(t, x_front, lw=lw, ls=ls, color=col)
    plt.xlabel(r'$t\sqrt{g^{*}/H_0}$')
    plt.ylabel(r'$x_f/H_0}$')

    #
    # figure pos dim

    if any(ttest) & (d.variables['phi'][:].data < 0.1) & (np.sin(alpha*np.pi/180) > 0.3) & (np.sin(alpha*np.pi/180) > 0.0) & ((i) % 1 == 0):
        print(perr)
        ax10.plot(tn/t_ad, xn/H, lw=1, ls='-', color=col)
        ax10.plot(ttest/t_ad, (aopt*ttest - bopt /
                  2*ttest*ttest + copt)/H, '--k')
        plt.xlabel(r'$t$ (s)', fontdict=font)
        plt.ylabel(r'$x_f$ (m)', fontdict=font)
        ax10.set(xlim=(0, 50), ylim=(0, 20))
        # print(ttest)
    #ax2.set_xscale("log", base=10)
    #ax2.set_yscale("log", base=10)

    RHO_F.append(H)

plt.show()

PHI = [d.variables['phi'][:].data for d in datasets]


# RHO_P = [d.variables['rho_p'][:].data for d in datasets]
# ALPHA = [d.variables['alpha'][:].data for d in datasets]
# RHO_A = [d.variables['rho_a']
#          [:].data for d in datasets if 'rho_f' in d.variables.keys()]

# RHO_F = [d.variables['rho_f']
#          [:].data for d in datasets if 'rho_f' in d.variables.keys()]

theta = np.array(theta)
aspect = np.array(aspect)
Re = np.array(Re)
At = np.array(At)
Stokes = np.array(Stokes)
vit = np.array(vit)
vitdim = np.array(vitdim)
atedim = np.array(atedim)
atead = np.array(atead)
ph = np.array(ph)
vsed = np.array(vsed)
vstokes = np.array(vstokes)
mmanip = np.array(mmanip)
cmanip = np.array(cmanip)
cmmanip = np.array(cmmanip)
jmanip = np.array(jmanip)
jumanip = np.array(jumanip)


anglem = theta[(mmanip == 1) & (vit > 0)]
vm = vit[(mmanip == 1) & (vit > 0)]
anglec = theta[(cmanip == 1) & (vit > 0)]
vc = vit[(cmanip == 1) & (vit > 0)]
anglecm = theta[(cmmanip == 1) & (vit > 0)]
vcm = vit[(cmmanip == 1) & (vit > 0)]
anglej = theta[(jmanip == 1) & (vit > 0)]
vj = vit[(jmanip == 1) & (vit > 0)]
angleju = theta[(jumanip == 1) & (vit > 0)]
vju = vit[(jumanip == 1) & (vit > 0)]

phim = ph[(mmanip == 1) & (vit > 0)]
phic = ph[(cmanip == 1) & (vit > 0)]
phicm = ph[(cmmanip == 1) & (vit > 0)]
phij = ph[(jmanip == 1) & (vit > 0)]
phiju = ph[(jumanip == 1) & (vit > 0)]
Stm = Stokes[(mmanip == 1) & (vit > 0)]
Stc = Stokes[(cmanip == 1) & (vit > 0)]
Stcm = Stokes[(cmmanip == 1) & (vit > 0)]
Stj = Stokes[(jmanip == 1) & (vit > 0)]
Stju = Stokes[(jumanip == 1) & (vit > 0)]

mph = max(ph)
am = ph/mph
alpham = (phim)/mph
alphac = (phic)/mph
alphacm = (phicm)/mph
alphaj = (phij)/mph
alphaju = (phiju)/mph


# fig vitesse fonction theta

fig, ax = plt.subplots(1, 1, figsize=(8, 7))
plt.subplots_adjust(left=.18, bottom=.15, right=0.9,
                    top=0.9, wspace=0, hspace=0)

for i, a in enumerate(alpham):
    # if (phim[i]>0.012) & (phim[i]<0.018):
    if (Stm[i] > 4e-2) & (Stm[i] < 6e-2):
        alpha = 1
    else:
        alpha = 0.2
    ax.plot(np.sin(anglem[i]*np.pi/180), vm[i], '.',
            markersize=16, color='tab:green', alpha=alpha)
ax.plot(np.sin(anglem[i]*np.pi/180), vm[i], '.',
        markersize=16, color='tab:green', alpha=1, label=r'LEGI')
for i, a in enumerate(alphac):
    # if (phic[i]>0.012) & (phic[i]<0.018):
    if (Stc[i] > 4e-2) & (Stc[i] < 6e-2):
        alpha = 1
    else:
        alpha = 0.2
    ax.plot(np.sin(anglec[i]*np.pi/180), vc[i], '.',
            markersize=16, color='tab:cyan', alpha=alpha)
ax.plot(np.sin(anglec[i]*np.pi/180), vc[i], '.',
        markersize=16, color='tab:cyan', alpha=1, label=r'IMFT')
for i, a in enumerate(alphacm):
    # if (phicm[i]>0.012) & (phicm[i]<0.018):
    if (Stcm[i] > 4e-2) & (Stcm[i] < 6e-2):
        alpha = 1
    else:
        alpha = 0.2
    ax.plot(np.sin(anglecm[i]*np.pi/180), vcm[i], '.',
            markersize=16, color='tab:gray', alpha=alpha)
ax.plot(np.sin(anglecm[i]*np.pi/180), vcm[i], '.',
        markersize=16, color='tab:gray', alpha=1, label=r'LEGI-IMFT')
for i, a in enumerate(alphaj):
    if (phij[i] > 0.012) & (phij[i] < 0.018):
        # if (Stj[i]>3e-2) & (Stj[i]<8e-2):
        alpha = 1
    else:
        alpha = 0.2
    ax.plot(np.sin(anglej[i]*np.pi/180), vj[i], '.',
            markersize=16, color='tab:red', alpha=alpha)
ax.plot(np.sin(anglej[i]*np.pi/180), vj[i], '.',
        markersize=16, color='tab:red', alpha=1, label=r'LEMTA')
for i, a in enumerate(alphaju):
    # if (phiju[i]>0.012) & (phiju[i]<0.018):
    if (Stju[i] > 4e-2) & (Stju[i] < 6e-2):
        alpha = 1
    else:
        alpha = 0.2
    ax.plot(np.sin(angleju[i]*np.pi/180), vju[i], 's',
            markersize=8, color='black', alpha=alpha)
ax.plot(np.sin(angleju[i]*np.pi/180), vju[i], 's',
        markersize=8, color='black', alpha=1, label=r'LEGI-Num')

# ax.plot([0,0.15],[0.44,0.44],'--k')
# ax.plot([0.07,0.75],[0.35,1.15],'--k')

ax.legend()
ax.legend(loc=2, prop={'size': 16})
plt.xlabel(r'$\sin{\alpha}$', fontdict=font)
plt.ylabel(r'$\mathcal{F}r$', fontdict=font)
ax.set(xlim=(-0.01, .8), ylim=(0, 1.5))
#ax.set_xscale("log", base=10)
#ax.set_yscale("log", base=10)

# figure vitesse fonction vs

anglec = 45
fig2, ax2 = plt.subplots(1, 1, figsize=(6, 3))
plt.subplots_adjust(left=.24, bottom=.25, right=0.9,
                    top=0.9, wspace=0, hspace=0)

ax2.plot([0, 1], [0.44, 0.44], '--k')

ax2.plot(vstokes[(cmanip == 1) & (theta > anglec-2) & (theta < anglec+2)]*vit[(cmanip == 1) & (theta > anglec-2) & (theta < anglec+2)]/vitdim[(cmanip == 1) & (theta >
         anglec-2) & (theta < anglec+2)], vit[(cmanip == 1) & (theta > anglec-2) & (theta < anglec+2)], '.', markersize=16, color='tab:cyan', label=r'IMFT $\alpha\approx 7^o$')
ax2.plot(vstokes[(cmmanip == 1) & (theta > anglec-2) & (theta < anglec+2)]*vit[(cmmanip == 1) & (theta > anglec-2) & (theta < anglec+2)]/vitdim[(cmmanip == 1) & (theta > anglec-2)
         & (theta < anglec+2)], vit[(cmmanip == 1) & (theta > anglec-2) & (theta < anglec+2)], '.', markersize=16, color='tab:gray', label=r'IMFT-LEGI $\alpha\approx 7^o$')
ax2.plot(vstokes[(mmanip == 1) & (theta > anglec-2) & (theta < anglec+2)]*vit[(mmanip == 1) & (theta > anglec-2) & (theta < anglec+2)]/vitdim[(mmanip == 1) & (theta >
         anglec-2) & (theta < anglec+2)], vit[(mmanip == 1) & (theta > anglec-2) & (theta < anglec+2)], '.', markersize=16, color='tab:green', label=r'LEGI $\alpha\approx 7^o$')
ax2.plot(vstokes[(jumanip == 1) & (theta > anglec-2) & (theta < anglec+2)]*vit[(jumanip == 1) & (theta > anglec-2) & (theta < anglec+2)]/vitdim[(jumanip == 1) & (theta >
         anglec-2) & (theta < anglec+2)], vit[(jumanip == 1) & (theta > anglec-2) & (theta < anglec+2)], 's', markersize=8, color='black', label=r'LEGI-Num $\alpha\approx 7^o$')

# ax2.legend()
#ax2.legend(loc=2, prop={'size': 16})
plt.xlabel(r'$St$', fontdict=font)
plt.ylabel(r'$\mathcal{F}r$', fontdict=font)
ax2.set(xlim=(0.02, 200), ylim=(0., 2))
ax2.set_xscale("log", base=10)
#ax2.set_yscale("log", base=10)

# figure attenuation vs vs

fig3, ax3 = plt.subplots(1, 1, figsize=(6, 3))
plt.subplots_adjust(left=.24, bottom=.25, right=0.9,
                    top=0.9, wspace=0, hspace=0)

ax3.plot([0, .04], [0., 0.], '--k')
x = np.logspace(-5, 0, 100)
ax3.plot(x-0.0, np.sqrt(x-0.04)/10, '--k')

ax3.plot(vstokes[(cmanip == 1) & (theta > anglec-2) & (theta < anglec+2)]*vit[(cmanip == 1) & (theta > anglec-2) & (theta < anglec+2)]/vitdim[(cmanip == 1) & (theta > anglec-2)
         & (theta < anglec+2)], (atead[(cmanip == 1) & (theta > anglec-2) & (theta < anglec+2)]), '.', markersize=16, color='tab:cyan', label=r'IMFT $\alpha\approx 7^o$')
ax3.plot(vstokes[(cmmanip == 1) & (theta > anglec-2) & (theta < anglec+2)]*vit[(cmmanip == 1) & (theta > anglec-2) & (theta < anglec+2)]/vitdim[(cmmanip == 1) & (theta > anglec-2)
         & (theta < anglec+2)], (atead[(cmmanip == 1) & (theta > anglec-2) & (theta < anglec+2)]), '.', markersize=16, color='tab:gray', label=r'IMFT-LEGI $\alpha\approx 7^o$')
ax3.plot(vstokes[(mmanip == 1) & (theta > anglec-2) & (theta < anglec+2)]*vit[(mmanip == 1) & (theta > anglec-2) & (theta < anglec+2)]/vitdim[(mmanip == 1) & (theta > anglec-2)
         & (theta < anglec+2)], (atead[(mmanip == 1) & (theta > anglec-2) & (theta < anglec+2)]), '.', markersize=16, color='tab:green', label=r'LEGI $\alpha\approx 7^o$')
ax3.plot(vstokes[(jumanip == 1) & (theta > anglec-2) & (theta < anglec+2)]*vit[(jumanip == 1) & (theta > anglec-2) & (theta < anglec+2)]/vitdim[(jumanip == 1) & (theta >
         anglec-2) & (theta < anglec+2)], (atead[(jumanip == 1) & (theta > anglec-2) & (theta < anglec+2)]), 's', markersize=8, color='black', label=r'LEGI-Num $\alpha\approx 7^o$')

# ax3.legend()
#ax3.legend(loc=2, prop={'size': 16})
plt.xlabel(r'$St$', fontdict=font)
plt.ylabel(r'$\tilde{\lambda}$', fontdict=font)
ax3.set(xlim=(0.02, 200), ylim=(-0.1, .1))
ax3.set_xscale("log", base=10)
#ax3.set_yscale("log", base=10)

# figure vitess vs phi

fig4, ax4 = plt.subplots(1, 1, figsize=(6, 6))
plt.subplots_adjust(left=.24, bottom=.15, right=0.9,
                    top=0.9, wspace=0, hspace=0)

ax4.plot(ph[(mmanip == 1) & (theta > anglec-2) & (theta < anglec+2) & (vit > 0)], vit[(mmanip == 1) & (theta >
         anglec-2) & (theta < anglec+2) & (vit > 0)], '.', markersize=14, color='tab:green', label=r'LEGI $\alpha=44^o$')
#ax4.plot(ph[(cmanip==1)],vit[(cmanip==1)],'.',markersize=14,color='tab:blue', label=r'IMFT')
ax4.plot(ph[(jmanip == 1)], vit[(jmanip == 1)], '.', markersize=14,
         color='tab:blue', label=r'LEMTA $\alpha=0$')

ax4.legend()
ax4.legend(loc=1, prop={'size': 16})
plt.xlabel(r'$\phi$', fontdict=font)
plt.ylabel(r'$U_c/U_0$', fontdict=font)
ax4.set(xlim=(0.01, 1), ylim=(0., 2))
ax4.set_xscale("log", base=10)
#ax4.set_yscale("log", base=10)

# figure atte vs phi

fig5, ax5 = plt.subplots(1, 1, figsize=(6, 6))
plt.subplots_adjust(left=.24, bottom=.15, right=0.9,
                    top=0.9, wspace=0, hspace=0)

ax5.plot(ph[(mmanip == 1) & (theta > anglec-2) & (theta < anglec+2) & (vit > 0)], atead[(mmanip == 1) & (theta >
         anglec-2) & (theta < anglec+2) & (vit > 0)], '.', markersize=14, color='tab:green', label=r'LEGI $\alpha=44^o$')
#ax4.plot(ph[(cmanip==1)],vit[(cmanip==1)],'.',markersize=14,color='tab:blue', label=r'IMFT')
ax5.plot(ph[(jmanip == 1)], atead[(jmanip == 1)], '.',
         markersize=14, color='tab:blue', label=r'LEMTA $\alpha=0$')

ax5.legend()
ax5.legend(loc=1, prop={'size': 16})
plt.xlabel(r'$\phi$', fontdict=font)
plt.ylabel(r'$\lambda/g^{,}$', fontdict=font)
ax5.set(xlim=(0.01, 1), ylim=(-0.1, 0.1))
ax5.set_xscale("log", base=10)
#ax4.set_yscale("log", base=10)

# figure test

fig6, ax6 = plt.subplots(1, 1, figsize=(6, 6))
plt.subplots_adjust(left=.24, bottom=.15, right=0.9,
                    top=0.9, wspace=0, hspace=0)

ax6.plot(Stokes, '.', markersize=14, color='tab:green',
         label=r'LEGI $\alpha=44^o$')
ax6.set_yscale("log", base=10)
