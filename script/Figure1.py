import os
import numpy as np
from astropy import units as u
from astropy import constants as const
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
from matplotlib import rc

plt.style.use('figures.mplstyle')
plt.rcParams.update({'figure.figsize': (12, 6)})

basepath = "../data/DIAD/"
files = [
    "fort14.a0.01.irr.grid_001.e1.0.amax0p25_amax1mm.dat",
    "fort14.a0.01.irr.grid_002.e1.0.amax0p25_amax1mm.dat",
    "fort14.a0.01.irr.grid_003.e1.0.amax0p25_amax1mm.dat",
    "fort14.a0.01.irr.grid_004.e1.0.amax0p25_amax1mm.dat",
]
Rstars = [1.055, 2.31, 2.125, 2.615]
Mstars = [0.1, 0.3, 0.5, 1.0]

fig = plt.figure()
gs = fig.add_gridspec(nrows=2, ncols=4, hspace=0, wspace=0)
ax = gs.subplots(sharex=False, sharey=True)

labels = ["\\texttt{0.1Msun}", "\\texttt{0.3Msun}", "\\texttt{0.5Msun}", "\\texttt{1Msun}"]

xticks = np.arange(0, 21, 5)
xtick_label = [[0,2.5,5,7.5,10],
          [0,8,16,24,32],
          [0,14,28,42,56],
          [0,25,50,75,100]]

for i in range(2):
    for j in range(4):
        ax[i, j].set_xlim(0, 24)
        ax[i, j].set_ylim(0, 24)
        ax[i, j].set_aspect("equal")

        Rg = 5.05 * Mstars[j]

        if i == 1:
            ax[i, j].set_xlabel(r"R [R$_g$]")
            ax[i, j].set_xticks(xticks)
        if i == 0:
            ax[i, j].set_title(labels[j], size=14.0)
            ax2 = ax[i,j].twiny()
            ax[i, j].set_xticks([])
            ax2.set_xlim(ax[i,j].get_xlim())
            ax2.set_xticks(np.asarray(xtick_label[j])/Rg)
            ax2.set_xticklabels(xtick_label[j])
            ax2.set_xlabel(r'R [au]')

        if j == 0:
            ax[i, j].set_ylabel(r"Z [R$_g$]", size=14, color="black")
        file = basepath + files[j]
        NR = 79
        NY = 2
        NZ = 480
        R = np.zeros((NR, NZ))
        z = np.zeros((NR, NZ))
        T = np.zeros((NR, NZ))
        d = np.zeros((NR, NZ))
        d1 = np.zeros((NR, NZ))
        Rstar = Rstars[j]
        Mstar = Mstars[j]
        l_unit = (Rstar * const.R_sun).to(u.cm).value

        with open(file) as fp:
            # skip first 2 lines
            fp.readline()
            nr = 0
            while nr < NR:
                nz = NZ - 1
                fp.readline()
                # read R from the 3rd
                line = fp.readline()
                R[nr, :] = line.split()[0]
                fp.readline()
                fp.readline()
                line = fp.readline()
                while line:
                    if nz > 0:
                        z[nr, nz] = line.split()[0]
                        d1[nr, nz] = line.split()[1]
                        T[nr, nz] = line.split()[2]
                        d[nr, nz] = line.split()[3]
                        line = fp.readline()
                        nz -= 1
                    else:
                        z[nr, nz] = line.split()[0]
                        d1[nr, nz] = line.split()[1]
                        T[nr, nz] = line.split()[2]
                        d[nr, nz] = line.split()[3]
                        break
                nr += 1
        fp.close()

        R *= l_unit
        z *= l_unit

        numr = 2000
        numz = 7000
        xval = np.linspace(min(R[:, 0]), max(R[:, 0]), numr, endpoint=True)
        zval = np.linspace(min(z[:, 1]), max(z[:, NZ - 1]), numz)
        Xc, Zc = np.meshgrid(xval, zval)

        method = "linear"
        Dc = griddata(
            (R.reshape(480 * 79), z.reshape(480 * 79)),
            d.reshape(480 * 79),
            (Xc.reshape(numz * numr), Zc.reshape(numz * numr)),
            method=method,
            fill_value=1.0e-40,
        )
        Dc = Dc.reshape(numz, numr)
        Tdustc = griddata(
            (R.reshape(480 * 79), z.reshape(480 * 79)),
            T.reshape(480 * 79),
            (Xc.reshape(numz * numr), Zc.reshape(numz * numr)),
            method=method,
            fill_value=10.0,
        )
        Tdustc = Tdustc.reshape(numz, numr)

        if i == 0:
            quantity = Tdustc
            labelplot = r"T [K]"
            valmin = 0
            valmax = 2500
        else:
            quantity = np.log10(Dc / 10 ** (-24))
            labelplot = r"$\log_{10}(\rho)$ [$10^{-24}$ g cm$^{-3}$]"
            valmin = 3
            valmax = 11

        plot = ax[i, j].pcolormesh(
            Xc * u.cm.to(u.AU) / Rg,
            Zc * u.cm.to(u.AU) / Rg,
            quantity,
            shading="auto",
            vmin=valmin,
            vmax=valmax,
            cmap=plt.cm.jet,
        )

        if j == 3:
            fig.subplots_adjust(right=0.9)
            cax = ax[i, j].inset_axes(
                [1.05, 0.1, 0.05, 0.8], transform=ax[i, j].transAxes
            )
            fig.colorbar(plot, ax=ax[i, j], cax=cax, label=labelplot)

fig.savefig("Figure1.png", bbox_inches="tight", dpi=400)
plt.close(fig)
