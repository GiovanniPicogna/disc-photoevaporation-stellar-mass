import os
import numpy as np
from astropy import units as u
from astropy import constants as const
from scipy.interpolate import griddata
import matplotlib.pyplot as plt

import matplotlib as mpl

mpl.use("pdf")

import matplotlib.pylab as pylab

params = {
    "legend.fontsize": "x-large",
    "figure.figsize": (12, 6),
    "axes.labelsize": "x-large",
    "axes.titlesize": "x-large",
    "xtick.labelsize": "x-large",
    "ytick.labelsize": "x-large",
}
pylab.rcParams.update(params)
from matplotlib import rc

rc("font", **{"family": "sans-serif"})
rc("text", usetex=True)
rc("xtick", labelsize=14.0)
rc("ytick", labelsize=14.0)
rc("axes", labelsize=14.0)

basepath = "../data/DIAD/"
files = [
    "fort14.a0.01.irr.grid_001.e1.0.amax0p25_amax1mm.dat",
    "fort14.a0.01.irr.grid_002.e1.0.amax0p25_amax1mm.dat",
    "fort14.a0.01.irr.grid_003.e1.0.amax0p25_amax1mm.dat",
    "fort14.a0.01.irr.grid_004.e1.0.amax0p25_amax1mm.dat",
]
Rstars = [1.055, 2.31, 2.125, 2.615]
Mstars = [0.1, 0.3, 0.5, 1.0]

width = 12
height = 6
fig = plt.figure()
gs = fig.add_gridspec(nrows=2, ncols=4, hspace=0, wspace=0)
ax = gs.subplots(sharex=True, sharey=True)

base_path = "/e/arch/users/picogna/Photoevaporation/Mass/"
dirs = ["Final/01Msun", "Final/03Msun", "Final/05Msun", "Final/1Msun"]
labels = ["\\texttt{01Msun}", "\\texttt{03Msun}", "\\texttt{05Msun}", "\\texttt{1Msun}"]
Mstar = [0.1, 0.3, 0.5, 1.0]
index_max = [415, 470, 452, 553]

xticks1 = np.arange(0, 21, 5)
xticks = np.arange(0, 21, 5)
print(xticks)

for i in range(2):
    for j in range(4):
        ax[i, j].set_xlim(0, 24)
        ax[i, j].set_ylim(0, 24)
        ax[i, j].set_aspect("equal")
        if i == 1:
            ax[i, j].set_xlabel(r"R [R$_g$]", size=14, color="black")
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

        Rg = 5.05 * Mstars[j]

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

        if i == 0:
            ax[i, j].set_title(labels[j], size=14.0)
            ax[i, j].spines["top"].set_visible(True)

        plot = ax[i, j].pcolormesh(
            Xc * u.cm.to(u.AU) / Rg,
            Zc * u.cm.to(u.AU) / Rg,
            quantity,
            shading="auto",
            vmin=valmin,
            vmax=valmax,
            cmap=plt.cm.jet,
        )
        ax[i, j].set_xticks(xticks)

        if j == 3:
            fig.subplots_adjust(right=0.9)
            cax = ax[i, j].inset_axes(
                [1.05, 0.1, 0.05, 0.8], transform=ax[i, j].transAxes
            )
            fig.colorbar(plot, ax=ax[i, j], cax=cax, label=labelplot)

plt.tight_layout(pad=0.0)
fig.set_size_inches(width, height)
fig.savefig("Figure1.png", bbox_inches="tight", dpi=400)
plt.close(fig)
