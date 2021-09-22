import numpy as np
import matplotlib as mpl

mpl.use("pdf")
import matplotlib.pyplot as plt
import os
from astropy import constants as const
from astropy import units as u

import h5py as h

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


def getFilenames():
    print(os.getcwd())
    filenames = [x for x in os.listdir("./") if ".h5" in x]
    filenames.sort()
    return filenames


def getVar(filename, step, variable):
    h5 = h.File(filename, "r")
    returnData = h5["Timestep_" + str(step) + "/vars"][variable][:]
    h5.close()
    return returnData


def getGridCell(filename=None, all=1):
    if not (filename):
        filename = getFilenames()[0]
    h5 = h.File(filename, "r")
    if all:
        x = h5["cell_coords"]["X"][:]
        y = h5["cell_coords"]["Y"][:]
        z = h5["cell_coords"]["Z"][:]
    else:
        x = h5["cell_coords"]["X"]
        y = h5["cell_coords"]["Y"]
        z = h5["cell_coords"]["Z"]
    x = x.astype("float64")
    y = y.astype("float64")
    z = z.astype("float64")
    return x, y, z


rscale = 10.0 * u.AU
mscale = 1.0 * u.solMass
vscale = np.sqrt(const.G * mscale / rscale) / 2.0 / np.pi
rhoscale = mscale / rscale ** 3
mu = 2.35

width = 12
height = 6
fig = plt.figure()
gs = fig.add_gridspec(nrows=2, ncols=4, hspace=0, wspace=0)
ax = gs.subplots(sharex=True, sharey=True)

dirs = [
    "01Msun",
    "03Msun",
    "05Msun",
    "1Msun",
]
labels = ["\\texttt{01Msun}", "\\texttt{03Msun}", "\\texttt{05Msun}", "\\texttt{1Msun}"]
Mstar = [0.1, 0.3, 0.5, 1.0]
index_max = [415, 470, 452, 553]

xticks1 = np.arange(0, 21, 5)
xticks = np.arange(0, 21, 5)
os.chdir("../data/PLUTO/")
for i in range(2):
    for j in range(4):
        os.chdir(dirs[j])
        ax[i, j].set_xlim(0, 24.99)
        ax[i, j].set_ylim(0, 24.99)
        ax[i, j].set_aspect("equal")
        if i == 1:
            ax[i, j].set_xlabel(r"R [R$_g$]", size=14, color="black")
        if j == 0:
            ax[i, j].set_ylabel(r"Z [R$_g$]", size=14, color="black")

        files = getFilenames()[:]
        Rg = 5.05 * Mstar[j]
        step = 0
        xcell, ycell, zcell = (getGridCell() * rscale).to(u.cm)
        X = xcell
        Z = ycell
        D = (getVar(files[step], step, "rho") * rhoscale).to(u.g / u.cm ** 3)
        Pr = (getVar(files[step], step, "prs") * vscale ** 2 * rhoscale).to(u.barye)
        T = ((Pr * mu * const.m_p) / (const.k_B * D)).to(u.K).value
        vr = (getVar(files[step], step, "vx1") * vscale).to(u.cm / u.s)
        vth = (getVar(files[step], step, "vx2") * vscale).to(u.cm / u.s)
        cs2 = (const.G * Mstar[j] * const.M_sun / (2 * X)).to(u.cm ** 2 / u.s ** 2)
        vel2 = vr ** 2 + vth ** 2
        Rstart = []
        if i == 0:
            quantity = T
            labelplot = r"T [K]"
            valmin = 1.0e3
            valmax = 1.1e4
        else:
            quantity = np.log10(D.value / 10 ** (-24))
            labelplot = r"$\log_{10}(\rho)$ [$10^{-24}$ g cm$^{-3}$]"
            valmin = 3
            valmax = 11
        Cd = getVar(files[step], step, "cd")
        if i == 0:
            ax[i, j].set_title(labels[j], size=14.0)

        plot = ax[i, j].pcolormesh(
            (X).to(u.AU).value / Rg,
            (Z).to(u.AU).value / Rg,
            quantity,
            shading="auto",
            vmin=valmin,
            vmax=valmax,
            cmap=plt.cm.jet,
        )
        sonic = np.loadtxt("sonic1.dat")
        soundx = (np.trim_zeros(sonic[:, 0]) * u.cm).to(u.AU)
        soundz = (np.trim_zeros(sonic[:, 1]) * u.cm).to(u.AU)
        ax[i, j].plot(soundx / Rg, soundz / Rg, "r-")
        for pr in range(20):
            strlines = np.loadtxt("stream" + str(pr + 1) + ".dat")
            strx = ((strlines[:, 0]) * u.cm).to(u.AU)
            strz = ((strlines[:, 1]) * u.cm).to(u.AU)
            ax[i, j].plot(strx / Rg, strz / Rg, ls="--", color="lightcoral")

        vals = [2.0e22]
        contours = ax[i, j].contour(
            X.to(u.AU).value / Rg,
            Z.to(u.AU).value / Rg,
            Cd,
            vals,
            linestyles="dashed",
            colors="red",
        )
        ax[i, j].set_xticks(xticks)

        if j == 3:
            fig.subplots_adjust(right=0.9)
            cax = ax[i, j].inset_axes(
                [1.05, 0.1, 0.05, 0.8], transform=ax[i, j].transAxes
            )
            fig.colorbar(plot, ax=ax[i, j], cax=cax, label=labelplot)
        os.chdir("../")

plt.tight_layout(pad=0.0)
fig.set_size_inches(width, height)

os.chdir("../../script")
fig.savefig("Figure4.png", bbox_inches="tight", dpi=400)
plt.close(fig)
