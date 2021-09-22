import os
import numpy as np
import h5py as h
from astropy import constants as const
from astropy import units as u
import matplotlib as mpl

mpl.use("pdf")
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab

params = {
    "legend.fontsize": "x-large",
    "figure.figsize": (7, 4.5),
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


def MdotLxPicogna(Lx1, Lx2):
    AL = -2.7326
    BL = 3.3307
    CL = -2.9868e-3
    DL = -7.2580
    mdot1 = 10 ** (AL * np.exp(((np.log(np.log10(Lx1)) - BL) ** 2) / CL) + DL)
    mdot2 = 10 ** (AL * np.exp(((np.log(np.log10(Lx2)) - BL) ** 2) / CL) + DL)
    return mdot1 / mdot2


def Lx(*Mstar):
    for x in Mstar:
        maxLx = 10 ** (1.42 * np.log10(x) + 30.37)
        minLx = 10 ** (1.66 * np.log10(x) + 30.25)
        aveLx = 10 ** (1.54 * np.log10(x) + 30.31)
        return minLx, aveLx, maxLx


def sci_notation(num, decimal_digits=1, precision=None, exponent=None):
    """
    Returns a string representation of the scientific
    notation of the given number formatted for use with
    LaTeX or Mathtext, with specified number of significant
    decimal digits and precision (number of decimal digits
    to show). The exponent to be used can also be specified
    explicitly.
    """
    if exponent is None:
        exponent = int(np.floor(np.log10(abs(num))))
    coeff = round(num / float(10 ** exponent), decimal_digits)
    if precision is None:
        precision = decimal_digits

    return r"${0:.{2}f}\cdot10^{{{1:d}}}$".format(coeff, exponent, precision)


mdot01 = np.loadtxt("../data/PLUTO/01Msun/mdot_ave.dat")
mdot03 = np.loadtxt("../data/PLUTO/03Msun/mdot_ave.dat")
mdot05 = np.loadtxt("../data/PLUTO/05Msun/mdot_ave.dat")
mdot07 = np.loadtxt("../data/PLUTO/07Msun/mdot_ave.dat")
mdot10 = np.loadtxt("../data/PLUTO/1Msun/mdot_ave.dat")
Lx07 = 10 ** (1.54 * np.log10(0.7) + 30.31)

AL = -2.7326
BL = 3.3307
CL = -2.9868e-3
DL = -7.2580
frac = 10 ** (AL * np.exp(((np.log(np.log10(2e30)) - BL) ** 2.0) / CL) + DL) / 10 ** (
    AL * np.exp(((np.log(np.log10(Lx07)) - BL) ** 2.0) / CL) + DL
)
mdot07 = mdot07 / frac
init_idx = 0
fin_idx = -50
data = [
    (mdot01[fin_idx:]),
    (mdot03[fin_idx:]),
    (mdot05[fin_idx:]),
    (mdot07[fin_idx:]),
    (mdot10[fin_idx:]),
]

width = 7
height = 4.5
fig1, ax1 = plt.subplots()
ax1.boxplot(data, positions=[0.1, 0.3, 0.5, 0.7, 1.0])
x = np.linspace(0.1, 1.0, 1000)
Mstelle = np.asarray([0.1, 0.3, 0.5, 0.7, 1.0])
ax1.plot(x, 3.93e-8 * x, "k-.", label="$3.93 \cdot 10^{-8} \, \mathrm{M}_\star$")
ax1.fill_between(
    x,
    3.93e-8 * x * MdotLxPicogna(Lx(x)[0], Lx(x)[1]),
    3.93e-8 * x * MdotLxPicogna(Lx(x)[2], Lx(x)[1]),
    color="lightblue",
    alpha=0.6,
)

ax1.set_ylabel("$\log_{10}(\dot{M}$ [$M_\odot$/yr])", size=14)
ax1.set_xlabel("$M_\star$ [$M_\odot$]", size=14)
ax1.set_xticks([0.1, 0.3, 0.5, 0.7, 1.0])
ax1.set_xticklabels([0.1, 0.3, 0.5, 0.7, 1.0])
ax1.set_xlim(0.02, 1.1)
ax1.set_ylim(1.0e-10, 4.5e-8)
ax1.legend(loc=2, fontsize=14)
yticks = [4e-9, 1.2e-8, 2.0e-8, 2.8e-8, 3.6e-8, 4.4e-8]
ax1.set_yticks(yticks)
ax1.set_yticklabels([sci_notation(y, decimal_digits=1) for y in yticks])
# use LaTeX formatted labels
plt.tight_layout(pad=0.0)
fig1.set_size_inches(width, height)
plt.savefig("Figure6.pdf", bbox_inches="tight", dpi=400)
