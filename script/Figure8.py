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


def MdotOwen12(Mstar):
    Lx = 10 ** (1.54 * np.log10(Mstar) + 30.31)
    return 6.25e-9 * (Mstar ** (-0.068)) * ((Lx / 1.0e30) ** 1.14)


mdot01 = np.loadtxt("../data/PLUTO/test_mass/01Msun/mdot_ave.dat")
mdot03 = np.loadtxt("../data/PLUTO/test_mass/03Msun/mdot_ave.dat")
mdot05 = np.loadtxt("../data/PLUTO/test_mass/05Msun/mdot_ave.dat")
mdot10 = np.loadtxt("../data/PLUTO/test_mass/1Msun/mdot_ave.dat")

fin_idx = -50
data = [
    (mdot01[fin_idx:]),
    (mdot03[fin_idx:]),
    (mdot05[fin_idx:]),
    (mdot10[fin_idx:]),
]

width = 7
height = 4.5
fig1, ax1 = plt.subplots()
ax1.boxplot(data, positions=[0.1, 0.3, 0.5, 1.0])
x = np.linspace(0.1, 1.0, 1000)

ax1.plot(
    x,
    2.285123e-8 * x + 6.824824e-9,
    "k-.",
    label="$2.285\cdot 10^{-8} M_\star + 6.825\cdot 10^{-9}$",
)
ax1.plot(
    x,
    6.25e-9 * x ** (-0.068) * (7.02e29 / 1.0e30) ** (1.14),
    "k:",
    label="Owen et al., 2012",
)
ax1.fill_between(
    x,
    (2.285e-8 * x + 6.825e-9) * MdotLxPicogna(Lx(x)[0], Lx(x)[1]),
    (2.285e-8 * x + 6.825e-9) * MdotLxPicogna(Lx(x)[2], Lx(x)[1]),
    color="lightblue",
    alpha=0.6,
)
ax1.set_ylabel("$\log_{10}(\dot{\mathrm{M}}_w$ [$M_\odot$/yr])", size=14)
ax1.set_xlabel("$M_\star$ [$M_\odot$]", size=14)
ax1.set_xticks([0.1, 0.3, 0.5, 1.0])
ax1.set_xticklabels([0.1, 0.3, 0.5, 1.0])
ax1.set_xlim(0.02, 1.1)
ax1.legend(loc=2, fontsize=14)
yticks = [5e-9, 1e-8, 2e-8]
ax1.set_yticks(yticks)
ax1.set_yticklabels([sci_notation(y, decimal_digits=0) for y in yticks])
# use LaTeX formatted labels
plt.tight_layout(pad=0.0)
fig1.set_size_inches(width, height)
plt.savefig("Figure8.pdf", bbox_inches="tight", dpi=400)
