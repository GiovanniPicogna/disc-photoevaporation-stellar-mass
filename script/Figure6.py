import os
import numpy as np
import h5py as h
from astropy import constants as const
from astropy import units as u
import matplotlib.pyplot as plt
from matplotlib import rc
from functions import MdotLxPicogna,Lx,sci_notation

plt.style.use('figures.mplstyle')

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
fin_idx = -50
data = [
    (mdot01[fin_idx:]),
    (mdot03[fin_idx:]),
    (mdot05[fin_idx:]),
    (mdot07[fin_idx:]),
    (mdot10[fin_idx:]),
]

fig1, ax1 = plt.subplots()
ax1.boxplot(data, positions=[0.1, 0.3, 0.5, 0.7, 1.0])
x = np.linspace(0.1, 1.0, 1000)
ax1.plot(x, 3.93e-8 * x, "k-.", label="$3.93 \cdot 10^{-8} \, \mathrm{M}_\star$")
ax1.fill_between(
    x,
    3.93e-8 * x * MdotLxPicogna(Lx(x)[0], Lx(x)[1]),
    3.93e-8 * x * MdotLxPicogna(Lx(x)[2], Lx(x)[1]),
    color="lightblue",
    alpha=0.6,
)

ax1.set_ylabel("$\log_{10}(\dot{\mathrm{M}}_w$ [$M_\odot$/yr])", size=14)
ax1.set_xlabel("$M_\star$ [$M_\odot$]", size=14)
ax1.set_xticks([0.1, 0.3, 0.5, 0.7, 1.0])
ax1.set_xticklabels([0.1, 0.3, 0.5, 0.7, 1.0])
ax1.set_xlim(0.02, 1.1)
ax1.set_ylim(1.0e-10, 4.5e-8)
ax1.legend(loc=2, fontsize=14)
yticks = [4e-9, 1.2e-8, 2.0e-8, 2.8e-8, 3.6e-8, 4.4e-8]
ax1.set_yticks(yticks)
ax1.set_yticklabels([sci_notation(y, decimal_digits=1) for y in yticks])
plt.tight_layout(pad=0.0)
plt.savefig("Figure6.pdf", bbox_inches="tight", dpi=400)
