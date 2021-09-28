import os
import numpy as np
import h5py as h
from astropy import constants as const
from astropy import units as u
import matplotlib.pyplot as plt
from matplotlib import rc
from functions import sci_notation

plt.style.use('figures.mplstyle')

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
ax1.set_ylabel("$\log_{10}(\dot{\mathrm{M}}_w$ [$M_\odot$/yr])", size=14)
ax1.set_xlabel("$M_\star$ [$M_\odot$]", size=14)
ax1.set_xticks([0.1, 0.3, 0.5, 1.0])
ax1.set_xticklabels([0.1, 0.3, 0.5, 1.0])
ax1.set_xlim(0.02, 1.1)
ax1.legend(loc=2, fontsize=14)
yticks = [5e-9, 1e-8, 2e-8]
ax1.set_yticks(yticks)
ax1.set_yticklabels([sci_notation(y, decimal_digits=0) for y in yticks])

plt.tight_layout(pad=0.0)
plt.savefig("Figure8.pdf", bbox_inches="tight", dpi=400)
