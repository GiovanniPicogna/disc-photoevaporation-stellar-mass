import os
import numpy as np
import h5py as h
from astropy import constants as const
from astropy import units as u
import matplotlib.pyplot as plt
from matplotlib import rc
from functions import moving_average

plt.style.use('figures.mplstyle')

plt.rcParams.update({'figure.figsize': (7, 9)})

fig, ax = plt.subplots(2, 1, sharey=False, sharex=True)
fig.subplots_adjust(left=0.15, bottom=0.16, right=0.99, top=0.97)
dirs = [
    "data/PLUTO/01Msun",
    "data/PLUTO/03Msun",
    "data/PLUTO/05Msun",
    "data/PLUTO/1Msun",
]
Mstar = [0.1, 0.3, 0.5, 1.0]
labels = [r"\texttt{0.1Msun}", r"\texttt{0.3Msun}", r"\texttt{0.5Msun}", r"\texttt{1Msun}"]
colors = ["blue", "red", "orange", "green"]
mu = 1.4

window = 150
os.chdir("../")
for j in range(2):
    for s in range(4):
        os.chdir(dirs[s])
        Rg = 5.05 * Mstar[s]
        sonic = np.loadtxt("sonicsurf.dat")
        radius_ss = moving_average(sonic[:, 0], window)
        temperature_ss = moving_average(sonic[:, 1], window)
        density_ss = moving_average(sonic[:, 2], window)

        cs2_theo = (const.G * const.M_sun * Mstar[s]) / (
            2.0 * (radius_ss * u.AU).to(u.cm)
        )
        T_theo = ((mu * const.m_p * cs2_theo) / (const.k_B)).to(u.K)
        if j == 0:
            ax[j].loglog(
                radius_ss / Rg, density_ss, ".", color=colors[s], label=labels[s]
            )
        else:
            ax[j].loglog(
                radius_ss / Rg, temperature_ss, ".", color=colors[s], label=labels[s]
            )

        if s == 3 and j == 1:
            ax[j].loglog(radius_ss / Rg, T_theo, "k--", label="theo")

        if j == 0:
            labelplot = r"$\rho$ [g cm$^{-3}$]"
        else:
            ax[j].set_xlabel("R [R$_g$]")
            labelplot = r"T [K]"
        ax[j].set_ylabel(labelplot)
        ax[j].set_xlim(0.1, 20)
        if j == 0:
            ax[j].set_ylim(5.0e-20, 5.0e-14)
        else:
            ax[j].set_ylim(1e3, 1e4)
            ax[j].legend()
        os.chdir("../../../")
plt.tight_layout(pad=0.0)
plt.savefig("script/Figure9.pdf", bbox_inches="tight", dpi=400)