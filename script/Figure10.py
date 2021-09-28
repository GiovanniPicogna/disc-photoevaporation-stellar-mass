import os
import numpy as np
import h5py as h
from astropy import constants as const
from astropy import units as u
import matplotlib.pyplot as plt
from functions import moving_average

plt.style.use('figures.mplstyle')

fig, ax = plt.subplots(1, 1, sharey=True, sharex=True)
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
for s in range(4):
    os.chdir(dirs[s])
    Rg = 5.05 * Mstar[s] * Mstar[s]
    sonic = np.loadtxt("sonicsurf.dat")
    radius_ss = moving_average(sonic[:, 0], window)
    temperature_ss = moving_average(sonic[:, 1], window)
    cs2_theo = (const.G * const.M_sun * Mstar[s]) / (2.0 * (radius_ss * u.AU).to(u.cm))
    T_theo = ((mu * const.m_p * cs2_theo) / (const.k_B)).to(u.K)
    ax.plot(radius_ss / Rg, temperature_ss, ".", color=colors[s], label=labels[s])

    if s == 3:
        ax.plot(radius_ss / Rg, T_theo, "k--", label="theo")
    ax.set_xlabel("R [R$_g$/M$_\star$]")
    ax.set_ylabel("T [K]")
    ax.set_xlim(0.1, 40)
    ax.set_ylim(1000, 10000)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.legend(loc="lower left")
    os.chdir("../../../")

plt.savefig("script/Figure10.pdf", bbox_inches="tight", dpi=400)
