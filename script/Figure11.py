import os
import numpy as np
import h5py as h
from astropy import constants as const
from astropy import units as u
import pandas as pd
import matplotlib.pyplot as plt

plt.style.use('figures.mplstyle')

data = {
    "Stellar mass": [
        0.1,
        0.15,
        0.2,
        0.25,
        0.3,
        0.35,
        0.5,
        0.55,
        0.7,
        0.75,
        0.8,
        0.85,
        0.9,
        0.95,
        1,
    ],
    "Disk fraction": [
        0.77,
        0.54,
        0.44,
        0.41,
        0.35,
        0.34,
        0.34,
        0.32,
        0.32,
        0.31,
        0.30,
        0.30,
        0.29,
        0.28,
        0.27,
    ],
}

df = pd.DataFrame(data, columns=["Stellar mass", "Disk fraction"])


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


def life_time_Alcala(Mstar, LX, alpha1, beta1, alpha2, beta2, age, profile):
    # Mstar in Solar masses
    # age in yr
    if profile == "Picogna":
        Mdisc0 = 0.14
        Mdisc0 = Mdisc0 * Mstar ** 1.0
    elif profile == "Owen":
        Mdisc0 = 0.06
        Mdisc0 = Mdisc0 * Mstar ** 1.0
    else:
        print("profile not implemented")
        Mdisc0 = 0

    Mdot = np.where(
        Mstar >= 0.2,
        10 ** (alpha2 * np.log10(Mstar) + beta2),
        10 ** (alpha1 * np.log10(Mstar) + beta1),
    )
    Mdot0 = (0.5 * Mdisc0 / age) ** 3 / Mdot ** 2
    if profile == "Picogna":
        Mwind = 3.93e-8 * Mstar * MdotLxPicogna(LX, Lx(Mstar)[1])
    elif profile == "Owen":
        Mwind = 6.25e-9 * Mstar ** (-0.068) * (LX / 1.0e30) ** (1.14)
    else:
        print("profile not implemented")
        Mwind = 0
    tlife = 0.5 * Mdisc0 * pow(Mdot, -1 / 3) * pow(Mwind, -2 / 3)
    return tlife


Mstars = np.asarray([0.1, 0.3, 0.5, 0.7, 1.0])
LX = np.asarray([0.059e30, 0.32e30, 0.702e30, 1.179e30, 2.04e30])
Mstars_arr = np.asarray(np.linspace(0.1, 1.0, 10))
LX_arr = 10 ** (1.54 * np.log10(Mstars_arr) + 30.31)
idx_disk = 0.99

alpha1 = [3.9, 5.26]
alpha1mean = 4.58
alpha2 = [1.13, 1.61]
alpha2mean = 1.37

beta1 = [-5.5, -6.72]
beta1mean = -6.11
beta2 = [-8.35, -8.57]
beta2mean = -8.46

life = np.zeros((5, 2))
for i in range(5):
    for j in range(2):
        life[i, j] = life_time_Alcala(
            Mstars[i], LX[i], alpha1[j], beta1[j], alpha2[j], beta2[j], 1.5e6, "Owen"
        )
lifeTau = np.zeros(5)
for i in range(5):
    lifeTau[i] = life_time_Alcala(
        Mstars[i], LX[i], alpha1mean, beta1mean, alpha2mean, beta2mean, 1.5e6, "Owen"
    )
life1 = np.zeros((5, 2))
for i in range(5):
    for j in range(2):
        life1[i, j] = life_time_Alcala(
            Mstars[i], LX[i], alpha1[j], beta1[j], alpha2[j], beta2[j], 1.5e6, "Picogna"
        )
lifeLup = np.zeros(5)
for i in range(5):
    lifeLup[i] = life_time_Alcala(
        Mstars[i], LX[i], alpha1mean, beta1mean, alpha2mean, beta2mean, 1.5e6, "Picogna"
    )

fig1, ax1 = plt.subplots()

Mstars = np.asarray([0.1, 0.3, 0.5, 0.7, 1.0])
LX = np.asarray([0.059, 0.32, 0.702, 1.179, 2.04])
Mstars_arr = np.asarray(np.linspace(0.1, 1.0, 37))

LXmin_arr = 10 ** (1.42 * np.log10(Mstars_arr) + 30.37)
LXmax_arr = 10 ** (1.66 * np.log10(Mstars_arr) + 30.25)
LX_arr = 10 ** (1.54 * np.log10(Mstars_arr) + 30.31)

colors = ["lightcoral", "firebrick", "red"]
plt.plot(
    data["Stellar mass"],
    -(5e6 / np.log(data["Disk fraction"])),
    ".",
    ms=20,
    color="grey",
    label="$\lambda$ Ori [5 Myr, Bayo+ 2012]",
)

labels = [r"$\alpha = 0.99$", r"$\alpha = 1.13$", r"$\alpha = 1.10$"]

plt.fill_between(
    Mstars_arr,
    life_time_Alcala(
        Mstars_arr,
        LXmin_arr,
        alpha1[0],
        beta1[0],
        alpha2[0],
        beta2[0],
        1.5e6,
        "Picogna",
    ),
    life_time_Alcala(
        Mstars_arr,
        LXmax_arr,
        alpha1[1],
        beta1[1],
        alpha2[1],
        beta2[1],
        1.5e6,
        "Picogna",
    ),
    color="lightseagreen",
    label="this work, M$_\mathrm{disc,0}$ = 0.14 M$_\star$",
    alpha=0.6,
)
plt.plot(
    Mstars_arr,
    life_time_Alcala(
        Mstars_arr,
        LX_arr,
        alpha1mean,
        beta1mean,
        alpha2mean,
        beta2mean,
        1.5e6,
        "Picogna",
    ),
    color="seagreen",
    ls="--",
)

plt.fill_between(
    Mstars_arr,
    life_time_Alcala(
        Mstars_arr, LXmin_arr, alpha1[0], beta1[0], alpha2[0], beta2[0], 1.5e6, "Owen"
    ),
    life_time_Alcala(
        Mstars_arr, LXmax_arr, alpha1[1], beta1[1], alpha2[1], beta2[1], 1.5e6, "Owen"
    ),
    color="lightgreen",
    label="Owen et al. (2012), M$_\mathrm{disc,0}$ = 0.06 M$_\star$",
    alpha=0.6,
)
plt.plot(
    Mstars_arr,
    life_time_Alcala(
        Mstars_arr, LX_arr, alpha1mean, beta1mean, alpha2mean, beta2mean, 1.5e6, "Owen"
    ),
    color="seagreen",
    ls="--",
)

plt.ylabel("inner disc lifetime [yr]")
plt.xlabel("M$_\star$ [M$_\odot$]")
plt.legend()
plt.yscale("log")

plt.tight_layout(pad=0.0)
plt.savefig("Figure11.pdf", bbox_inches="tight", dpi=400)
