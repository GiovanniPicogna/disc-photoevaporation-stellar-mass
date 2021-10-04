import os
import numpy as np
import h5py as h
from astropy import constants as const
from astropy import units as u
import matplotlib.pyplot as plt
from matplotlib import rc
from functions import Sigmadot

plt.style.use('figures.mplstyle')

Mstar = [0.1, 0.3, 0.5, 1.0]

fig, ax = plt.subplots(ncols=1, nrows=1)

labels = [r"\texttt{0.1Msun}", r"\texttt{0.3Msun}", r"\texttt{0.5Msun}", r"\texttt{1Msun}"]
Mdot = [8.746e-09, 1.349e-08, 1.918e-08, 2.931e-08]

x = np.linspace(0.1, 100, 4000)
a01 = -0.93483826 
b01 =  5.56625743
c01 = -13.1970819
d01 =  16.1436840
e01 = -11.3750655
f01 =  5.55895808
g01 = -2.13751254
plt.loglog(
    x,
    (
        (
            Sigmadot(a01, b01, c01, d01, e01, f01, g01, x)
            * (Mdot[0] * u.M_sun / u.yr / u.au).to(u.g / u.s / u.au)
            / (2.0 * np.pi * x * u.au)
        ).to(u.g / u.cm / u.cm / u.s)
    ),
    "k:",
    label=labels[0],
)

x = np.linspace(0.1, 180, 4000)
a03 = -0.54742481 
b03 =  4.73555351 
c03 = -16.4134522 
d03 =  28.9407822 
e03 = -27.4476651 
f03 =  14.1897770 
g03 = -4.01282828
plt.loglog(
    x,
    (
        (
            Sigmadot(a03, b03, c03, d03, e03, f03, g03, x)
            * (Mdot[1] * u.M_sun / u.yr / u.au).to(u.g / u.s / u.au)
            / (2.0 * np.pi * x * u.au)
        ).to(u.g / u.cm / u.cm / u.s)
    ),
    "k-.",
    label=labels[1],
)

x = np.linspace(0.1, 600, 4000)
a05 = -1.23202374
b05 =  10.8504904
c05 = -38.6938600
d05 =  71.2489311
e05 = -71.4279364
f05 =  37.8706679
g05 = -9.35078265
plt.loglog(
    x,
    (
        (
            Sigmadot(a05, b05, c05, d05, e05, f05, g05, x)
            * (Mdot[2] * u.M_sun / u.yr / u.au).to(u.g / u.s / u.au)
            / (2.0 * np.pi * x * u.au)
        ).to(u.g / u.cm / u.cm / u.s)
    ),
    "k--",
    label=labels[2],
)

x = np.linspace(0.1, 600, 4000)
a1 = -0.46180984
b1 =  4.78877278
c1 = -20.3882637
d1 =  45.5660721
e1 = -56.5313405
f1 =  37.9238589
g1 = -11.9467093
plt.loglog(
    x,
    (
        (
            Sigmadot(a1, b1, c1, d1, e1, f1, g1, x)
            * (Mdot[3] * u.M_sun / u.yr / u.au).to(u.g / u.s / u.au)
            / (2.0 * np.pi * x * u.au)
        ).to(u.g / u.cm / u.cm / u.s)
    ),
    "k-",
    label=labels[3],
)

plt.xlabel(r"R [au]", fontsize=14)
plt.ylabel(r"Culmative Mass-loss Rate", fontsize=14)
plt.ylim(1.0e-15, 2.0e-11)
plt.xlim(0.3, 350)
plt.legend(loc=3)

plt.ylabel(r"$\dot{\Sigma}_w$ [g cm$^{-2}$ s$^{-1}$]")
plt.legend(loc=1)

plt.tight_layout(pad=0.0)

plt.savefig("Figure12.pdf", bbox_inches="tight", dpi=400)
plt.close(fig)
