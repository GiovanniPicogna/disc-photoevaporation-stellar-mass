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


def Sigmadot(a, b, c, d, e, f, g, Rau):
    logR = np.log10(Rau)
    lnx = np.log(Rau)
    ln10 = np.log(10)
    return (
        ln10
        * (
            6 * a * lnx ** 5 / (Rau * ln10 ** 6)
            + 5 * b * lnx ** 4 / (Rau * ln10 ** 5)
            + 4 * c * lnx ** 3 / (Rau * ln10 ** 4)
            + 3 * d * lnx ** 2 / (Rau * ln10 ** 3)
            + 2 * e * lnx / (Rau * ln10 ** 2)
            + f / (Rau * ln10)
        )
        * 10
        ** (
            a * logR ** 6
            + b * logR ** 5
            + c * logR ** 4
            + d * logR ** 3
            + e * logR ** 2
            + f * logR
            + g
        )
    )


Mstar = [0.1, 0.3, 0.5, 1.0]

width = 7
height = 4.5
fig, ax = plt.subplots(ncols=1, nrows=1)
labels = [r"\texttt{01Msun}", r"\texttt{03Msun}", r"\texttt{05Msun}", r"\texttt{1Msun}"]
Mdot = [
    pow(10, -8.448195192063903),
    pow(10, -7.918753594822812),
    pow(10, -7.7478864615527225),
    pow(10, -7.58791478942398),
    pow(10, -7.424759928546305),
]
idx = 1.0
mu = 2.35
x = np.linspace(0.1, 50, 4000)
a01 = -3.83368381
b01 = 22.9099733
c01 = -55.1282068
d01 = 67.8918937
e01 = -45.0137998
f01 = 16.2976639
g01 = -3.54258614
Mdot01 = 3.55e-09
Ts = (
    (0.125 * x ** (0.85) * mu * const.m_p * const.G * const.M_sun * Mstar[0])
    / (2.0 * const.k_B * ((x * u.AU).to(u.cm)))
).to(u.K)
Rg = 5.05 * 0.1
plt.loglog(
    x / Rg,
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
a03 = -1.32064160
b03 = 13.0474633
c03 = -53.6989939
d03 = 117.602662
e03 = -144.376859
f03 = 94.7854366
g03 = -26.7363296
Mdot03 = 9.63e-9
Ts = (
    (0.125 * x ** (0.85) * mu * const.m_p * const.G * const.M_sun * Mstar[1])
    / (2.0 * const.k_B * ((x * u.AU).to(u.cm)))
).to(u.K)
Rg = 5.05 * 0.3
plt.loglog(
    x / Rg,
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
b05 = 10.8504904
c05 = -38.6938600
d05 = 71.2489311
e05 = -71.4279364
f05 = 37.8706679
g05 = -9.35078265
Mdot05 = 1.7e-08
Ts = (
    (0.125 * x ** (0.85) * mu * const.m_p * const.G * const.M_sun * Mstar[2])
    / (2.0 * const.k_B * ((x * u.AU).to(u.cm)))
).to(u.K)
Rg = 5.05 * 0.5
plt.loglog(
    x / Rg,
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
a1 = -0.63444166
b1 = 6.35869647
c1 = -26.1444624
d1 = 56.4476652
e1 = -67.7402956
f1 = 43.9212155
g1 = -13.2315666
Mdot1 = 4.35e-08
Ts = (
    (0.125 * x ** (0.85) * mu * const.m_p * const.G * const.M_sun * Mstar[3])
    / (2.0 * const.k_B * ((x * u.AU).to(u.cm)))
).to(u.K)
Rg = 5.05
plt.loglog(
    x / Rg,
    (
        (
            Sigmadot(a1, b1, c1, d1, e1, f1, g1, x)
            * (Mdot[4] * u.M_sun / u.yr / u.au).to(u.g / u.s / u.au)
            / (2.0 * np.pi * x * u.au)
        ).to(u.g / u.cm / u.cm / u.s)
    ),
    "k-",
    label=labels[3],
)

plt.xlabel(r"R [R$_g$]", fontsize=14)
plt.ylabel(r"Culmative Mass-loss Rate", fontsize=14)
plt.ylim(1.0e-15, 2.0e-11)
plt.xlim(0.2, 100)
plt.legend(loc=3)

plt.ylabel(r"$\dot{\Sigma}_w$ [g cm$^{-2}$ s$^{-1}$]")
plt.legend(loc=1)

plt.tight_layout(pad=0.0)
fig.set_size_inches(width, height)

plt.savefig("Figure7b.pdf", bbox_inches="tight", dpi=400)
plt.close(fig)
