import os
import numpy as np
from astropy import constants as const
from astropy import units as u
import matplotlib.pyplot as plt
from matplotlib import rc
from functions import getFilenames, getVar, getGridCell

plt.style.use('figures.mplstyle')

os.chdir("../data/PLUTO/01Msun/")
files = getFilenames()[:]
step = 0
file = int(step)
rscale = 10.0 * u.AU
mscale = 1.0 * u.solMass
vscale = np.sqrt(const.G * mscale / rscale) / 2.0 / np.pi
rhoscale = mscale / rscale ** 3

xcell1, ycell1, zcell1 = (getGridCell(files[file]) * rscale).to(u.cm)

DIMT1 = xcell1.shape[0]
DIMR1 = xcell1.shape[1]

mu = 2.35

D01 = (getVar(files[file], step, "rho") * rhoscale).to(u.g / u.cm ** 3)
Cd01 = getVar(files[file], step, "cd")
vr01 = (getVar(files[file], step, "vx1") * vscale).to(u.cm / u.s).value
vth01 = (getVar(files[file], step, "vx2") * vscale).to(u.cm / u.s).value
vph01 = (getVar(files[file], step, "vx3") * vscale).to(u.cm / u.s).value
Pr01 = (getVar(files[file], step, "prs") * vscale ** 2 * rhoscale).to(u.barye)
Tgas01 = ((Pr01 * mu * const.m_p) / (const.k_B * D01)).to(u.K).value

os.chdir("../1Msun/")
files = getFilenames()[:]

xcell10, ycell10, zcell10 = (getGridCell(files[file]) * rscale).to(u.cm)

DIMT10 = xcell10.shape[0]
DIMR10 = xcell10.shape[1]

D10 = (getVar(files[file], step, "rho") * rhoscale).to(u.g / u.cm ** 3)
Cd10 = getVar(files[file], step, "cd")
vr10 = (getVar(files[file], step, "vx1") * vscale).to(u.cm / u.s).value
vth10 = (getVar(files[file], step, "vx2") * vscale).to(u.cm / u.s).value
vph10 = (getVar(files[file], step, "vx3") * vscale).to(u.cm / u.s).value
Pr10 = (getVar(files[file], step, "prs") * vscale ** 2 * rhoscale).to(u.barye)
Tgas10 = ((Pr10 * mu * const.m_p) / (const.k_B * D10)).to(u.K).value

os.chdir("../03Msun/")
files = getFilenames()[:]

xcell3, ycell3, zcell3 = (getGridCell(files[file]) * rscale).to(u.cm)

DIMT3 = xcell3.shape[0]
DIMR3 = xcell3.shape[1]

D03 = (getVar(files[file], step, "rho") * rhoscale).to(u.g / u.cm ** 3)
Cd03 = getVar(files[file], step, "cd")
vr03 = (getVar(files[file], step, "vx1") * vscale).to(u.cm / u.s).value
vth03 = (getVar(files[file], step, "vx2") * vscale).to(u.cm / u.s).value
vph03 = (getVar(files[file], step, "vx3") * vscale).to(u.cm / u.s).value
Pr03 = (getVar(files[file], step, "prs") * vscale ** 2 * rhoscale).to(u.barye)
Tgas03 = ((Pr03 * mu * const.m_p) / (const.k_B * D03)).to(u.K).value

os.chdir("../05Msun/")
files = getFilenames()[:]

xcell5, ycell5, zcell5 = (getGridCell(files[file]) * rscale).to(u.cm)

DIMT5 = xcell5.shape[0]
DIMR5 = xcell5.shape[1]

D05 = (getVar(files[file], step, "rho") * rhoscale).to(u.g / u.cm ** 3)
Cd05 = getVar(files[file], step, "cd")
vr05 = (getVar(files[file], step, "vx1") * vscale).to(u.cm / u.s).value
vth05 = (getVar(files[file], step, "vx2") * vscale).to(u.cm / u.s).value
vph05 = (getVar(files[file], step, "vx3") * vscale).to(u.cm / u.s).value
Pr05 = (getVar(files[file], step, "prs") * vscale ** 2 * rhoscale).to(u.barye)
Tgas05 = ((Pr05 * mu * const.m_p) / (const.k_B * D05)).to(u.K).value


def radius1d1(tvalue, rmin, rmax):
    rad = xcell1[tvalue, rmin:rmax].to(u.AU)
    return rad


def radius1d5(tvalue, rmin, rmax):
    rad = xcell5[tvalue, rmin:rmax].to(u.AU)
    return rad


def radius1d3(tvalue, rmin, rmax):
    rad = xcell3[tvalue, rmin:rmax].to(u.AU)
    return rad


def radius1d10(tvalue, rmin, rmax):
    rad = xcell10[tvalue, rmin:rmax].to(u.AU)
    return rad


labels = [r"\texttt{0.1Msun}", r"\texttt{0.3Msun}", r"\texttt{0.5Msun}", r"\texttt{1Msun}"]

fig, ax = plt.subplots()
r1 = radius1d1(DIMT1 - 1, 0, DIMR1 - 1)
gamma = 1.4
OmegaKep = np.sqrt(const.G * 0.1 * const.M_sun / r1.to(u.cm) ** 3).to(1.0 / u.s)
H1 = (
    np.sqrt(gamma * Pr01[DIMT1 - 1, 0 : DIMR1 - 1] / D01[DIMT1 - 1, 0 : DIMR1 - 1])
    * 1.0
    / (OmegaKep * (r1.to(u.cm)))
)
ax.plot(r1, H1, "k:", label=labels[0])

r3 = radius1d3(DIMT3 - 1, 0, DIMR3 - 1)
OmegaKep = np.sqrt(const.G * 0.3 * const.M_sun / r3.to(u.cm) ** 3).to(1.0 / u.s)
H3 = (
    np.sqrt(gamma * Pr03[DIMT3 - 1, 0 : DIMR3 - 1] / D03[DIMT3 - 1, 0 : DIMR3 - 1])
    * 1.0
    / (OmegaKep * (r3.to(u.cm)))
)
ax.plot(r3, H3, "k-.", label=labels[1])

r5 = radius1d5(DIMT5 - 1, 0, DIMR5 - 1)
OmegaKep = np.sqrt(const.G * 0.5 * const.M_sun / r5.to(u.cm) ** 3).to(1.0 / u.s)
H5 = (
    np.sqrt(gamma * Pr05[DIMT5 - 1, 0 : DIMR5 - 1] / D05[DIMT5 - 1, 0 : DIMR5 - 1])
    * 1.0
    / (OmegaKep * (r5.to(u.cm)))
)
ax.plot(r5, H5, "k--", label=labels[2])

r10 = radius1d10(DIMT10 - 1, 0, DIMR10 - 1)
OmegaKep = np.sqrt(const.G * const.M_sun / r10.to(u.cm) ** 3).to(1.0 / u.s)
H10 = (
    np.sqrt(gamma * Pr10[DIMT10 - 1, 0 : DIMR10 - 1] / D10[DIMT10 - 1, 0 : DIMR10 - 1])
    * 1.0
    / (OmegaKep * (r10.to(u.cm)))
)
ax.plot(r10, H10, "k-", label=labels[3])

ax.set_xlabel("R [au]", fontsize=14)
ax.set_ylabel("Scale Height H/R", fontsize=14)
ax.set_xlim((0.0, 40))
ax.set_ylim(0.03, 0.25)

plt.legend(prop={"size": 14}, loc="upper left")
os.chdir("../../../script")
plt.savefig("Figure5.pdf", bbox_inches="tight", dpi=400)
plt.close(fig)
