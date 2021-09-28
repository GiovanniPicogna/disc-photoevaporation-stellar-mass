import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc

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

data = [mdot01, mdot03, mdot05, mdot07, mdot10]

fig = plt.figure()
labels = [
    r"\texttt{0.1Msun}",
    r"\texttt{0.3Msun}",
    r"\texttt{0.5Msun}",
    r"\texttt{0.7Msun}",
    r"\texttt{1Msun}",
]
colors = ["black", "black", "black", "darkgrey", "black"]
linestyles = ["-", ":", "-.", "-", "--"]
fig = plt.figure()
for i in range(5):
    plt.plot(data[i], color=colors[i], linestyle=linestyles[i], label=labels[i])
plt.ylabel(r"$\log_{10}(\mathrm{\dot{M}}_w\,\, [M_\odot\, \textrm{yr}^{-1}])$")
plt.xlabel(r"t [orbits @ 10 au]")
plt.yscale("log")
plt.legend(ncol=2, loc="lower right")
plt.xlim(0, 300)

plt.savefig("Figure3.pdf", bbox_inches="tight", dpi=400)
plt.close(fig)
