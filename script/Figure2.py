import numpy as np
import matplotlib.pyplot as plt
import os
from lmfit import Model, Parameters
from matplotlib import rc

plt.style.use('figures.mplstyle')

from functions import xiT,interpPicog29,interpPicog30,interpPicog31

Te29 = np.zeros((34, 40))
xi29 = np.zeros((34, 40))
filepath = "../data/MOCASSIN/idloutput_baspec29_highres_moremean.dat"
with open(filepath, "r") as fp:
    for j in range(40):
        line = fp.readline()
        for i in range(34):
            line = fp.readline()
            x, y = line.split(",")
            Te29[i, j] = float(x)
            xi29[i, j] = float(y)

Te30 = np.zeros((34, 40))
xi30 = np.zeros((34, 40))
filepath = "../data/MOCASSIN/idloutput_baspec30_highres_moremean.dat"
with open(filepath, "r") as fp:
    for j in range(40):
        fp.readline()
        for i in range(34):
            line = fp.readline()
            x, y = line.split(",")
            Te30[i, j] = float(x)
            xi30[i, j] = float(y)

Te31 = np.zeros((34, 40))
xi31 = np.zeros((34, 40))
filepath = "../data/MOCASSIN/idloutput_baspec31_highres_moremean.dat"
with open(filepath, "r") as fp:
    for j in range(40):
        fp.readline()
        for i in range(34):
            line = fp.readline()
            x, y = line.split(",")
            Te31[i, j] = float(x)
            xi31[i, j] = float(y)

b = -26.5735278
c = -6.59234502
d = 4.03588632
m = 0.20863786
n = 29

bval29 = np.zeros(40)
cval29 = np.zeros(40)
dval29 = np.zeros(40)
mval29 = np.zeros(40)

for i in range(40):
    params = Parameters()
    if i < 4:
        b = -26
        params.add("b", value=b, min=-160.0, max=-10)
    else:
        b = -26.0
        params.add("b", value=b, min=-100.0, max=-16)
    params.add("c", value=c, min=-8.0, max=-3)
    params.add("d", value=d)  # , min=1.14, max=1.16)
    if i < 10:
        params.add("m", value=m, min=0.0, max=0.5)
    else:
        params.add("m", value=m, min=0.0, max=0.7)

    logXi = xi29[:, i]
    logTe = np.log10(Te29[:, i])

    we = [1.0 for j in range(np.size(logTe))]

    for j in range(np.size(logTe)):
        if logXi[j] > -6.5:
            we[j] = 5.0
        if logXi[j] > -5.5:
            we[j] = 10.0
        if logXi[j] < -5 and i > 18:
            we[j] = 0.25
        if logXi[j] < -4 and i > 35:
            we[j] = 0.25
        if logXi[j] > -4.6 and i > 3:
            we[j] = 20.0
        if logXi[j] > -3 and i > 3:
            we[j] = 50.0

    gmod = Model(interpPicog29)
    result = gmod.fit(logTe, x=logXi, params=params, weights=we)
    x = np.linspace(-8, -2, 100)
    y = interpPicog29(
        x,
        result.best_values["b"],
        result.best_values["c"],
        result.best_values["d"],
        result.best_values["m"],
    )

    bval29[i] = result.best_values["b"]
    cval29[i] = result.best_values["c"]
    dval29[i] = result.best_values["d"]
    mval29[i] = result.best_values["m"]

b = -26.5735278
c = -6.59234502
d = 4.03588632
m = 0.20863786
n = 30

bval30 = np.zeros(40)
cval30 = np.zeros(40)
dval30 = np.zeros(40)
mval30 = np.zeros(40)

for i in range(40):
    params = Parameters()
    params.add("b", value=b, min=-80, max=-16.0)
    params.add("c", value=c)
    params.add("d", value=d)
    if i < 10:
        params.add("m", value=m, min=0.0, max=0.5)
    else:
        params.add("m", value=m, min=0.0, max=0.75)

    logXi = xi30[:, i]
    logTe = np.log10(Te30[:, i])

    we = [1.0 for j in range(np.size(logTe))]

    for j in range(np.size(logTe)):
        if logXi[j] > -6.5:
            we[j] = 5.0
        if logXi[j] > -5.5:
            we[j] = 10.0
        if logXi[j] > -5.1 and i > 3:
            we[j] = 20.0
        if logXi[j] > -4.8 and i > 34:
            we[j] = 30.0
        if logXi[j] > -3 and i > 3:
            we[j] = 40.0
        if logXi[j] < -5.8 and i > 3:
            we[j] = 0.25
        if logXi[j] < -4 and i > 35 and logTe[j] > 3.0:
            we[j] = 0.0
        if logXi[j] < -5 and i > 35 and logTe[j] > 2.0:
            we[j] = 0.0

    gmod = Model(interpPicog30)
    result = gmod.fit(logTe, x=logXi, params=params, weights=we)
    x = np.linspace(-8, -2, 100)
    y = interpPicog30(
        x,
        result.best_values["b"],
        result.best_values["c"],
        result.best_values["d"],
        result.best_values["m"],
    )

    bval30[i] = result.best_values["b"]
    cval30[i] = result.best_values["c"]
    dval30[i] = result.best_values["d"]
    mval30[i] = result.best_values["m"]

b = -26.5735278
c = -6.59234502
d = 4.03588632
m = 0.20863786
n = 31

bval31 = np.zeros(40)
cval31 = np.zeros(40)
dval31 = np.zeros(40)
mval31 = np.zeros(40)

for i in range(40):

    params = Parameters()

    if i < 9:
        params.add("b", value=b, min=-100, max=-10.0)
    else:
        params.add("b", value=b, min=-60, max=-10.0)
    params.add("c", value=c, min=-15.0, max=0)
    params.add("d", value=d)  # , min=1.14, max=1.16)
    if i < 10:
        params.add("m", value=m, min=0.0, max=0.4)
    else:
        params.add("m", value=m, min=0.0, max=0.5)

    logXi = xi31[:, i]
    logTe = np.log10(Te31[:, i])

    we = [1.0 for j in range(np.size(logTe))]

    for j in range(np.size(logTe)):

        if logXi[j] > -6.0 and i > 3:
            we[j] = 10.0
        if logXi[j] > -7.0 and i < 6:
            we[j] = 10.0
        if logXi[j] > -8.0 and logXi[j] < -7 and i < 1:
            we[j] = 10.0
        if logXi[j] > -5.5 and i < 1:
            we[j] = 100.0
        if logXi[j] > -4 and i > 6:
            we[j] = 30.0
        if logXi[j] > -3.5 and i < 6:
            we[j] = 10.0

    gmod = Model(interpPicog31)
    result = gmod.fit(logTe, x=logXi, params=params, weights=we)
    x = np.linspace(-8, -2, 100)
    y = interpPicog31(
        x,
        result.best_values["b"],
        result.best_values["c"],
        result.best_values["d"],
        result.best_values["m"],
    )

    bval31[i] = result.best_values["b"]
    cval31[i] = result.best_values["c"]
    dval31[i] = result.best_values["d"]
    mval31[i] = result.best_values["m"]

xi = np.linspace(-8, -2, 100)
ba = [-53.6260, -33.4958, -36.0245]
ca = [-5.8336, -5.4189, -5.2231]
da = [3.9421, 3.9088, 3.8720]
ma = [0.1220, 0.2039, 0.1527]
legend = [
    "5e20",
    "1e21",
    "1.5e21",
    "2e21",
    "2.5e21",
    "3e21",
    "3.5e21",
    "4.e21",
    "4.5e21",
    "5e21",
    "5.5e21",
    "6e21",
    "6.5e21",
    "7e21",
    "7.5e21",
    "8e21",
    "8.5e21",
    "9e21",
    "9.5e21",
    "1e22",
    "1.05e22",
    "1.1e22",
    "1.15e22",
    "1.2e22",
    "1.25e22",
    "1.3e22",
    "1.35e22",
    "1.4e22",
    "1.45e22",
    "1.5e22",
    "1.55e22",
    "1.6e22",
    "1.65e22",
    "1.7e22",
    "1.75e22",
    "1.8e22",
    "1.85e22",
    "1.9e22",
    "1.95e22",
    "2e22",
]
lines = ["-", "--", ":"]

fig = plt.figure()
k = 0
for i in [0, 1, 15]:
    plt.plot(
        xi,
        xiT(xi, ba[k], ca[k], da[k], ma[k]),
        c="0.5",
        linestyle=lines[k],
        label=legend[i],
    )
    plt.plot(
        xi,
        interpPicog29(xi, bval29[i], cval29[i], dval29[i], mval29[i]),
        "b",
        linestyle=lines[k],
    )
    plt.plot(
        xi,
        interpPicog30(xi, bval30[i], cval30[i], dval30[i], mval30[i]),
        "r",
        linestyle=lines[k],
    )
    plt.plot(
        xi,
        interpPicog31(xi, bval31[i], cval31[i], dval31[i], mval31[i]),
        "k",
        linestyle=lines[k],
    )
    k = k + 1
plt.ylim(1.3, 4.1)

# Add first legend:  only labeled data is included
leg1 = plt.legend(loc="upper left", prop={"size": 14})
# Add second legend for the maxes and mins.
# leg1 will be removed from figure
leg2 = plt.legend(
    ["Ercolano+2008", "\\texttt{Spec29}", "\\texttt{Spec30}", "\\texttt{Spec31}"],
    loc="lower right",
    prop={"size": 14},
)
plt.xlabel("$\\log_{10}(\\xi)$", size=14.0)
plt.ylabel("$\\log_{10}$(T [K])", size=14.0)
# Manually add the first legend back
fig.add_artist(leg1)
plt.savefig("Figure2.pdf", bbox_inches="tight", dpi=400)
