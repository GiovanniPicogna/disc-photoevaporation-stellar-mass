import os
import h5py as h 
import numpy as np

#############################################
# Functions to calculate the Temperature from the ionization parameter
#############################################

def xiT(x, b, c, d, m):
    T = d + (1.5 - d) / (1 + (x / c) ** b) ** m
    return T


def interpPicog29(x, b, c, d, m):
    return d + (1.5 - d) / (1.0 + (x / c) ** b) ** m


def interpPicog30(x, b, c, d, m):
    return d + (1.6 - d) / (1.0 + (x / c) ** b) ** m


def interpPicog31(x, b, c, d, m):
    return d + (1.7 - d) / (1.0 + (x / c) ** b) ** m

#############################################
# Functions to read PLUTO outputs
#############################################

def getFilenames():
    filenames = [x for x in os.listdir("./") if ".h5" in x]
    filenames.sort()
    return filenames


def getVar(filename, step, variable):
    h5 = h.File(filename, "r")
    returnData = h5["Timestep_" + str(step) + "/vars"][variable][:]
    h5.close()
    return returnData


def getGridCell(filename=None, all=1):
    if not (filename):
        filename = getFilenames()[0]
    h5 = h.File(filename, "r")
    if all:
        x = h5["cell_coords"]["X"][:]
        y = h5["cell_coords"]["Y"][:]
        z = h5["cell_coords"]["Z"][:]
    else:
        x = h5["cell_coords"]["X"]
        y = h5["cell_coords"]["Y"]
        z = h5["cell_coords"]["Z"]
    x = x.astype("float64")
    y = y.astype("float64")
    z = z.astype("float64")
    return x, y, z

#############################################
# Functions defining the mass-loss rates
#############################################

def MdotLxPicogna(Lx1, Lx2):
    AL = -2.7326
    BL = 3.3307
    CL = -2.9868e-3
    DL = -7.2580
    mdot1 = 10 ** (AL * np.exp(((np.log(np.log10(Lx1)) - BL) ** 2) / CL) + DL)
    mdot2 = 10 ** (AL * np.exp(((np.log(np.log10(Lx2)) - BL) ** 2) / CL) + DL)
    return mdot1 / mdot2

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

def Lx(*Mstar):
    for x in Mstar:
        maxLx = 10 ** (1.42 * np.log10(x) + 30.37)
        minLx = 10 ** (1.66 * np.log10(x) + 30.25)
        aveLx = 10 ** (1.54 * np.log10(x) + 30.31)
        return minLx, aveLx, maxLx


#############################################
# Miscellaneous
#############################################

def sci_notation(num, decimal_digits=1, precision=None, exponent=None):
    """
    Returns a string representation of the scientific
    notation of the given number formatted for use with
    LaTeX or Mathtext, with specified number of significant
    decimal digits and precision (number of decimal digits
    to show). The exponent to be used can also be specified
    explicitly.
    """
    if exponent is None:
        exponent = int(np.floor(np.log10(abs(num))))
    coeff = round(num / float(10 ** exponent), decimal_digits)
    if precision is None:
        precision = decimal_digits

    return r"${0:.{2}f}\cdot10^{{{1:d}}}$".format(coeff, exponent, precision)\

def moving_average(x, w):
    return np.convolve(x, np.ones(w), "valid") / w

