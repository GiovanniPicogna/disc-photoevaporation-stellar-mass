
Stellar parameters were obtained from the Siess et al. 2000 website with age = 1.00E+06 and metallicity : Z = 0.02 , no overshooting.

grid_001:
 ST       L         R        Teff     Mass  Disk_mass  Outer_rad 
 M6   8.559E-02  1.055E+00   2928.   0.100  2.67e-2    400 au
 
grid_002:
 ST       L         R        Teff     Mass  Disk_mass  Outer_rad
 M5   6.887E-01  2.310E+00   3360.   0.300  2.96e-2    400 au
 
grid_003:
 ST       L         R        Teff     Mass  Disk_mass  Outer_rad
 M1   9.288E-01  2.125E+00   3771.   0.500  3.69e-2    400 au
 
grid_004:
 ST       L         R        Teff     Mass  Disk_mass  Outer_rad
 K6   2.335E+00  2.615E+00   4278.   1.000  4.50e-2    400 au

Disk model parameters:
- Mass accretion rate: 1.0e-8
- alpha (Shakura & Sunyaev viscosity parameter) = 0.01 
- inclination = 60 degrees
- minimum grain size in disk atmosphere: 0.005 microns
- maximum grain size in disk atmosphere: 0.25 microns
- minimum grain size in disk midplane: 0.005 microns
- maximum grain size in disk midplane: 1 mm 
- epsilon (dust settling parameter) = 1 i.e., well-mixed dust

Here is a typical disk structure file name:
fort14.a0.01.irr.grid_001.e1.0.amax0p25_amax1mm.dat

The naming convention is as follows
fort14 - Fortran file name
a0.01 - alpha
irr - label we use to denoted this is the irradiated structure
grid_001 - label to specify set of model parameters
e1.0 - epsilon
amax0p25_amax1mm - maximum grain size in atmosphere and midplane, respectively

The structure is calculated for 79 points in radius and 480 points in height.  Within each file
there is a pattern that repeats 79 times for each of the 79 points in radius as follows:
You can ignore the first line.
In the second line, the first number (irad) tells you the number of the iteration (from 1-79).
In the third line, the first number tells you the radius in stellar radii the irad number corresponds to.
Then there are two more lines you can ignore.
Then there are 6 columns.  The first is z (height in the disk in stellar radii) and the third is temperature (in K).
Then it repeats.




