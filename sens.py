
#!/usr/local/bin/python3

#
# Author: Evan Keane
# Date: 11/05/2020
# Description: Determines MeerKAT sensitivity curves etc. for various 
#              user-defined input configurations
#
# TODO
# 1. tidy up, add opacity, zenith angle and temp transforms etc.
# 2. read in array configuration CSV files
# 3. read in Haslam map for specific and average gl,gb values and ranges
# 4. add a check that the zenith angle is possible for the given sky coords!
#

# Load some useful packages
import argparse
import sys
import math as m
import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
from functions import *

h_over_k = 4.8e-11 # Planck's constant divided by Boltzmann's constant in s*K
kB = 1380.0        # Boltzmann's constant in Jy*m^2*K^-1
nfreqs = 2000      # number of points at which to sample frequency for output plots etc.

# Parse command line arguments
parser = argparse.ArgumentParser()
parser.add_argument('-radius', type=float, dest='radius', help='choose distance from the array centre, in km, for chosen sub-array (default: entire array, 4.2km)', default=4.2)
parser.add_argument('-glgb', nargs=2, type=float, dest='coord', help='enter specific Galactic coordinates to use (default: gl=180.0, gb=-90.0) NOT DONE YET', default=[180.0,-90.0])
parser.add_argument('-gallos', dest='gal', help='choose either 10th, 50th or 90th percentile value for the galaxy contribution to the sky temperature (low/medium/high, default: low)', default='low')
parser.add_argument('-pwv', dest='pwv', help='choose either 5mm, 10mm or 20mm for the PWV value for choosing (a) the zenith opacity, and (b) the atmospheric temperature contribution to the sky temperature (low/medium/high, default: low)', default="low")
#parser.add_argument('-tel', dest='tel', help='Which telescopes to include - options are all, low, mid, ska (meaning SKA1-Mid dishes only) or mk (default: all)', default="all")
parser.add_argument('-nelements', type=int, dest='nelements', help='choose the inner nelements elements (default: entire array)', default=64)
parser.add_argument('-o', dest='output', help='choose the type of output - plot, file or both (default: plot)', default="plot")
parser.add_argument('-zenith', type=float, dest='zenith', help='choose a zenith angle in degrees (default: 0.0)', default=0.0)
parser.add_argument('--version', action='version', version='%(prog)s 0.0.1')
args = parser.parse_args()

# Set values from command line inputs
#tel = args.tel
gl = args.coord[0]
gb = args.coord[1]
gal = args.gal
output = args.output
plot = True
if output != "plot":
    plot = False
zenith = args.zenith
pwv = args.pwv
radius = args.radius
nelements = args.nelements

# Get the effective collecting area
#Aeff_SKA = get_aeff("SKA",plot)
Aeff_MK  = get_aeff("MeerKAT",plot)
#Aeff_Eff  = get_aeff("Effelsberg",plot)

# System Temperature
#get_tsys("SKA")
#get_tsys("MeerKAT")
#Tsys_SKA, f = get_tsys("SKA",gal,pwv,zenith,plot)
Tsys_MK, f  = get_tsys("MeerKAT",gal,pwv,zenith,plot)
#Tsys_Eff, f = get_tsys("Effelsberg",gal,pwv,zenith,plot)

# Gain - single dish
# in m^2/K (Tsys, i.e. LoS dependent)
f = np.logspace(np.log10(0.35),np.log10(50),200)
plt.grid(True)
#plt.semilogx(f,Aeff_SKA(f)/Tsys_SKA(f)) #(Trcv(f)+Tspill(f)+Tsky(f)))
plt.semilogx(f,Aeff_MK(f)/Tsys_MK(f)) #(Trcv(f)+Tspill(f)+Tsky(f)))
plt.title("Gain - single Dish")
plt.ylabel("Aeff/Tsys (m^2/K)")
plt.xlabel("Frequency (GHz)")
plt.show()

# in K/Jy (LoS independent but no indication of Tsys impact on performance)
f = np.logspace(np.log10(0.35),np.log10(50),200)
plt.grid(True)
#plt.semilogx(f,Aeff_SKA(f)/(2*kB))
plt.semilogx(f,Aeff_MK(f)/(2*kB))
plt.title("Gain - single Dish")
plt.ylabel("Aeff/(2*kB) (K/Jy)")
plt.xlabel("Frequency (GHz)")
plt.show()

# Gain - user-requested sub-array
# need to read in the configs, for now just do whole array
# 
# Also compare some actually relevant telescopes, maybe a different flag for imaging- and NIP-relevant ones to show
f = np.logspace(np.log10(0.35),np.log10(50),200)
plt.grid(True)
plt.loglog(f,64.0*(Aeff_MK(f)/Tsys_MK(f)), label='my numbers')
plt.title("Gain - entire MeerKAT array")
plt.ylabel("Aeff/Tsys (m^2/K)")
plt.xlabel("Frequency (GHz)")

plt.show()

Nmk  = 64
if ((radius < 4.2) or (nelements < 64)):

    # Read in the full 64-dish array configuration
    array = np.genfromtxt("./Configuration/MK_dist_metres.txt",dtype=[('name','S6'),('xm','f8'),('ym','f8')],skip_header=0,usecols=(1,2,3)) # NB this is an array of tuples, not a 2-D array, due to carrying the name label for each dish
    dist = np.zeros(np.size(array))    # work out the distance from centre in km
    for i in range(0,np.size(array)):
        dist[i] = 0.001*np.sqrt(array[i][1]*array[i][1]+array[i][2]*array[i][2])
    # Sort by distance from array centre
    array = array[np.argsort(dist)]
    dist = dist[np.argsort(dist)]
    print(dist)
    # Choose the sub-array of interest, first using the nelements flag
    # Implicitly assumed that the innermost nelements are wanted
    subarray = array[0:nelements]
    subdist  = dist[0:nelements]
    if (radius < 4.2):
        # if radius is also specified re-work out nelements, 
        # i.e. if both set the smaller sub-array is taken
        # if radius > subdist[nelements-1]    --> DO NOTHING
        if (radius < subdist[nelements-1]): # --> Determine a smaller nelements
            subarray = array[np.where( dist < radius)]
            subdist = dist[np.where( dist < radius)]
            nelements = subarray.shape[0]
#    Nmk = 0
#    for i in range(0,nelements):
#        if subarray[i][0][0] == 'M':
#            Nmk +=1
    print ("Considering a radius of %.1f km"%(subdist[nelements-1]))
    if nelements % 4 != 0: # CBF constraint, can add 64, 60, 56, ... dishes
        print ("%d %d %d"%(nelements, (nelements/4), (nelements%4)))
        print ("There are %d dishes within that radius."%(nelements)) 
        nelements = (nelements//4)*4
        print ("But Beamformer adds in fours so will only use %d dishes."%(nelements))
    print ("Considering %d MeerKAT dishes"%(nelements))
    plt.grid(True)
    plt.loglog(f,nelements*(Aeff_MK(f)/Tsys_MK(f)))
    plt.title("Gain - subarray radius %.1f km - %d MeerKAT dishes"%(radius,nelements))
    plt.ylabel("Aeff/Tsys (m^2/K)")
    plt.xlabel("Frequency (GHz)")
    plt.show()
    plt.grid(True)
    plt.loglog(f,nelements*(Aeff_MK(f)/(2*kB)))
    plt.title("Gain - subarray radius %.1f km - %d MeerKAT"%(radius,nelements))
    plt.ylabel("Aeff/(2*kB) (K/Jy)")
    plt.xlabel("Frequency (GHz)")
    plt.show()

# Prints out sensitivity as a function of frequency
for i in range (0,f.size):
    print (f[i], Nmk*(Aeff_MK(f[i])/Tsys_MK(f[i])))

sys.exit()
