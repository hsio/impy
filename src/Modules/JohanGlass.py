# Write a time history of implosion areal density
# A. Zylstra 2013/10/25

from Implosion import *
from Resources.IO import *
from Resources.Constants import *
from numpy import arange
import numpy as np
import math
import csv
import os

# helper methods
def calcRhoR(impl, it, ir1, ir2):
    """Calculate the areal density between indices ir1 and ir2 at time index it in impl. Return in mg/cm2"""
    rR = 0
    for ir in range(ir1, ir2):
        rR += impl.rhoR(ir,it)
    return rR*1e3

def shellTempMethod(impl, it, method='simple'):
    """Calculate the unablated shell material based on Te less than 'average'.
    'simple' method calculates average T as a direct average of min/max temps
    'mass' uses a mass-averaged Te
    """
    i1 = impl.ir_fuel()+1
    i2 = impl.ir_max()-1

    rho = np.arange(i2-i1, dtype=np.float64)
    vol = np.arange(i2-i1, dtype=np.float64)
    Te = np.arange(i2-i1, dtype=np.float64)
    for ir in range(i1,i2):
        rho[ir-i1] = impl.rho(ir, it)
        vol[ir-i1] = impl.vol(ir, it)
        Te[ir-i1] = impl.Te(ir, it)
    #print(Te)

    # calculate total mass:
    mass = np.dot(vol, rho)
    # calculate max, min, temperature:
    MaxTe = np.max(Te)
    MinTe = np.min(Te)
    if method == 'simple':
        AvgTe = (MaxTe-MinTe)/2.
    elif method == 'massAvg':
        AvgTe = np.dot(Te, np.multiply(rho,vol)) / mass

    # calculate mass below the 'average' Te
    unablatedMass = 0.
    for i in range(len(rho)):
        if Te[i] < AvgTe:
            unablatedMass += rho[i]*vol[i]

    return unablatedMass / mass

def critDensityMethod(impl, it):
    """Calculate the unablated shell material based on ne>nc."""
    i1 = impl.ir_fuel()+1
    i2 = impl.ir_max()

    rho = np.arange(i2-i1, dtype=np.float64)
    vol = np.arange(i2-i1, dtype=np.float64)
    ne = np.arange(i2-i1, dtype=np.float64)
    for ir in range(i1,i2):
        rho[ir-i1] = impl.rho(ir, it)
        vol[ir-i1] = impl.vol(ir, it)
        ne[ir-i1] = impl.ne(ir, it)

    # calculate total mass:
    mass = np.dot(vol, rho)

    # critical density for 3w light:
    nc = 1.1e21 / 0.351**2.  # 1/cc

    # calculate mass with ne>nc
    unablatedMass = 0.
    for i in range(len(rho)):
        if ne[i] > nc:
            unablatedMass += rho[i]*vol[i]

    return unablatedMass / mass
    
def run(impl):
    # input sanity check:
    if not isinstance(impl,Implosion):
        print("WARNING: invalid input.")
        return
    
    File = csv.writer(open(os.path.join(OutputDir,'JohanGlass.csv'),'w'))
    File.writerow( ["t (s)", "Fuel rR (mg/cm2)", "Shell rR (mg/cm2)", "Total rR (mg/cm2)", "Shell Remaining (AvgTe method)", "Shell Remaining (MassAvgTe method)", "Shell Remaining (ne>nc method)"] )
    
    for it in range(impl.it_min(), impl.it_max()):
        ir1 = impl.ir_min()
        ir2 = impl.ir_fuel()
        ir3 = impl.ir_max()
        fuel = calcRhoR(impl, it, ir1, ir2)
        shell = calcRhoR(impl, it, ir2, ir3)

        # shell 'remaining' calculations
        AvgTe = shellTempMethod(impl, it)
        MassAvgTe = shellTempMethod(impl, it, method='massAvg')
        critDen = critDensityMethod(impl, it)

        # do output:
        t = impl.t(it)

        File.writerow( [ t , fuel , shell , (fuel+shell), AvgTe, MassAvgTe, critDen] )
    