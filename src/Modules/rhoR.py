# Write a time history of implosion areal density
# A. Zylstra 2012/08/16

from Implosion import *
from Resources.IO import *
from Resources.Constants import *
from numpy import arange
import math
import csv
import os

dt = 50e-12
dr = 1e-4 #1um

# helper method
def calcRhoR(impl, t, r1, r2):
    """Calculate the areal density between r1 and r2 at time t in impl. Return in mg/cm2"""
    rR = 0
    for r in list(arange(r1, r2, dr)):
        rR += impl.rho(r+dr/2,t)*dr
    return rR*1e3
    
def run(impl):
    # input sanity check:
    if not isinstance(impl,Implosion):
        print("WARNING: invalid input.")
        return
    
    File = csv.writer(open(os.path.join(OutputDir,'rhoR.csv'),'w'))
    File.writerow( ["t (ns)", "Fuel rR (mg/cm2)", "Shell rR (mg/cm2)", "Total rR (mg/cm2)"] )
    
    for t in list(arange(impl.tmin(), impl.tmax(), dt)):
        r1 = impl.rmin(t)
        r2 = impl.rfuel(t)
        r3 = impl.rmax(t)
        fuel = calcRhoR(impl, t, r1, r2)
        shell = calcRhoR(impl, t, r2, r3)
        File.writerow( [ t , fuel , shell , (fuel+shell)] )
    