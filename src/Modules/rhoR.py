# Write a time history of implosion areal density
# A. Zylstra 2013/04/13

from Implosion import *
from Resources.IO import *
from Resources.Constants import *
from numpy import arange
import math
import csv
import os

# helper method
def calcRhoR(impl, it, ir1, ir2):
    """Calculate the areal density between indices ir1 and ir2 at time index it in impl. Return in mg/cm2"""
    rR = 0
    for ir in range(ir1, ir2):
        rR += impl.rhoR(ir,it)
    return rR*1e3
    
def run(impl):
    # input sanity check:
    if not isinstance(impl,Implosion):
        print("WARNING: invalid input.")
        return
    
    File = csv.writer(open(os.path.join(OutputDir,'rhoR.csv'),'w'))
    File.writerow( ["t (s)", "Fuel rR (mg/cm2)", "Shell rR (mg/cm2)", "Total rR (mg/cm2)"] )
    
    for it in range(impl.it_min(), impl.it_max()):
        ir1 = impl.ir_min()
        ir2 = impl.ir_fuel()
        ir3 = impl.ir_max()
        fuel = calcRhoR(impl, it, ir1, ir2)
        shell = calcRhoR(impl, it, ir2, ir3)
        t = impl.t(it)
        File.writerow( [ t , fuel , shell , (fuel+shell)] )
    