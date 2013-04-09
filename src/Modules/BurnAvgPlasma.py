# Calculate burn-averaged plasma parameters
# A. Zylstra 2012/09/19

from Implosion import *
from Resources.IO import *
from Resources.Plasma import *
from Resources import Fusion
from numpy import arange
from numpy import zeros
import math
import csv
import os

# integration step sizes
dt = 10e-12 #10ps
dr = 5e-4 #5um

# global implosion
impl = 0

# weighting function
WeightFunction = Fusion.D3He
A1 = 2
A2 = 3
Z1 = 1
Z2 = 2

# ------------------------------------
# Fuel fraction helpers
# ------------------------------------
def f(r, t, A, Z):
    """Find fraction for fuel ion with A and Z."""
    if Z == 1:
        if A == 2:
            return fD(r,t)
        if A == 3:
            return fT(r,t)
    if Z == 2:
        if A == 3:
            return f3He(r,t)
    return 0
def fD(r, t):
    """Calculate D fraction in implosion impl."""
    for i in range( len(impl.IonF(r,t)) ):
        if impl.IonA(r,t)[i] == 2 and impl.IonZ(r,t)[i] == 1:
            return impl.IonF(r,t)[i]
    return 0
def f3He(r, t):
    """Calculate 3He fraction in implosion impl."""
    for i in range( len(impl.IonF(r,t)) ):
        if impl.IonA(r,t)[i] == 3 and impl.IonZ(r,t)[i] == 2:
            return impl.IonF(r,t)[i]
    return 0
def fT(r, t):
    """Calculate 3Te fraction in implosion impl."""
    for i in range( len(impl.IonF(r,t)) ):
        if impl.IonA(r,t)[i] == 3 and impl.IonZ(r,t)[i] == 1:
            return impl.IonF(r,t)[i]
    return 0
    
# ------------------------------------
# Burn rate calculators
# ------------------------------------
def rate(t):
    """Calculate the 'rate' at time t (s). Returns [lD , ND, PlasmaF , MFP , tau_i , R_cs , R_vt ]"""
    # return values
    Yield = lD = ND = PlasmaF = MFP = tau_i = R_cs = R_vt = 0
    rfuel = impl.rfuel(t)
    
    # iterate over radius
    for r in arange( impl.rmin(t) , impl.rmax(t) , dr):
        r1 = r + dr/2
        f1 = f(r1,t,A1,Z1)
        f2 = f(r1,t,A2,Z2)
        dblcount = 1
        if (A1 == A2) and (Z1 == Z2):
            dblcount = 2 #account for factor of 2 if reactants are identical

        # plasma conditions:
        ni = impl.ni(r1,t)
        Ti = impl.Ti(r1,t)
        ne = impl.ne(r1,t)
        Te = impl.Te(r1,t)
        A = impl.Abar(r1,t)
        Z = impl.Zbar(r1,t)
        
        # weight function
        w = WeightFunction(Ti)*pow(ni,2)*(f1*f2/dblcount)*4*math.pi*pow(r1,2)*dr
        
        # add weighted incremental values
        Yield += w
        lD += w*LambdaD(ne,Te)
        ND += w*PlasmaParameter(ne,Te)
        PlasmaF += w*fpe(ne)
        MFP += w*IonMFP(ni,Ti,Z,A)
        tau_i += w*Taui(ni,Ti,Z,A)
        R_cs += w * rfuel / impl.c(r,t)
        R_vt += w * rfuel / uTherm(Ti,A)
    
    # return array
    return [Yield, lD , ND, PlasmaF , MFP , tau_i , R_cs , R_vt ]
        
 
# ------------------------------------
# Main method
# ------------------------------------
def run(i):
    """Calculate total yield."""
    # input sanity check:
    if not isinstance(i,Implosion):
        print("WARNING: invalid input.")
        return
    global impl
    impl = i # global implosion variable for this module
        
    # things to calculate
    Yield = 0
    lD = 0 # debye length in cm
    ND = 0 # plasma parameter
    PlasmaF = 0 # plasma frequency
    MFP = 0 # ion MFP
    tau_i = 0 # ion-ion collision time
    R_vt = 0 # ion transit time R/vt
    R_cs = 0 # radius / sound speed
    
    # output files
    file = csv.writer(open(os.path.join(OutputDir,'BurnAvgPlasma.csv'),'w'))
    
    #iterate over all time:
    for t in list(arange(impl.tmin(), impl.tmax(), dt)):
        dRate = rate(t)
        Yield += dRate[0]
        lD += dRate[1]
        ND += dRate[2]
        PlasmaF += dRate[3]
        MFP += dRate[4]
        tau_i += dRate[5]
        R_vt += dRate[6]
        R_cs += dRate[7]
        
    #normalization
    lD = lD/Yield
    ND = ND/Yield
    PlasmaF = PlasmaF / Yield
    MFP = MFP / Yield
    tau_i = tau_i / Yield
    R_vt = R_vt / Yield
    R_cs = R_cs / Yield
    
    file.writerow([ " Debye length",lD,"cm"])
    file.writerow(["Plasma parameter",ND])
    file.writerow([ "Plasma freq." , PlasmaF , "Hz" ])
    file.writerow([ "Mean Free Path" , MFP , "cm" ])
    file.writerow([ "Ion collision time" , tau_i , "s" ])
    file.writerow([ "R/vt" , R_vt , "s" ])
    file.writerow([ "R/cs" , R_cs , "s" ])
