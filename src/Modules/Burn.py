# Calculate yields, Ti
# A. Zylstra 2012/08/14

from Implosion import *
from Resources.IO import *
from Resources.Fusion import *
from numpy import arange
import math
import csv
import os

# integration step sizes
dt = 10e-12 #10ps
dr = 5e-4 #5um

# ------------------------------------
# Fuel fraction helpers
# ------------------------------------
def fD(impl, r, t):
    """Calculate D fraction in implosion impl."""
    for i in range( len(impl.IonF(r,t)) ):
        if impl.IonA(r,t)[i] == 2 and impl.IonZ(r,t)[i] == 1:
            return impl.IonF(r,t)[i]
    return 0
def f3He(impl, r, t):
    """Calculate 3He fraction in implosion impl."""
    for i in range( len(impl.IonF(r,t)) ):
        if impl.IonA(r,t)[i] == 3 and impl.IonZ(r,t)[i] == 2:
            return impl.IonF(r,t)[i]
    return 0
    
# ------------------------------------
# Burn rate calculators
# ------------------------------------
def DDrate(impl, t):
    """Calculate the DD burn rate at time t (s)."""
    f1 = 0 # D fraction (atomic)
    ret = 0
    
    for r in arange( impl.rmin(t) , impl.rfuel(t) , dr ):
        r1 = r + dr/2
        f1 = fD(impl,r1,t)
        ret += DD(impl.Ti(r1,t))*pow(impl.ni(r1,t),2)*(f1*f1/2)*4*math.pi*pow(r1,2)*dr
    return ret
def D3Herate(impl, t):
    """Calculate the D3He burn rate at time t (s)."""
    f1 = 0 # D fraction (atomic)
    f2 = 0 # 3He fraction (atomic)
    ret = 0

    for r in arange( impl.rmin(t) , impl.rfuel(t) , dr ):
        r1 = r + dr/2
        f1 = fD(impl,r1,t)
        f2 = f3He(impl,r1,t)
        ret += D3He(impl.Ti(r1,t))*pow(impl.ni(r1,t),2)*(f1*f2)*4*math.pi*pow(r1,2)*dr
    return ret

def HeHerate(impl, t):
    """Calculate the 3He3He burn rate at time t (s)."""
    f1 = 0 # 3he fraction (atomic)
    ret = 0

    for r in arange( impl.rmin(t) , impl.rfuel(t) , dr ):
        r1 = r + dr/2
        f1 = f3He(impl,r1,t)
        ret += HeHe(impl.Ti(r1,t))*pow(impl.ni(r1,t),2)*(f1*f1/2)*4*math.pi*pow(r1,2)*dr
    return ret


# ------------------------------------
# Burn-averaged Ti "rate" calculators
# ------------------------------------
def DDTirate(impl, t):
    """Calculate the DD burn rate at time t (s)."""
    f1 = 0 # D fraction (atomic)
    ret = 0
    
    for r in arange( impl.rmin(t) , impl.rfuel(t) , dr ):
        r1 = r + dr/2
        f1 = fD(impl,r1,t)
        ret += impl.Ti(r1,t)*DD(impl.Ti(r1,t))*pow(impl.ni(r1,t),2)*(f1*f1/2)*4*math.pi*pow(r1,2)*dr
    return ret
def D3HeTirate(impl, t):
    """Calculate the D3He burn rate at time t (s)."""
    f1 = 0 # D fraction (atomic)
    f2 = 0 # 3He fraction (atomic)
    ret = 0
    
    for r in arange( impl.rmin(t) , impl.rfuel(t) , dr ):
        r1 = r + dr/2
        f1 = fD(impl,r1,t)
        f2 = f3He(impl,r1,t)
        ret += impl.Ti(r1,t)*D3He(impl.Ti(r1,t))*pow(impl.ni(r1,t),2)*(f1*f2)*4*math.pi*pow(r1,2)*dr
    return ret

def HeHeTirate(impl, t):
    """Calculate the 3He3He burn rate at time t (s)."""
    f1 = 0 # 3He fraction (atomic)
    ret = 0

    for r in arange( impl.rmin(t) , impl.rfuel(t) , dr ):
        r1 = r + dr/2
        f1 = f3He(impl,r1,t)
        ret += impl.Ti(r1,t)*HeHe(impl.Ti(r1,t))*pow(impl.ni(r1,t),2)*(f1*f1/2)*4*math.pi*pow(r1,2)*dr
    return ret
  
# ------------------------------------
# Main method
# ------------------------------------
def run(impl):
    """Calculate total yield."""
    # input sanity check:
    if not isinstance(impl,Implosion):
        print("WARNING: invalid input.")
        return
        
    
    # Yields
    YDD = 0
    YD3He = 0
    YHeHe = 0
    # ion temps (burn-averaged)
    TiDD = 0
    TiD3He = 0
    TiHeHe = 0
    # Bang (peak emission) times
    BTDD = 0
    BTD3He = 0
    BTHeHe = 0
    PeakRateDD = 0
    PeakRateD3He = 0
    PeakRateHeHe = 0
    
    # output files
    rateFile = csv.writer(open(os.path.join(OutputDir,'BurnRate.csv'),'w'))
    rateFile.writerow( ["t (s)", "DD (1/s)", "D3He (1/s)", "3He3He (1/s)"] )
    yieldFile = csv.writer(open(os.path.join(OutputDir,'Yield.csv'),'w'))
    TiFile = csv.writer(open(os.path.join(OutputDir,'Ti.csv'),'w'))
    BTFile = csv.writer(open(os.path.join(OutputDir,'BangTime.csv'),'w'))
    
    #iterate over all time:
    for t in list(arange(impl.tmin(), impl.tmax(), dt)):
        #DD
        dDD = DDrate(impl,t)
        YDD += dDD*dt
        TiDD += DDTirate(impl,t)*dt
        if dDD > PeakRateDD:
            BTDD = t
            PeakRateDD = dDD
        #D3He
        dD3He = D3Herate(impl,t)
        YD3He += dD3He*dt
        TiD3He += D3HeTirate(impl,t)*dt
        if dD3He > PeakRateD3He:
            BTD3He = t
            PeakRateD3He = dD3He
        #3He3He
        d3He3He = HeHerate(impl,t)
        YHeHe += d3He3He*dt
        TiHeHe += HeHeTirate(impl,t)*dt
        if d3He3He > PeakRateHeHe:
            BTHeHe = t
            PeakRateHeHe = d3He3He
        #output
        rateFile.writerow( [t, dDD, dD3He, d3He3He] )
        
    #If there is yield for a species, do output:
    if YDD > 0:
        TiDD = TiDD / YDD
        print("DD yield = " + '{:.2e}'.format(YDD))
        print("DD Ti = " + '{:.2f}'.format(TiDD))
        print("DD BT = " + '{:.2e}'.format(BTDD))
        yieldFile.writerow( ["DD",YDD] )
        TiFile.writerow( ["DD",TiDD] )
        BTFile.writerow( ["DD",BTDD] )
    if YD3He > 0:
        TiD3He = TiD3He / YD3He
        print("D3He yield = " + '{:.2e}'.format(YD3He))
        print("D3He Ti = " + '{:.2f}'.format(TiD3He))
        print("D3He BT = " + '{:.2e}'.format(BTD3He))
        yieldFile.writerow( ["D3He",YD3He] )
        TiFile.writerow( ["D3He",TiD3He] )
        BTFile.writerow( ["D3He",BTD3He] )
    if YHeHe > 0:
        TiHeHe = TiHeHe / YHeHe
        print("3He3He yield = " + '{:.2e}'.format(YHeHe))
        print("3He3He Ti = " + '{:.2f}'.format(TiHeHe))
        print("3He3He BT = " + '{:.2e}'.format(BTHeHe))
        yieldFile.writerow( ["3He3He",YHeHe] )
        TiFile.writerow( ["3He3He",TiHeHe] )
        BTFile.writerow( ["3He3He",BTHeHe] )
        