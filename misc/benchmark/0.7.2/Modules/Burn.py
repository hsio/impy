# Calculate yields, Ti
# A. Zylstra 2012/08/17

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
def fT(impl, r, t):
    """Calculate 3Te fraction in implosion impl."""
    for i in range( len(impl.IonF(r,t)) ):
        if impl.IonA(r,t)[i] == 3 and impl.IonZ(r,t)[i] == 1:
            return impl.IonF(r,t)[i]
    return 0
    
# ------------------------------------
# Burn rate calculators
# ------------------------------------
def DDrate(impl, t):
    """Calculate the DD burn rate at time t (s)."""
    f1 = 0 # D fraction (atomic)
    ret = 0
    ret2 = 0 #Ti 'rate'
    
    for r in arange( impl.rmin(t) , impl.rfuel(t) , dr ):
        r1 = r + dr/2
        f1 = fD(impl,r1,t)
        temp = DD(impl.Ti(r1,t))*pow(impl.ni(r1,t),2)*(f1*f1/2)*4*math.pi*pow(r1,2)*dr
        ret += temp
        ret2 += temp*impl.Ti(r1,t)
    return [ret , ret2]
def D3Herate(impl, t):
    """Calculate the D3He burn rate at time t (s)."""
    f1 = 0 # D fraction (atomic)
    f2 = 0 # 3He fraction (atomic)
    ret = 0
    ret2 = 0 #Ti 'rate'

    for r in arange( impl.rmin(t) , impl.rfuel(t) , dr ):
        r1 = r + dr/2
        f1 = fD(impl,r1,t)
        f2 = f3He(impl,r1,t)
        temp = D3He(impl.Ti(r1,t))*pow(impl.ni(r1,t),2)*(f1*f2)*4*math.pi*pow(r1,2)*dr
        ret += temp
        ret2 += temp*impl.Ti(r1,t)
    return [ret , ret2]

def HeHerate(impl, t):
    """Calculate the 3He3He burn rate at time t (s)."""
    f1 = 0 # 3he fraction (atomic)
    ret = 0
    ret2 = 0 #Ti 'rate'

    for r in arange( impl.rmin(t) , impl.rfuel(t) , dr ):
        r1 = r + dr/2
        f1 = f3He(impl,r1,t)
        temp = HeHe(impl.Ti(r1,t))*pow(impl.ni(r1,t),2)*(f1*f1/2)*4*math.pi*pow(r1,2)*dr
        ret += temp
        ret2 += temp*impl.Ti(r1,t)
    return [ret , ret2]

def DTrate(impl, t):
    """Calculate the DT burn rate at time t (s)."""
    f1 = 0 # D fraction (atomic)
    f2 = 0 # T fraction (atomic)
    ret = 0
    ret2 = 0 #Ti 'rate'

    for r in arange( impl.rmin(t) , impl.rfuel(t) , dr ):
        r1 = r + dr/2
        f1 = fD(impl,r1,t)
        f2 = fT(impl,r1,t)
        temp = DT(impl.Ti(r1,t))*pow(impl.ni(r1,t),2)*(f1*f2)*4*math.pi*pow(r1,2)*dr
        ret += temp
        ret2 += temp*impl.Ti(r1,t)
    return [ret , ret2]

def TTrate(impl, t):
    """Calculate the TT burn rate at time t (s)."""
    f1 = 0 # T fraction (atomic)
    ret = 0
    ret2 = 0 #Ti 'rate'

    for r in arange( impl.rmin(t) , impl.rfuel(t) , dr ):
        r1 = r + dr/2
        f1 = fT(impl,r1,t)
        temp = TT(impl.Ti(r1,t))*pow(impl.ni(r1,t),2)*(f1*f1/2)*4*math.pi*pow(r1,2)*dr
        ret += temp
        ret2 += temp*impl.Ti(r1,t)
    return [ret , ret2]
 
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
    YDT = 0
    YTT = 0
    # ion temps (burn-averaged)
    TiDD = 0
    TiD3He = 0
    TiHeHe = 0
    TiDT = 0
    TiTT = 0
    # Bang (peak emission) times
    BTDD = 0
    BTD3He = 0
    BTHeHe = 0
    BTDT = 0
    BTTT = 0
    PeakRateDD = 0
    PeakRateD3He = 0
    PeakRateHeHe = 0
    PeakRateDT = 0
    PeakRateTT = 0
    
    # output files
    rateFile = csv.writer(open(os.path.join(OutputDir,'BurnRate.csv'),'w'))
    rateFile.writerow( ["t (s)", "DD (1/s)", "D3He (1/s)", "3He3He (1/s)", "DT (1/s)", "TT (1/s)"] )
    yieldFile = csv.writer(open(os.path.join(OutputDir,'Yield.csv'),'w'))
    TiFile = csv.writer(open(os.path.join(OutputDir,'Ti.csv'),'w'))
    BTFile = csv.writer(open(os.path.join(OutputDir,'BangTime.csv'),'w'))
    
    #iterate over all time:
    for t in list(arange(impl.tmin(), impl.tmax(), dt)):
        #DD
        [dDD, dTiDD] = DDrate(impl,t)
        YDD += dDD*dt
        TiDD += dTiDD*dt
        if dDD > PeakRateDD:
            BTDD = t
            PeakRateDD = dDD
        #D3He
        [dD3He, dTiD3He] = D3Herate(impl,t)
        YD3He += dD3He*dt
        TiD3He += dTiD3He*dt
        if dD3He > PeakRateD3He:
            BTD3He = t
            PeakRateD3He = dD3He
        #3He3He
        [d3He3He, dTi3He3He] = HeHerate(impl,t)
        YHeHe += d3He3He*dt
        TiHeHe += dTi3He3He*dt
        if d3He3He > PeakRateHeHe:
            BTHeHe = t
            PeakRateHeHe = d3He3He
        #DD
        [dDT, dTiDT] = DTrate(impl,t)
        YDT += dDT*dt
        TiDT += dTiDT*dt
        if dDT > PeakRateDT:
            BTDT = t
            PeakRateDT = dDT
        #DD
        [dTT, dTiTT] = TTrate(impl,t)
        YTT += dTT*dt
        TiTT += dTiTT*dt
        if dTT > PeakRateTT:
            BTTT = t
            PeakRateTT = dTT
        #output
        rateFile.writerow( [t, dDD, dD3He, d3He3He, dDT, dTT] )
    
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
    if YDT > 0:
        TiDT = TiDT / YDT
        print("DT yield = " + '{:.2e}'.format(YDT))
        print("DT Ti = " + '{:.2f}'.format(TiDT))
        print("DT BT = " + '{:.2e}'.format(BTDT))
        yieldFile.writerow( ["DT",YDT] )
        TiFile.writerow( ["DT",TiDT] )
        BTFile.writerow( ["DT",BTDT] )
    if YTT > 0:
        TiTT = TiTT / YTT
        print("TT yield = " + '{:.2e}'.format(YTT))
        print("TT Ti = " + '{:.2f}'.format(TiTT))
        print("TT BT = " + '{:.2e}'.format(BTTT))
        yieldFile.writerow( ["TT",YTT] )
        TiFile.writerow( ["TT",TiTT] )
        BTFile.writerow( ["TT",BTTT] )
        