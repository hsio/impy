# Calculate yields, Ti, BT using Molvig reduced reactivit
# A. Zylstra 2012/08/17

from Implosion import *
from Resources.IO import *
from Resources.Fusion import *
from Resources.Constants import *
import numpy
from scipy.integrate import quad
import scipy.interpolate
import math
import csv
import os

# integration step sizes
dt = 10e-12 #10ps
dr = 5e-4 #5um

#some interpolators
LeffInt = 0
NkArr = numpy.array( [] , float)

# -----------------------------------------------
#      Molvig-Knudsen reactivities
# -----------------------------------------------
#Helper functions:
def phi(theta, r, R):
    """Helper function for Molvig L calculations."""
    return ( math.pi/2.0 - theta - math.asin(r*math.sin(theta)/R) )
def L(theta, r, R):
    """spherical geometry L helper function."""
    return math.sqrt( math.pow(r-R*math.cos(phi(theta,r,R)),2) + math.pow(R*math.sin(phi(theta,r,R)),2) )
def integrand1(theta,r,R):
    """integrand for Leff calculation"""
    return ( math.sin(theta) / math.pow(L(theta,r,R),2) )
def LeffInit(impl):
    """Precompute Leff vs R,r and set up interpolation."""
    rR = []
    LR = []
    for i in list(numpy.arange(0, 1, 1e-3)):
        rR.append( i )
        LR.append( math.sqrt( 2. / quad(integrand1, 0, math.pi, args=(i,1))[0] ) )
    global LeffInt
    LeffInt = scipy.interpolate.InterpolatedUnivariateSpline(rR, LR)
def Leff(R, r, t):
    """Effective scale length in spherical geometry. R is shell radius."""
    if r >= R:
        return 0
    return R*LeffInt(r/R)
    
# Knudsen number:
def NkInit(impl):
    """Precompute Nk vs r,t."""
    global NkArr
    ir = math.ceil(max( impl.rfuel(impl.tmin()) , impl.rfuel(impl.tmax()) ) / dr)
    it = math.ceil(( impl.tmax() - impl.tmin() ) / dt)
    NkArr = numpy.ndarray( shape=(ir, it) , dtype=float )
    for i in range(int(ir)):
        for j in range(int(it)):
            r = i*dr
            t = impl.tmin() + j*dt
            lD = 743*math.sqrt(1000*impl.Te(r,t) / impl.ne(r,t))
            b90 = 2*math.pow(impl.Zbar(r,t)*e,2)/(3*kB*1000*11600*impl.Ti(r,t))
            LL = max( math.log(lD/b90) , 1 )
            #LL = max( LogL(impl.ni(r,t) , impl.Ti(r,t) , impl.Zbar(r,t) , impl.Abar(r,t) ) , 1)
            L = Leff(impl.rfuel(t),r,t)
            if L > 0:
                temp = (math.sqrt(0.26)/math.pi)*pow(1000*11600*kB*impl.Ti(r,t),2) / (impl.ni(r,t)*LL*L*math.pow(e,4))
            else:
                temp = 1e9
            if temp < 0:
                print("UH OH negative Nk!!!")
                print(LL)
                print(L)
                print(impl.Ti(r,t))
                print(impl.ni(r,t))
                #input("foo")
            NkArr[i,j] = temp 
def Nk(impl, r, t, Z1, Z2):
    """Knudsen number for fusion with ion charges Z1,Z2."""
    ir = int( r / dr )
    it = int( ( t - impl.tmin() ) / dt )
    return  ( NkArr[ir,it] / (pow(Z1*Z2,2)) )
    
# Knudsen distribution function:
def fK(impl, En,Nk,r,t):
    """Knudsen distribution function."""
    #normalized energy
    eps = En / impl.Ti(r,t)
    #evaluated as prefactor and exponential
    p1 = 2 / math.sqrt( math.pi + Nk*math.pow(eps,1.5) )
    p2 = math.exp( -1.0*(eps+0.8*Nk*math.pow(eps,2.5)+0.32*math.pow(Nk*eps*eps,2))/(1 + 0.8*Nk*math.pow(eps,1.5)) )
    return p1*p2
    
# DD reactivity:
def integrand2(En, r, t, impl):
    """Integrand for Molvig-Knudsen DD reactivity."""
    NkDD = Nk(impl, r, t, 1, 1)
    return sigmaDDn(En)*c*math.sqrt(En/(1000*938))*fK(impl,En,NkDD,r,t)*math.sqrt(En)
def svDDMolvig(impl, r,t):
    """Calculate DD reactivity (Molvig-Knudsen)."""
    R = impl.rfuel(t)
    if r > R or impl.Ti(r,t) < 0.5:
        return 0
    return quad(integrand2, 0, 1000, args=(r,t,impl))[0]
# D3He reactivity:
def integrand3(En, r, t, impl):
    """Integrand for Molvig-Knudsen D3He reactivity."""
    NkDHe = Nk(impl, r, t, 1, 2)
    return sigmaD3He(En)*c*math.sqrt(En/(1250*938))*fK(impl,En,NkDHe,r,t)*math.sqrt(En)
def svD3HeMolvig(impl, r,t):
    """Calculate D3He reactivity (Molvig-Knudsen)."""
    R = impl.rfuel(t)
    if r > R or impl.Ti(r,t) < 0.5:
        return 0
    return quad(integrand3, 0, 1000, args=(r,t,impl))[0]
# 3He3He reactivity:
def integrand4(En, r, t, impl):
    """Integrand for Molvig-Knudsen 3He3He reactivity."""
    NkHeHe = Nk(impl, r, t, 2, 2)
    return sigmaHeHe(En)*c*math.sqrt(En/(1500*938))*fK(impl,En,NkHeHe,r,t)*math.sqrt(En)
def svHeHeMolvig(impl, r,t):
    """Calculate 3He3He reactivity (Molvig-Knudsen)."""
    R = impl.rfuel(t)
    if r > R or impl.Ti(r,t) < 0.5:
        return 0
    return quad(integrand4, 0, 1000, args=(r,t,impl))[0]

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
    ret = 0 # burn rate
    ret2 = 0 # Ti 'rate'
    
    for r in numpy.arange( impl.rmin(t) , impl.rfuel(t) , dr ):
        r1 = r + dr/2
        f1 = fD(impl,r1,t)
        temp = svDDMolvig(impl,r1,t)*pow(impl.ni(r1,t),2)*(f1*f1/2)*4*math.pi*pow(r1,2)*dr
        ret += temp
        ret2 += temp*impl.Ti(r1,t)
    return [ret , ret2]
def D3Herate(impl, t):
    """Calculate the D3He burn rate at time t (s)."""
    f1 = 0 # D fraction (atomic)
    f2 = 0 # 3He fraction (atomic)
    ret = 0 # burn rate
    ret2 = 0 # Ti 'rate'

    for r in numpy.arange( impl.rmin(t) , impl.rfuel(t) , dr ):
        r1 = r + dr/2
        f1 = fD(impl,r1,t)
        f2 = f3He(impl,r1,t)
        temp = svD3HeMolvig(impl,r1,t)*pow(impl.ni(r1,t),2)*(f1*f2)*4*math.pi*pow(r1,2)*dr
        ret += temp
        ret2 += temp*impl.Ti(r1,t)
    return [ret , ret2]
def HeHerate(impl, t):
    """Calculate the 3He3He burn rate at time t (s)."""
    f1 = 0 # 3he fraction (atomic)
    ret = 0 # burn rate
    ret2 = 0 # Ti 'rate'

    for r in numpy.arange( impl.rmin(t) , impl.rfuel(t) , dr ):
        r1 = r + dr/2
        f1 = f3He(impl,r1,t)
        temp = svHeHeMolvig(impl,r1,t)*pow(impl.ni(r1,t),2)*(f1*f1/2)*4*math.pi*pow(r1,2)*dr
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
        
    # Do precomputation
    LeffInit(impl)
    NkInit(impl)
    
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
    rateFile = csv.writer(open(os.path.join(OutputDir,'BurnRate_Molvig.csv'),'w'))
    rateFile.writerow( ["t (s)", "DD (1/s)", "D3He (1/s)", "3He3He (1/s)"] )
    yieldFile = csv.writer(open(os.path.join(OutputDir,'Yield_Molvig.csv'),'w'))
    TiFile = csv.writer(open(os.path.join(OutputDir,'Ti_Molvig.csv'),'w'))
    BTFile = csv.writer(open(os.path.join(OutputDir,'BangTime_Molvig.csv'),'w'))
    
    #iterate over all time:
    for t in list(numpy.arange(impl.tmin(), impl.tmax(), dt)):
        #DD
        [dDD,dTiDD] = DDrate(impl,t)
        YDD += dDD*dt
        TiDD += dTiDD*dt
        if dDD > PeakRateDD:
            BTDD = t
            PeakRateDD = dDD
        #D3He
        [dD3He,dTiD3He] = D3Herate(impl,t)
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
        #output
        rateFile.writerow( [t, dDD, dD3He, d3He3He] )
        
    #If there is yield for a species, do output:
    if YDD > 0:
        TiDD = TiDD / YDD
        print("Molvig DD yield = " + '{:.2e}'.format(YDD))
        print("Molvig DD Ti = " + '{:.2f}'.format(TiDD))
        print("Molvig DD BT = " + '{:.2e}'.format(BTDD))
        yieldFile.writerow( ["DD",YDD] )
        TiFile.writerow( ["DD",TiDD] )
        BTFile.writerow( ["DD",BTDD] )
    if YD3He > 0:
        TiD3He = TiD3He / YD3He
        print("Molvig D3He yield = " + '{:.2e}'.format(YD3He))
        print("Molvig D3He Ti = " + '{:.2f}'.format(TiD3He))
        print("Molvig D3He BT = " + '{:.2e}'.format(BTD3He))
        yieldFile.writerow( ["D3He",YD3He] )
        TiFile.writerow( ["D3He",TiD3He] )
        BTFile.writerow( ["D3He",BTD3He] )
    if YHeHe > 0:
        TiHeHe = TiHeHe / YHeHe
        print("Molvig 3He3He yield = " + '{:.2e}'.format(YHeHe))
        print("Molvig 3He3He Ti = " + '{:.2f}'.format(TiHeHe))
        print("Molvig 3He3He BT = " + '{:.2e}'.format(BTHeHe))
        yieldFile.writerow( ["3He3He",YHeHe] )
        TiFile.writerow( ["3He3He",TiHeHe] )
        BTFile.writerow( ["3He3He",BTHeHe] )
        