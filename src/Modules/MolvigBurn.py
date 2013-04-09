# Calculate yields, Ti, BT using Molvig reduced reactivit
# A. Zylstra 2012/08/30

from Implosion import *
from Resources.IO import *
from Resources.Fusion import *
from Resources.Constants import *
import numpy
from scipy.integrate import quad
from scipy.optimize import fminbound
import scipy.interpolate
import math
import csv
import os

# integration step sizes
dt = 20e-12 #10ps
dr = 10e-4 #5um
# for energy integrals of x sections
Emin = 0.5
Emax = 500

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
def Leff(R, r):
    """Effective scale length in spherical geometry. R is shell radius."""
    if r >= R:
        return 0
    return R*LeffInt(r/R)
    
# Knudsen number:
def CalcNk(Ti, Te, ne, Zbar, R, r):
    lD = 743*math.sqrt(1000*Te / ne)
    b90 = 2*math.pow(Zbar*e,2)/(3*kB*1000*11600*Ti)
    LL = max( math.log(lD/b90) , 1 )
    L = Leff(R,r)
    if L == 0:
        return 1e9
    Nk = (math.sqrt(0.33)/math.pi)*pow(1000*11600*kB*Ti,2) / ((ne/Zbar)*LL*L*math.pow(e,4))
    if Nk < 0 or L < 0:
        return 0.
    return Nk

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
            temp = CalcNk( impl.Ti(r,t) , impl.Te(r,t) , impl.ne(r,t) , impl.Zbar(r,t) , impl.rfuel(t) , r)
            NkArr[i,j] = temp 
def Nk(impl, r, t, Z1, Z2):
    """Knudsen number for fusion with ion charges Z1,Z2."""
    ir = int( r / dr )
    it = int( ( t - impl.tmin() ) / dt )
    return  ( NkArr[ir,it] / (pow(Z1*Z2,2)) )
    
# Knudsen distribution function:
def fK(En, Nk, Ti):
    """Knudsen distribution function."""
    #normalized energy
    eps = En / Ti
    #evaluated as prefactor and exponential
    #p1 = 1 / math.sqrt( math.pi + Nk*math.pow(eps,1.5) )
    #p2 = math.exp( -1.0*(eps+0.4*Nk*math.pow(eps,2.5)) )
    p1 = 1 / math.sqrt( math.pi + Nk*math.pow(eps,1.5) )
    p2 = math.exp( -1.*(eps+0.8*Nk*math.pow(eps,2.5)+0.32*math.pow(Nk*eps*eps,2))/(1+0.8*Nk*math.pow(eps,1.5)) )
    return p1*p2
    
# DD reactivity:
def integrandDD(En, NkDD, Ti):
    """Integrand for Molvig-Knudsen DD reactivity."""
    # see Brysk POP (1973)
    return math.sqrt((2+2)/4)*c*math.pow(Ti/(2*1e3*938),0.5)*math.pow(2/Ti,2)*sigmaDDn(En)*fK(En,NkDD,Ti)*En
def integrandDDcm(En, NkDD, Ti):
    """Integrand for Molvig-Knudsen DD reactivity CM energy."""
    return En*integrandDD(En, NkDD, Ti)
def svDDMolvig(impl, r,t):
    """Calculate DD reactivity (Molvig-Knudsen). Returns [sigmav , Ecm]"""
    R = impl.rfuel(t)
    Ti = impl.Ti(r,t)
    if r > R or Ti < 0.5:
        return [0,0]
    E1 = max(Ti/2, Emin)
    E2 = min(10*Ti, Emax)
    NkDD = Nk(impl, r, t, 1, 1)
    Molvig = quad(integrandDD, E1, E2, args=(NkDD,Ti), epsrel = 1e-4, epsabs = 0)[0]
    if Molvig <= 0:
        return [ DD(Ti) , Eg(Ti,1,1,2,2)]
    Ecm = quad(integrandDDcm, E1, E2, args=(NkDD,Ti), epsrel = 1e-4, epsabs = 0)[0] / Molvig
    return [Molvig , Ecm]
# D3He reactivity:
def integrandD3He(En, NkD3He, Ti):
    """Integrand for Molvig-Knudsen D3He reactivity."""
    # see Brysk POP (1973)
    return c*math.pow(Ti/(2*1.25e3*938),0.5)*math.pow(2/Ti,2)*sigmaD3He(En)*fK(En,NkD3He,Ti)*En
def integrandD3Hecm(En, NkD3He, Ti):
    """Integrand for Molvig-Knudsen D3He reactivity CM energy."""
    return En*integrandD3He(En, NkD3He, Ti)
def svD3HeMolvig(impl, r,t):
    """Calculate D3He reactivity (Molvig-Knudsen). Returns [sigmav , Ecm]"""
    R = impl.rfuel(t)
    Ti = impl.Ti(r,t)
    if r > R or Ti < 0.5:
        return [0,0]
    E1 = max(Ti/2, Emin)
    E2 = min(15*Ti, Emax)
    NkD3He = Nk(impl, r, t, 1, 2)
    Molvig = quad(integrandD3He, E1, E2, args=(NkD3He,Ti), epsrel = 1e-4, epsabs = 0)[0]
    if Molvig <= 0:
        return [ D3He(Ti) , Eg(Ti,1,2,2,3)]
    Ecm = quad(integrandD3Hecm, E1, E2, args=(NkD3He,Ti), epsrel = 1e-4, epsabs = 0)[0] / Molvig
    return [Molvig , Ecm]
# 3He3He reactivity:
def integrandHeHe(En, NkHeHe, Ti):
    """Integrand for Molvig-Knudsen 3He3He reactivity."""
    return c*math.pow(Ti/(2*1.5e3*938),0.5)*math.pow(2/Ti,2)*sigmaHeHe(2*En)*fK(En,NkHeHe,Ti)*En
def integrandHeHecm(En, NkHeHe, Ti):
    """Integrand for Molvig-Knudsen 3He3He reactivity."""
    return En*integrandHeHe(En, NkHeHe, Ti)
def svHeHeMolvig(impl, r,t):
    """Calculate 3He3He reactivity (Molvig-Knudsen). Returns [sigmav , Ecm]"""
    R = impl.rfuel(t)
    Ti = impl.Ti(r,t)
    if r > R or Ti < 0.5:
        return [0,0]
    E1 = max(Ti/2, Emin)
    E2 = min(20*Ti, Emax)
    NkHeHe = Nk(impl, r, t, 2, 2)
    Molvig = quad(integrandHeHe, E1, E2, args=(NkHeHe,Ti), epsrel = 1e-4, epsabs = 0)[0]
    if Molvig <= 0:
        return [ HeHe(Ti) , Eg(Ti,2,2,3,3)]
    Ecm = quad(integrandHeHecm, E1, E2, args=(NkHeHe,Ti), epsrel = 1e-4, epsabs = 0)[0] / Molvig
    return [Molvig , Ecm]

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
        temp = svDDMolvig(impl,r1,t)
        temp2 = pow(impl.ni(r1,t),2)*(f1*f1/2)*4*math.pi*pow(r1,2)*dr
        ret += temp[0]*temp2
        ret2 += temp[1]*temp[0]*temp2
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
        temp = svD3HeMolvig(impl,r1,t)
        temp2 = pow(impl.ni(r1,t),2)*(f1*f2)*4*math.pi*pow(r1,2)*dr
        ret += temp[0]*temp2
        ret2 += temp[1]*temp[0]*temp2
    return [ret , ret2]
def HeHerate(impl, t):
    """Calculate the 3He3He burn rate at time t (s)."""
    f1 = 0 # 3he fraction (atomic)
    ret = 0 # burn rate
    ret2 = 0 # Ti 'rate'

    for r in numpy.arange( impl.rmin(t) , impl.rfuel(t) , dr ):
        r1 = r + dr/2
        f1 = f3He(impl,r1,t)
        temp = svHeHeMolvig(impl,r1,t)
        temp2 = pow(impl.ni(r1,t),2)*(f1*f1/2)*4*math.pi*pow(r1,2)*dr
        ret += temp[0]*temp2
        ret2 += temp[1]*temp[0]*temp2
    return [ret , ret2]

# test Molvig reduced reactivity
def checkCalc(verb=False):
    if verb:
        print("9keV plasma, 6g/cc D, 10um scale length")
        Ti = 9
        ni = 6 / (2*mp) # 6 g/cc D plasma
        L = 10e-4 # 10um scale length
        Nk = CalcNk(Ti, Ti, ni, 1, L, 0)
        print(Nk)
        Molvig = quad(integrandDD, Ti/2, 10*Ti, args=(Nk,Ti), epsrel = 1e-4, epsabs = 0)[0]
        Thermal = DD( Ti )
        print(Molvig)
        print(Thermal)
        ratio = Molvig / Thermal
        print(ratio)
    # check that Molvig reactivities reduce to Maxwellians
    Ti = 10
    Nk = 0
    MolvigDD = quad(integrandDD, Ti/2, 10*Ti, args=(Nk,Ti), epsrel = 1e-4, epsabs = 0)[0]
    MolvigD3He = quad(integrandD3He, Ti/2, 15*Ti, args=(Nk,Ti), epsrel = 1e-4, epsabs = 0)[0]
    MolvigHeHe = quad(integrandHeHe, Ti/2, 20*Ti, args=(Nk,Ti), epsrel = 1e-4, epsabs = 0)[0]
    AccErr = 0.1 # accept 10% errors
    if (MolvigDD - DD(Ti))/MolvigDD > AccErr:
        return False
    if (MolvigD3He - D3He(Ti))/MolvigD3He > AccErr:
        return False
    if (MolvigHeHe - HeHe(Ti))/MolvigHeHe > AccErr:
        return False
    return True
    
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
    # Quick self-check
    if not checkCalc():
        print("ERROR: Molvig reactivities do not properly reduce to Maxwellian. Aborting calculation.")
        return
    
    # Yields
    YDD = 0
    YD3He = 0
    YHeHe = 0
    # ion temps (burn-averaged)
    EcmDD = 0
    EcmD3He = 0
    EcmHeHe = 0
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
    EcmFile = csv.writer(open(os.path.join(OutputDir,'Ecm_Molvig.csv'),'w'))
    BTFile = csv.writer(open(os.path.join(OutputDir,'BangTime_Molvig.csv'),'w'))
    
    #iterate over all time:
    for t in list(numpy.arange(impl.tmin(), impl.tmax(), dt)):
        #DD
        [dDD,dEcmDD] = DDrate(impl,t)
        YDD += dDD*dt
        EcmDD += dEcmDD*dt
        if dDD > PeakRateDD:
            BTDD = t
            PeakRateDD = dDD
        #D3He
        [dD3He,dEcmD3He] = D3Herate(impl,t)
        YD3He += dD3He*dt
        EcmD3He += dEcmD3He*dt
        if dD3He > PeakRateD3He:
            BTD3He = t
            PeakRateD3He = dD3He
        #3He3He
        [d3He3He, dEcm3He3He] = HeHerate(impl,t)
        YHeHe += d3He3He*dt
        EcmHeHe += dEcm3He3He*dt
        if d3He3He > PeakRateHeHe:
            BTHeHe = t
            PeakRateHeHe = d3He3He
        #output
        rateFile.writerow( [t, dDD, dD3He, d3He3He] )
        
    #If there is yield for a species, do output:
    if YDD > 0:
        EcmDD = EcmDD / YDD
        print("Molvig DD yield = " + '{:.2e}'.format(YDD))
        print("Molvig DD Ecm = " + '{:.2f}'.format(EcmDD))
        print("Molvig DD BT = " + '{:.2e}'.format(BTDD))
        yieldFile.writerow( ["DD",YDD] )
        EcmFile.writerow( ["DD",EcmDD] )
        BTFile.writerow( ["DD",BTDD] )
    if YD3He > 0:
        EcmD3He = EcmD3He / YD3He
        print("Molvig D3He yield = " + '{:.2e}'.format(YD3He))
        print("Molvig D3He Ecm = " + '{:.2f}'.format(EcmD3He))
        print("Molvig D3He BT = " + '{:.2e}'.format(BTD3He))
        yieldFile.writerow( ["D3He",YD3He] )
        EcmFile.writerow( ["D3He",EcmD3He] )
        BTFile.writerow( ["D3He",BTD3He] )
    if YHeHe > 0:
        EcmHeHe = EcmHeHe / YHeHe
        print("Molvig 3He3He yield = " + '{:.2e}'.format(YHeHe))
        print("Molvig 3He3He Ecm = " + '{:.2f}'.format(EcmHeHe))
        print("Molvig 3He3He BT = " + '{:.2e}'.format(BTHeHe))
        yieldFile.writerow( ["3He3He",YHeHe] )
        EcmFile.writerow( ["3He3He",EcmHeHe] )
        BTFile.writerow( ["3He3He",BTHeHe] )
        