# Calculate yields, Ti, BT using Molvig reduced reactivity
# A. Zylstra 2012/09/22

from Implosion import *
from Resources.IO import *
from Resources import Fusion
from Resources.Constants import *
import numpy
from scipy.integrate import quad
from scipy.optimize import fminbound
import scipy.interpolate
import math
import csv
import os

# integration step sizes
dt = 10e-12 #10ps
dr = 5e-4 #5um

# cutoffs in Ti for calculations
Ti_Min = 0.5
Ti_Max = 50
# for energy integrals of x sections
Emin = 0.5
Emax = 500

# reactions
reactions = []
# syntax: [ name , A1, Z1, A2, Z2, cross section fn , reactivity fn ]
reactions.append( [ "DD" , 2, 1, 2, 1, Fusion.sigmaDDn , Fusion.DD ] )
reactions.append( [ "DT" , 2, 1, 3, 1, Fusion.sigmaDT , Fusion.DT ] )
reactions.append( [ "TT" , 3, 1, 3, 1, Fusion.sigmaTT , Fusion.TT ] )
reactions.append( [ "D3He" , 2, 1, 3, 2, Fusion.sigmaD3He , Fusion.D3He ] )
reactions.append( [ "3He3He" , 3, 2, 3, 2, Fusion.sigmaHeHe , Fusion.HeHe ] )

# global implosion
impl = 0

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
def LeffInit():
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

def NkInit():
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
def Nk(r, t, Z1, Z2):
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
# Reactivity calculators
# ------------------------------------
def integrand(En, rxn, Nk, Ti):
    """Integrand for Molvig-Knudsen reactivity."""
    # see Brysk POP (1973)
    A = (rxn[1] + rxn[3] )/2
    return c*math.pow(Ti/(A*1e3*938),0.5)*math.pow(2/Ti,2)*rxn[5](En)*fK(En,Nk,Ti)*En
def integrandcm(En, rxn, Nk, Ti):
    """Integrand for Molvig-Knudsen reactivity CM energy."""
    return En*integrand(En, rxn, Nk, Ti)
def svMolvig(rxn, r,t):
    """Calculate reactivity (Molvig-Knudsen) for given reaction. Returns [sigmav , Ecm, Ti]"""
    R = impl.rfuel(t)
    Ti = impl.Ti(r,t)
    if (r > R) or (Ti < Ti_Min) or (Ti > Ti_Max):
        return [0,0,0,0]
    E1 = max(Ti/2, Emin)
    E2 = min(20*Ti, Emax)
    Nkx = Nk(r, t, rxn[2], rxn[4])
    Molvig = quad(integrand, E1, E2, args=(rxn,Nkx,Ti), epsrel = 1e-4, epsabs = 0)[0]
    if Molvig <= 0:
        return [ rxn[6](Ti) , Fusion.Eg(Ti,rxn[2],rxn[4],rxn[1],rxn[3]), Ti]
    Ecm = quad(integrandcm, E1, E2, args=(rxn,Nkx,Ti), epsrel = 1e-4, epsabs = 0)[0] / Molvig
    return [Molvig , Ecm, Ti, Nkx]
    
# ------------------------------------
# Burn rate calculators
# ------------------------------------
def rate(rxn, t):
    """Calculate the DD burn rate at time t (s)."""
    # fuel info
    A1 = rxn[1]
    Z1 = rxn[2]
    A2 = rxn[3]
    Z2 = rxn[4]
    dblcount = 1
    if (A1 == A2) and (Z1 == Z2):
        dblcount = 2 #account for factor of 2 if reactants are identical
    #fuel fractions
    f1 = 0
    f2 = 0
    # return values
    Yield = 0
    Ecm = 0 # Ecm 'rate'
    Ti = 0 # Ti 'rate'
    Nk = 0 # Nk 'rate'
    
    for r in numpy.arange( impl.rmin(t) , impl.rfuel(t) , dr ):
        r1 = r + dr/2
        f1 = f(r1,t,A1,Z1)
        f2 = f(r1,t,A2,Z2)
        temp = svMolvig(rxn,r,t)
        temp2 = pow(impl.ni(r1,t),2)*(f1*f2/dblcount)*4*math.pi*pow(r1,2)*dr
        Yield += temp[0]*temp2
        Ecm += temp[1]*temp[0]*temp2
        Ti += temp[2]*temp[0]*temp2
        Nk += temp[3]*temp[0]*temp2
    return [Yield , Ecm, Ti, Nk]

# ------------------------------------
# Test functions
# ------------------------------------
def checkCalc(verb=False):
    if verb:
        print("no verbosity implemented yet")
    
    # check that Molvig reactivities reduce to Maxwellians
    AccErr = 0.2 # accept 10% errors
    Ti = 10
    Nkx = 0
    for rxn in reactions:
        Molvig = quad(integrand, Ti/2, Ti*20, args=(rxn,Nkx,Ti), epsrel = 1e-4, epsabs = 0)[0]
        Thermal = rxn[6](Ti)
        if math.fabs(Molvig-Thermal)/Thermal >= AccErr:
            print("Failed!")
            print(rxn)
            print(Molvig)
            print(Thermal)
            print(rxn[5](20))
            return False
    return True
    
def calc_Molvig_Thermal(Nk):
    """ Calculate the ratio of Molvig to Thermal reactivity for given Nk."""
    file = csv.writer(open(os.path.join(OutputDir,'Molvig_Thermal.csv'),'w')) # Ratio vs Ti for fixed Nk
    file2 = csv.writer(open(os.path.join(OutputDir,'Molvig_Thermal2.csv'),'w')) # Ratio vs Ti for Nk normalized to D3He @ 15keV
    file3 = csv.writer(open(os.path.join(OutputDir,'Molvig_Thermal3.csv'),'w'))
    Tmin = 1
    Tmax = 30
    dT = 0.25
    # for second file output, normalized to conditions:
    Z1_0 = 1
    Z2_0 = 2 # D3He
    A1_0 = 2
    A2_0 = 3
    Ti_0 = 15
    # for third file output, normalized to conditions:
    rhoL = 1e-3 # 1 mg/cm2
    ni = 2e22 # 1/cc
    # construct header for rate file:
    file.writerow( ["Nk" , Nk] )
    file2.writerow( ["Nk" , Nk, "T0", Ti_0])
    file3.writerow( ["rhoL" , rhoL, "ni", ni] )
    rateFileHeader = ["Ti (keV)"]
    for i in reactions:
        rateFileHeader.append( i[0] )
    file.writerow( rateFileHeader )
    file2.writerow( rateFileHeader )
    file3.writerow( rateFileHeader )
    
    for Ti in numpy.arange(Tmin,Tmax,dT):
        line = [Ti] # line for output
        line2 = [Ti] # line for output
        line3 = [Ti] # line for output
        for rxn in reactions:
            # first file:
            Molvig = quad(integrand, Ti/2, Ti*20, args=(rxn,Nk,Ti), epsrel = 1e-4, epsabs = 0)[0]
            Thermal = rxn[6](Ti)
            line.append( Molvig/Thermal )
            
            # second file:
            Z1 = rxn[2]
            Z2 = rxn[4]
            A1 = rxn[1]
            A2 = rxn[3]
            NkNorm = Nk*math.pow(Ti/Ti_0,2)*math.pow(Z1_0*Z2_0/(Z1*Z2),2)*(A1+A2)/(A1_0+A2_0) # normalized Knudsen number
            Molvig = quad(integrand, Ti/2, Ti*20, args=(rxn,NkNorm,Ti), epsrel = 1e-4, epsabs = 0)[0]
            line2.append( Molvig/Thermal )
            
            # third file
            Zbar = (Z1+Z2)/2
            Abar = (A1+A2)/2
            ne = Zbar*ni
            Te = Ti
            lD = 743*math.sqrt(1000*Te / ne)
            b90 = 2*math.pow(Zbar*e,2)/(3*kB*1000*11600*Ti)
            rho = (ni/Na)*Abar
            L = rhoL / rho
            LL = max( math.log(lD/b90) , 1 )
            Nk = (math.sqrt(0.33)/math.pi)*pow(1000*11600*kB*Ti,2) / ((ne/Zbar)*LL*L*math.pow(e,4))
            Molvig = quad(integrand, Ti/2, Ti*20, args=(rxn,Nk,Ti), epsrel = 1e-4, epsabs = 0)[0]
            line3.append( Molvig/Thermal )
            
            
        file.writerow( line )
        file2.writerow( line2 )
        file3.writerow( line3 )
            
    
    
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
        
    # Do precomputation
    LeffInit()
    NkInit()
    # Quick self-check
    if not checkCalc():
        print("ERROR: Molvig reactivities do not properly reduce to Maxwellian. Aborting calculation.")
        return
    
    # Yields
    Y = numpy.zeros( len(reactions) )
    # ion temps (burn-averaged)
    Ti = numpy.zeros( len(reactions) )
    Ecm = numpy.zeros( len(reactions) )
    Nk = numpy.zeros( len(reactions) )
    # Bang (peak emission) times
    BT = numpy.zeros( len(reactions) )
    PeakRate = numpy.zeros( len(reactions) )
    
    # output files
    rateFile = csv.writer(open(os.path.join(OutputDir,'BurnRate_Molvig.csv'),'w'))
    yieldFile = csv.writer(open(os.path.join(OutputDir,'Yield_Molvig.csv'),'w'))
    EcmFile = csv.writer(open(os.path.join(OutputDir,'Ecm_Molvig.csv'),'w'))
    TiFile = csv.writer(open(os.path.join(OutputDir,'Ti_Molvig.csv'),'w'))
    NkFile = csv.writer(open(os.path.join(OutputDir,'Nk_Molvig.csv'),'w'))
    BTFile = csv.writer(open(os.path.join(OutputDir,'BangTime_Molvig.csv'),'w'))
    # construct header for rate file:
    rateFileHeader = ["t (s)"]
    for i in reactions:
        rateFileHeader.append( i[0] + " (1/s)" )
    rateFile.writerow( rateFileHeader )
    
    #iterate over all time:
    for t in list(numpy.arange(impl.tmin(), impl.tmax(), dt)):
        rateFileRow = [t]
        
        #iterate over reactions
        for i in range(len(reactions)):
            [dY, dEcm, dTi, dNk] = rate(reactions[i],t)
            Y[i] += dY*dt
            Ecm[i] += dEcm*dt
            Ti[i] += dTi*dt
            Nk[i] += dNk*dt
            rateFileRow.append(dY)
            if dY > PeakRate[i]:
                BT[i] = t
                PeakRate[i] = dY

        #output
        rateFile.writerow( rateFileRow )
        
    #If there is yield for a species, do output:
    #iterate over reactions
    for i in range(len(reactions)):
        if Y[i] > 0:
            Ti[i] = Ti[i] / Y[i]
            Nk[i] = Nk[i] / Y[i]
            Ecm[i] = Ecm[i] / Y[i]
            print(reactions[i][0]+" yield = " + '{:.2e}'.format(Y[i]))
            print(reactions[i][0]+" Ecm = " + '{:.2f}'.format(Ecm[i]))
            print(reactions[i][0]+" Nk = " + '{:.2f}'.format(Nk[i]))
            print(reactions[i][0]+" Ti = " + '{:.2f}'.format(Ti[i]))
            print(reactions[i][0]+" BT = " + '{:.2e}'.format(BT[i]))
            yieldFile.writerow( [reactions[i][0],Y[i]] )
            EcmFile.writerow( [reactions[i][0],Ecm[i]] )
            NkFile.writerow( [reactions[i][0],Nk[i]] )
            TiFile.writerow( [reactions[i][0],Ti[i]] )
            BTFile.writerow( [reactions[i][0],BT[i]] )
        