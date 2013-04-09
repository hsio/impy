# Plasma parameter calculators
# A. Zylstra 2012/08/08

import math
from Resources.Constants import *

# units for inputs
# [T] = keV
# [n] = 1/cc

def IonMFP(ni,Z,A):
    """Calculate the ion mean free path (cm)."""
    if ni == 0:
        return 0
    return (2.28e7)*(1./Z)*math.sqrt( A/ni )
    
def nuii(ni, Ti, Z, A):
    """Calculate the ion-ion collision rate (s)."""
    if ni == 0:
        return 0
    return (4.80e-8)*pow(Z,4)*math.sqrt(1/(A*pow(1000*Ti,3)))*ni*LogL(ni,Ti,Z,A)
    
def Tauii(ni, Ti, Z, A):
    """Calculate the ion-ion collision time (s)."""
    if ni == 0:
        return 0
    nu = nuii(ni, Ti, Z, A)
    return (1./nu)
    
def LambdaD(ne, Te):
    """Calculate the Debye length (cm)."""
    if ne == 0:
        return 0
    return (7.43e2)*math.sqrt( 1000*Te / ( ne ) )
    
def uTherm(Ti, Z):
    """Calculate the ion thermal velocity (cm/s)."""
    ret = (9.79e5)*math.sqrt( 1000*Ti/Z )
    return max(1,ret) #causes problems if u = 0
    
def LogL(ni, Ti, Z, A): #See C.K. Li PRL 1993
    """Calculate the Coulomb logarithm."""
    if ni == 0:
        return 0
    mr = mp*A/2.
    pperp = pow(e*Z,2) / (mr * pow(uTherm(Ti,Z),2) )
    pmin = math.sqrt( pow(pperp,2) + pow(hbar/(2*mr*uTherm(Ti,Z)),2) )
    return math.log( LambdaD(ni*Z,Ti) / pmin )
    
def PFermi(ne):
    """Calculate the Fermi pressure for given electron density."""
    Pf = ( pow(h,2) / (20*me) )*pow(3/3.1415,2/3)*pow(ne,5/3) #dyn
    return Pf * 1e-6 * 1e-9 #Gbar