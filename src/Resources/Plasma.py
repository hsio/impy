# Plasma parameter calculators
# A. Zylstra 2012/09/21

import math
from Resources.Constants import *

# units for inputs
# [T] = keV
# [n] = 1/cc
# Z,A in atomic units

def PlasmaParameter(ne, Te):
    """Calculate plasma parameter."""
    return (1.72e9)*math.pow(1000*Te,3/2)/math.sqrt(ne)
    
def IonMFP(ni, Ti, Z, A):
    """Calculate the ion mean free path (cm)."""
    return uTherm(Ti,A)*Taui(ni,Ti,Z,A)
    
def Taui(ni, Ti, Z, A): # Atzeni 10.134
    """Calculate the ion-ion collision time (s)."""
    if ni == 0:
        return 0
    return (6.60e-10)*math.pow(Ti,1.5)*math.pow(A,0.5)/((ni/1e21)*math.pow(Z,4)*LogLi(ni,Ti,Z,A))
    
def Taue(ne, Te, Z): # Atzeni 10.134
    """Calculate electron collision time (s)."""
    if ne == 0:
        return 0
    return (1.09e-11)*math.pow(Te,1.5)/( (ne/(1e21*Z)) * math.pow(Z,2) * LogLe(ne,Te) )
    
def Tauei(ne, Te, Z): # Atzeni 10.135
    """Calculate electron-ion collision time (s)."""
    return (Z*mp/me)*Taue(ne,Te,Z)
    
def LambdaD(ne, Te):
    """Calculate the Debye length (cm)."""
    if ne == 0:
        return 0
    return (7.43e2)*math.sqrt( 1000*Te / ( ne ) )
    
def fpe(ne):
    """ Plasma frequency (Hz)."""
    return (5.64e4/(2*math.pi))*math.sqrt(ne)
    
def uTherm(Ti, A):
    """Calculate the ion thermal velocity (cm/s)."""
    ret = (9.79e5)*math.sqrt( 1000*Ti/A )
    return max(1,ret) #causes problems if u = 0
    
def LogLi(ni, Ti, Z, A): #Formulary
    """Calculate the Coulomb logarithm for ion-ion collisions."""
    if ni == 0 or Ti == 0:
        return 0
    # assuming one ion species or an average-ion Zbar
    return 33.36 - math.log(Z*Z) - 0.5*math.log(ni*Z*Z) + 1.5*math.log(Ti)
    
def LogLe(ne, Te): # Atzeni 10.136
    """Calculate the Coulomb logarithm for electron-ion collisions."""
    if ne == 0 or Te == 0:
        return 0
    if Te <= 0.01:
        print("Warning: LogLe called in invalid regime.")
    return 7.1 - 0.5*math.log(ne/1e21) + math.log(Te)
    
def PFermi(ne):
    """Calculate the Fermi pressure for given electron density."""
    Pf = ( pow(h,2) / (20*me) )*pow(3/3.1415,2/3)*pow(ne,5/3) #dyn
    return Pf * 1e-6 * 1e-9 #Gbar