""" Plasma parameter calculators. Functions can take scalar quantities or numpy.ndarray
:author: Alex Zylstra
:date: 2014-01-03
"""

__author__ = 'Alex Zylstra'
__date__ = '2014-01-03'
__version__ = '0.1.0'

import numpy as np
from impy.resources.constants import *

# units for inputs
# [T] = keV
# [n] = 1/cc
# Z,A in atomic units

def PlasmaParameter(ne, Te):
    """Calculate plasma parameter.

    :param ne: Plasma electron number density [1/cc]
    :param Te: Plasma electron temperature [keV]
    """
    return (1.72e9)*np.power(1000*Te,3/2)/np.sqrt(ne)
    
def IonMFP(ni, Ti, Z, A):
    """Calculate the ion mean free path [cm]. Defined simply as thermal velocity times ion-ion collision time.

    :param ni: Ion number density [1/cc]
    :param Ti: Ion temperature [keV]
    :param Z: Ion atomic number [e]
    :param A: Ion mass [AMU]
    """
    return uTherm(Ti,A)*Taui(ni,Ti,Z,A)
    
def Taui(ni, Ti, Z, A):
    """Calculate the ion-ion collision time [s]. Uses Atzeni Eq 10.134.

    :param ni: Ion number density [1/cc]
    :param Ti: Ion temperature [keV]
    :param Z: Ion atomic number [e]
    :param A: Ion mass [AMU]
    """
    if ni <= 0:
        return 0
    return (6.60e-10)*np.power(Ti,1.5)*np.power(A,0.5)/((ni/1e21)*np.power(Z,4)*LogLi(ni,Ti,Z,A))
    
def Taue(ne, Te, Z):
    """Calculate electron collision time [s]. Uses Atzeni Eq 10.134.

    :param ni: Electron number density [1/cc]
    :param Ti: Electron temperature [keV]
    :param Z: Ion atomic number [e]
    """
    if ne == 0:
        return 0
    return (1.09e-11)*np.power(Te,1.5)/( (ne/(1e21*Z)) * np.power(Z,2) * LogLe(ne,Te) )
    
def Tauei(ne, Te, Z):
    """Calculate electron-ion collision time [s]. Uses Atzeni Eq 10.135.

    :param ni: Ion number density [1/cc]
    :param Ti: Ion temperature [keV]
    :param Z: Ion atomic number [e]
    :param A: Ion mass [AMU]
    """
    return (Z*mp/me)*Taue(ne,Te,Z)
    
def LambdaD(ne, Te):
    """Calculate the Debye length [cm].

    :param ne: Electron number density [1/cc]
    :param Te: Electron temperature [keV]
    """
    if ne <= 0:
        return 0
    return (7.43e2)*np.sqrt( 1000*Te / ( ne ) )
    
def fpe(ne):
    """ Plasma frequency [Hz].

    :param ne: Electron density [1/cc]
    """
    return (5.64e4/(2*np.pi))*np.sqrt(ne)
    
def uTherm(Ti, A):
    """Calculate the ion thermal velocity [cm/s].

    :param Ti: Ion temperature [keV]
    :param A: Ion mass [AMU]
    """
    ret = (9.79e5)*np.sqrt( 1000*Ti/A )
    return max(1,ret) #causes problems if u = 0
    
def LogLi(ni, Ti, Z, A):
    """Calculate the Coulomb logarithm for ion-ion collisions. Uses Atzeni Eq. 10.136.

    :param ni: Ion number density [1/cc]
    :param Ti: Ion temperature [keV]
    :param Z: Ion charge. Must use average ion for mixtures [e]
    :param A: Ion mass. Must use average ion for mixtures [AMU]
    """
    if ni == 0 or Ti == 0:
        return 0
    # assuming one ion species or an average-ion Zbar
    return 33.36 - np.log(Z*Z) - 0.5*np.log(ni*Z*Z) + 1.5*np.log(Ti)
    
def LogLe(ne, Te):
    """Calculate the Coulomb logarithm for electron collisions. Atzeni Eq. 10.136.

    :param ne: Ion number density [1/cc]
    :param Te: Ion temperature [keV]
    """
    if ne == 0 or Te == 0:
        return 0
    if Te <= 0.01:
        print("Warning: LogLe called in invalid regime.")
    return 7.1 - 0.5*np.log(ne/1e21) + np.log(Te)
    
def PFermi(ne):
    """Calculate the Fermi pressure in Gbar for given electron density.

    :param ne: Electron number density [1/cc]
    """
    Pf = ( pow(h,2) / (20*me) )*pow(3/3.1415,2/3)*pow(ne,5/3)  #dyn
    return Pf * 1e-6 * 1e-9  #Gbar