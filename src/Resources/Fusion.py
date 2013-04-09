# Fusion rate/xsection calculators
# A. Zylstra 2012/08/22

#All reactivities use units cm^3/s
import math
from Resources.Constants import *
from Resources.Plasma import *

def T9(Ti):
    """Convert from ion temperature in keV to T9."""
    return Ti*11600*1000*pow(10,-9)

# -----------------------------------------------
#      Maxwell-averaged Reactivities
# -----------------------------------------------
def HD(Ti):
    """HD reactivity, Ti in keV."""
    if Ti <= 0.1:
        return 0
    T = T9(Ti)
    if T <= 0.11:
        return (1/Na)*1.81e3*pow(T,-2/3)*math.exp(-3.721*pow(T,-1/3))*(1+14.3*T-90.5*pow(T,2)+395*pow(T,3))
    return (1/Na)*2.58e3*pow(T,-2/3)*math.exp(-3.721*pow(T,-1/3))*(1+3.96*T+.116*pow(T,2))

def DD(Ti):
    """DD reactivity, Ti in keV."""
    if Ti <= 0.1:
        return 0
    T = T9(Ti)
    return (1/Na)*4.67e8*pow(T,-2/3)*math.exp(-4.259*pow(T,-1/3))*(1+1.079*T-0.1124*pow(T,2)+5.68e-3*pow(T,3))

def D3He(Ti):
    """D3He reactivity, Ti in keV."""
    if Ti <= 0.1:
        return 0
    return 4.98e-16*math.exp(-0.152*pow(math.fabs(math.log(Ti/802.6)),2.65))

def HeHe(Ti):
    """3He3He reactivity, Ti in keV."""
    if Ti <= 0.1:
        return 0
    T = T9(Ti)
    return (1/Na)*5.59e10*pow(T,-2/3)*math.exp(-12.277*pow(T,-1/3))*(1-0.135*T+.0254*pow(T,2)-1.29e-3*pow(T,3))

# -----------------------------------------------
#      cross sections
# -----------------------------------------------
def sigmaDDn(En):
    """ DD cross section as a function of CM energy. [E] = keV, [sigma] = cm^2."""
    return (1.e-24) * (372 / (1+pow(1.220-4.36e-4*En,2))) / (En*(math.exp(46.097/math.sqrt(En)) - 1))
def sigmaD3He(En):
    """ D3He cross section as a function of CM energy. [E] = keV, [sigma] = cm^2."""
    return (1.e-24) * (647 + 50200/(1+pow(1.076 - 3.98e-2*En,2))) / (En*(math.exp(89.27/math.sqrt(En))-1))
def sigmaHeHe(En):
    """ 3He3He cross section as a function of CM energy. [E] = keV, [sigma] = cm^2."""
    return (1.e-24) * (5000 / En) * math.exp(-31.29*4.0*math.sqrt(3/En))