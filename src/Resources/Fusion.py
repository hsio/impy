# Fusion rate/xsection calculators
# A. Zylstra 2012/08/28

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
#      Gamow
# -----------------------------------------------
def Eg(Ti,Z1,Z2,A1,A2):
    """Gamow energy, in keV."""
    xi = 6.2696*pow(Z1*Z2,2/3)*pow(A1*A2/(A1+A2),1/3)
    return xi * pow(Ti,2/3)
    
# -----------------------------------------------
#      cross sections
#   from Bosch and Hale, Nucl. Fusion 32, 5 (1992)
# -----------------------------------------------
def S(En):
    """S-factor, using fit formalism in Bosch and Hale."""
    return 0
def sigmaDDn(En):
    """ DD cross section as a function of CM energy. [E] = keV, [sigma] = cm^2."""
    A1 = 5.3701e4
    A2 = 3.3027e2
    A3 = -1.2706e-1
    A4 = 2.9327e-5
    A5 = -2.5151e-9
    SE = (A1 + En*(A2 + En*(A3 + En*(A4 + En*A5))))
    BG = 31.3970
    return (1.e-27) * SE / ( En*math.exp(BG/math.sqrt(En)) )
def sigmaD3He(En):
    """ D3He cross section as a function of CM energy. [E] = keV, [sigma] = cm^2."""
    A1 = 5.7501e6
    A2 = 2.5226e3
    A3 = 4.5566e1
    B1 = -3.1995e-3
    B2 = -8.5530e-6
    B3 = 5.9014e-8
    SE = (A1 + En*(A2+En*A3)) / (1 + En*(B1+En*(B2+En*B3)))
    BG = 68.7508
    return (1.e-27) * SE / ( En*math.exp(BG/math.sqrt(En)) )
def sigmaHeHe(En):
    """ 3He3He cross section as a function of CM energy. [E] = keV, [sigma] = cm^2."""
    #return (1.e-24) * (5000 / En) * math.exp(-31.29*4.0*math.sqrt(3/En))
    #Adelberger et al, 2010
    Sa = 5.32 #MeV b
    Sb = -6.44 #b
    Sc = 30.7 #b / MeV
    S = Sa + Sb*En*1e-3 + Sc*math.pow(En*1e-3,2) # MeV b
    return (1.e-24) * (S/(1e-3*En)) * math.exp(-31.29*4.0*math.sqrt(3/En))