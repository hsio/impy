# Fusion rate/xsection calculators
# A. Zylstra 2012/09/05

#All reactivities use units cm^3/s
import math
from Resources.Constants import *
from Resources.Plasma import *
from scipy.integrate import quad

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
    
def DT(Ti):
    """DT reactivity, Ti in keV."""
    # from Bosch and Hale
    C1 = 1.17302e-9
    C2 = 1.51361e-2
    C3 = 7.51886e-2
    C4 = 4.60643e-3
    C5 = 1.35e-2
    C6 = -1.0675e-4
    C7 = 1.366e-5
    BG = 34.3827
    mc2 = 1124656
    theta = Ti / (1-Ti*(C2+Ti*(C4+Ti*C6))/(1+Ti*(C3+Ti*(C5*Ti*C7))))
    #print("----")
    #print(theta)
    #print(Ti)
    #print(BG*BG/(4*theta))
    xi = math.pow( BG*BG/(4*theta) , 0.333 )
    return C1*theta*math.sqrt(xi/(mc2*math.pow(Ti,3)))*math.exp(-3*xi)

def TTintegrand(En, Ti):
    """ for integrating TT cross section."""
    return c*math.sqrt(Ti/(math.pi*3e3*938))*math.pow(2/Ti,2)*sigmaTT(En)*math.exp(-En/Ti)*En
def TT(Ti):
    """TT reactivity, Ti in keV."""
    E1 = max(Ti/2, 0.5)
    E2 = max(10*Ti, 200)
    ret = quad(TTintegrand, E1, E2, args=(Ti), epsrel=1e-4, epsabs=0)
    return ret[0]
    
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
    En = 2*En # convert energy from CM to lab = x (A1+A2)/A2
    #Junker et al, 1998
    Sa = 5.3 #MeV b
    Sb = -3.7 #b
    Sc = 1.9 #b / MeV
    S = Sa + Sb*En*1e-3 + Sc*math.pow(En*1e-3,2) # MeV b
    return (1.e-24) * (S/(1e-3*En)) * math.exp(-31.29*4.0*math.sqrt(3/En))
def sigmaDT(En):
    """ DT cross section as a function of CM energy. [E] = keV, [sigma] = cm^2."""
    # Bosch and Hale
    A1 = 6.927e4
    A2 = 7.454e8
    A3 = 2.050e6
    A4 = 5.2002e4
    B1 = 6.38e1
    B2 = -9.95e-1
    B3 = 6.981e-5
    B4 = 1.728e-4
    SE = (A1 + En*(A2 + En*(A3 + En*A4))) / (1 + En*(B1+En*(B2+En*(B3+En*B4))))
    BG = 34.3827
    return (1.e-27) * SE / ( En*math.exp(BG/math.sqrt(En)) )
def sigmaTT(En):
    """ TT cross section as a function of CM energy. [E] = keV, [sigma] = cm^2."""
    # S factor from S. Winkler et al, J. Phys G Nud. Pan. Phys. 18 (1992)
    En = 2*En # convert energy from CM to lab = x (A1+A2)/A2
    SE = 0.20 - 0.32*En*1e-3 + 0.476*math.pow( En*1e-3 , 2) # MeV b
    BG = 54.342
    return (1.e-24) * (SE/(1e-3*En)) * math.exp(-BG / math.sqrt(En))