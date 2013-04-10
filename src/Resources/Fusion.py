# Fusion rate/xsection calculators
# A. Zylstra 2013/03/13

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
    # Angulo 1999
    if Ti <= 0.1:
        return 0
    T = T9(Ti)
    if T <= 0.11:
        return (1/Na)*1.81e3*pow(T,-2/3)*math.exp(-3.721*pow(T,-1/3))*(1+14.3*T-90.5*pow(T,2)+395*pow(T,3))
    return (1/Na)*2.58e3*pow(T,-2/3)*math.exp(-3.721*pow(T,-1/3))*(1+3.96*T+.116*pow(T,2))

def DD(Ti):
    """DD reactivity, Ti in keV."""
    # from Bosch and Hale
    C1 = 5.43360e-12
    C2 = 5.85778e-3
    C3 = 7.68222e-3
    C4 = 0
    C5 = -2.964e-6
    C6 = 0
    C7 = 0
    BG = 31.3970
    mc2 = 937814
    theta = Ti / (1-Ti*(C2+Ti*(C4+Ti*C6))/(1+Ti*(C3+Ti*(C5+Ti*C7))))
    xi = math.pow( BG*BG/(4*theta) , 0.333 )
    return C1*theta*math.sqrt(xi/(mc2*math.pow(Ti,3)))*math.exp(-3*xi)

def D3He(Ti):
    """D3He reactivity, Ti in keV."""
    # from Bosch and Hale
    C1 = 5.51036e-10
    C2 = 6.41918e-3
    C3 = -2.02896e-3
    C4 = -1.91080e-5
    C5 = 1.35776e-4
    C6 = 0
    C7 = 0
    BG = 68.7508
    mc2 = 1124572
    theta = Ti / (1-Ti*(C2+Ti*(C4+Ti*C6))/(1+Ti*(C3+Ti*(C5+Ti*C7))))
    xi = math.pow( BG*BG/(4*theta) , 0.333 )
    return C1*theta*math.sqrt(xi/(mc2*math.pow(Ti,3)))*math.exp(-3*xi)

def HeHe(Ti):
    """3He3He reactivity, Ti in keV."""
    # Angulo 1999
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
    theta = Ti / (1-Ti*(C2+Ti*(C4+Ti*C6))/(1+Ti*(C3+Ti*(C5+Ti*C7))))
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

def p11B(Ti):
    """p11B reactivity, Ti in keV."""
    # from Angulo
    T9 = Ti*11600*1000*1e-9
    ret = 1
    if T9 <= 2:
        ret = (2.68e12*math.pow(T9,-2/3))*math.exp(-12.097*math.pow(T9,-1/3))
        ret = ret * (1 + 1.62*T9 - 1.31*T9**2 + 0.26*T9**3)
        ret = ret + (2.12e6)*math.pow(T9,-3/2)*math.exp(-1.724/T9)
    else:
        ret = (5.84e11*math.pow(T9,-2/3)*math.exp(-12.097*math.pow(T9,-1/3)))
        ret = ret / ( (math.pow(T9,-2/3) - 1.47)**2 + 0.187 )
    return (1/Na)*ret

def p15N(Ti):
    """p15N reactivity, Ti in keV."""
    # from Angulo
    T9 = Ti*11600*1000*1e-9
    ret = 1
    if T9 > 2.5:
        ret = 4.17e7*math.pow(T9,0.917)*math.exp(-3.292/T9)
    else:
        ret = 1.12e12*math.pow(T9,-2/3)*math.exp(-15.253*math.pow(T9,-1/3) -(T9/.28)**2)*(1+4.95*T9+143*T9**2)
        ret = ret + 1.01e8*math.pow(T9,-3/2)*math.exp(-3.643/T9)
        ret = ret + 1.19e9*math.pow(T9,-3/2)*math.exp(-7.406/T9)
    return (1/Na)*ret

    
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
    if En <= 1:
        return 0
    # Bosch and Hale
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
    if En <= 1:
        return 0
    # Bosch and Hale
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
    if En <= 1:
        return 0
    #return (1.e-24) * (5000 / En) * math.exp(-31.29*4.0*math.sqrt(3/En))
    En = 2*En # convert energy from CM to lab = x (A1+A2)/A2
    # Angulo
    Sa = 5.18 #MeV b
    Sb = -2.22 #b
    Sc = 0.804 #b / MeV
    S = Sa + Sb*En*1e-3 + Sc*math.pow(En*1e-3,2) # MeV b
    return (1.e-24) * (2*S/(1e-3*En)) * math.exp(-31.29*4.0*math.sqrt(3/En))
def sigmaDT(En):
    """ DT cross section as a function of CM energy. [E] = keV, [sigma] = cm^2."""
    if En <= 1:
        return 0
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
    if En <= 1:
        return 0
    # S factor from S. Winkler et al, J. Phys G Nud. Pan. Phys. 18 (1992)
    En = 2*En # convert energy from CM to lab = x (A1+A2)/A2
    SE = 0.20 - 0.32*En*1e-3 + 0.476*math.pow( En*1e-3 , 2) # MeV b
    BG = 54.342
    return (1.e-24) * (SE/(1e-3*En)) * math.exp(-BG / math.sqrt(En))