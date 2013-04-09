# Calculate energy and entropy in gas
# A. Zylstra 2012/08/15

from Implosion import *
from Resources.IO import *
from Resources.Constants import *
from Resources.Plasma import PFermi
from numpy import arange
from numpy import dot
import math
import csv
import os

# integration step sizes
dt = 10e-12 #10ps
dr = 5e-4 #5um

# ------------------------------------
# Helper methods
# ------------------------------------
def ThermalEnergy(impl,t):
    """Calculate the total thermal energy in the gas at time t (s). Returns Joules."""
    Energy = 0.
    for r in list(arange(impl.rmin(t), impl.rfuel(t), dr)):
        Vol = 4*math.pi*pow(r+dr/2,2)*dr
        #ion thermal energy
        Energy += 1.5*impl.ni(r,t)*Vol*kB*impl.Ti(r,t)*11600.*1000.
        #electron thermal energy
        Energy += 1.5*impl.ne(r,t)*Vol*kB*impl.Te(r,t)*11600.*1000.
    return Energy*1e-7
    
def alpha(impl,t):
    """Calculate mass-weighted ratio of hydro pressure to Fermi pressure at time t."""
    TotalMass = 1e-12 # small but non-zero
    alpha = 0.
    for r in list(arange(impl.rmin(t), impl.rfuel(t), dr)):
        Vol = 4*math.pi*pow(r+dr/2,2)*dr
        Mass = Vol*impl.ni(r,t)*impl.Abar(r,t)*mp
        Pf = max(PFermi(impl.ne(r,t)) , 1e-9)
        alpha += (impl.P(r,t) / Pf ) * Mass
        TotalMass += Mass
    alpha = alpha / TotalMass
    return alpha

def run(impl):
    # input sanity check:
    if not isinstance(impl,Implosion):
        print("WARNING: invalid input.")
        return
    
    File = csv.writer(open(os.path.join(OutputDir,'EnergyEntropy.csv'),'w'))
    E = []
    A = []
    T = []
    maxE = 0. #max energy dumped into gas after shock collapse
    maxA = 0. #max alpha in gas after shock collapse
    for t in list(arange(impl.tmin(), impl.tmax(), dt)):
        tempE = ThermalEnergy(impl,t)
        tempA = alpha(impl,t)
        if tempE > maxE:
            maxE = tempE
        if tempA > maxA:
            maxA = tempA
        T.append(t)
        E.append(tempE)
        A.append(tempA)
    
    # Summary output
    File.writerow( ["Thermal energy in gas after shock = " + str(maxE) + "J"] )
    File.writerow( ["Max gas P / Pf = " + str(maxA)] )
    File.writerow( [] )
    
    # Do output
    File.writerow( ["t (s)", "E (J)", "alpha"] )
    for i in range(len(E)):
        File.writerow( [ T[i] , E[i] , A[i] ] )