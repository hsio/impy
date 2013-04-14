# Calculate energy and entropy in gas
# A. Zylstra 2013/04/13

from Implosion import *
from Resources.IO import *
from Resources.Constants import *
from Resources.Plasma import PFermi
from numpy import arange
from numpy import dot
import math
import csv
import os

# ------------------------------------
# Helper methods
# ------------------------------------
def ThermalEnergy(impl,it):
    """Calculate the total thermal energy in the gas at time index it. Returns Joules."""
    Energy = 0.
    for ir in range(impl.ir_min(), impl.ir_fuel()):
        Vol = impl.vol(ir,it)
        #ion thermal energy
        Energy += 1.5*impl.ni(ir,it)*Vol*kB*impl.Ti(ir,it)*11600.*1000.
        #electron thermal energy
        Energy += 1.5*impl.ne(ir,it)*Vol*kB*impl.Te(ir,it)*11600.*1000.
    return Energy*1e-7
    
def alpha(impl,it):
    """Calculate mass-weighted ratio of hydro pressure to Fermi pressure at time index it."""
    TotalMass = 1e-15 # small but non-zero
    alpha = 0.
    for ir in range(impl.ir_min(), impl.ir_fuel()):
        Vol = impl.vol(ir,it)
        Mass = Vol*impl.ni(ir,it)*impl.Abar(ir,it)*mp
        Pf = max(PFermi(impl.ne(ir,it)) , 1e-9)
        alpha += (impl.P(ir,it) / Pf ) * Mass
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
    for it in range(impl.it_min(), impl.it_max()):
        tempE = ThermalEnergy(impl,it)
        tempA = alpha(impl,it)
        if tempE > maxE:
            maxE = tempE
        if tempA > maxA:
            maxA = tempA
        T.append( impl.t(it) )
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