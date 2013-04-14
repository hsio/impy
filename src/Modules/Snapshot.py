# Write a snapshot of the implosion
# A. Zylstra 2013/04/13

from Implosion import *
from Resources.IO import *
from Resources.Constants import *
from Resources.Plasma import *
from numpy import arange
from numpy import dot
import math
import csv
import os

tsnap = 1e-9 # time to output snapshot

def run(impl):
    # input sanity check:
    if not isinstance(impl,Implosion):
        print("WARNING: invalid input.")
        return

    # find index closest to requested time:
    it_snap = 0
    delta = math.fabs(impl.t(it_snap) - tsnap)
    for it in range( impl.it_min() , impl.it_max() ):
        if math.fabs(impl.t(it) - tsnap) < delta:
            it_snap = it
            delta = math.fabs(impl.t(it_snap) - tsnap)
    
    File = csv.writer(open(os.path.join(OutputDir,'Snapshot.csv'),'w'))
    File.writerow( ["Snapshot time = ", tsnap] )
    File.writerow( ["r (cm)", "u (cm/s)", "cs (cm/s)", "rho (g/cc)", "ni (1/cc)", "Ti (keV)", "Te (keV)", "P (GBar)", "Abar", "Zbar", "Ion MFP (cm)", "tau_ii (s)"] )
    
    for ir in range(impl.ir_min(), impl.ir_max()):
        u = impl.u(ir,it_snap)
        cs = impl.c(ir,it_snap)
        rho = impl.rho(ir,it_snap)
        ni = impl.ni(ir,it_snap)
        Ti = impl.Ti(ir,it_snap)
        Te = impl.Te(ir,it_snap)
        P = impl.P(ir,it_snap)
        Zbar = impl.Zbar(ir,it_snap)
        Abar = impl.Abar(ir,it_snap)
        MFP = IonMFP(ni, Ti, Zbar , Abar)
        tau = Taui(ni, Ti, Zbar, Abar)
        r = impl.r(ir,it)
        File.writerow( [r, u, cs, rho, ni, Ti, Te, P, Abar, Zbar, MFP, tau] )
    