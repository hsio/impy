# Write a snapshot of the implosion
# A. Zylstra 2012/08/15

from Implosion import *
from Resources.IO import *
from Resources.Constants import *
from Resources.Plasma import *
from numpy import arange
from numpy import dot
import math
import csv
import os

tsnap = 1e-9
dr = 5e-4 #5um

def run(impl):
    # input sanity check:
    if not isinstance(impl,Implosion):
        print("WARNING: invalid input.")
        return
    if (tsnap < impl.tmin()) or (tsnap > impl.tmax()):
        print("WARNING: invalid snapshot time.")
        return
    
    File = csv.writer(open(os.path.join(OutputDir,'Snapshot.csv'),'w'))
    File.writerow( ["Snapshot time = ", tsnap] )
    File.writerow( ["r (cm)", "u (cm/s)", "cs (cm/s)", "rho (g/cc)", "ni (1/cc)", "Ti (keV)", "Te (keV)", "P (GBar)", "Abar", "Zbar", "Ion MFP (cm)", "tau_ii (s)"] )
    
    for r in list(arange(impl.rmin(tsnap), impl.rmax(tsnap), dr)):
        u = impl.u(r,tsnap)
        cs = impl.c(r,tsnap)
        rho = impl.rho(r,tsnap)
        ni = impl.ni(r,tsnap)
        Ti = impl.Ti(r,tsnap)
        Te = impl.Te(r,tsnap)
        P = impl.P(r,tsnap)
        Zbar = impl.Zbar(r,tsnap)
        Abar = impl.Abar(r,tsnap)
        MFP = IonMFP(ni, Ti, Zbar , Abar)
        tau = Taui(ni, Ti, Zbar, Abar)
        File.writerow( [r, u, cs, rho, ni, Ti, Te, P, Abar, Zbar, MFP, tau] )
    