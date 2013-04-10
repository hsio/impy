# make various trajectories for plotting
# A. Zylstra 2012/10/26

from Implosion import *
from Resources.IO import *
from numpy import arange
import math
import csv
import os

# integration step sizes
dt = 10e-12 #10ps
dr = 25e-4 #25um

# ------------------------------------
# Make Lagrange plots
# ------------------------------------
def LagrangePlots(impl):
    """Make Lagrange plots of initial material."""
    File = csv.writer(open(os.path.join(OutputDir,'Lagrange.csv'),'w'))
    #steps in radius:
    for r in list(arange(impl.rmin(impl.tmin()),impl.rmax(impl.tmin()), dr)):
        #steps in time:
        rLast = r
        vLast = 0
        for t in list(arange(impl.tmin(), impl.tmax(), dt)):
            File.writerow( [ t , rLast+vLast*dt ] )
            rLast = rLast+vLast*dt
            vLast = impl.u(rLast, t)
        File.writerow([])
        
# ------------------------------------
# Spit out shell radius
# ------------------------------------
def ShellTrajectory(impl):
    """Write the shell radius versus time."""
    File = csv.writer(open(os.path.join(OutputDir,'ShellR.csv'),'w'))
    File.writerow( [ "Time (ns)" , "Fuel Radius (cm)" ] )
    for t in list(arange(impl.tmin(), impl.tmax(), dt)):
        File.writerow( [ t , impl.rfuel(t) ] )

def run(impl):
    # input sanity check:
    if not isinstance(impl,Implosion):
        print("WARNING: invalid input.")
        return
        
    LagrangePlots(impl)
    ShellTrajectory(impl)