# make various trajectories for plotting
# A. Zylstra 2013/04/13

from Implosion import *
from Resources.IO import *
from numpy import arange
import math
import csv
import os

zones = 10 # Lagrance for every 10th zone

# ------------------------------------
# Make Lagrange plots
# ------------------------------------
def LagrangePlots(impl):
    """Make Lagrange plots of initial material."""
    File = csv.writer(open(os.path.join(OutputDir,'Lagrange.csv'),'w'))
    #steps in radius:
    for ir in range(impl.ir_min() , impl.ir_max() , zones ):
        for it in range(impl.it_min(), impl.it_max() ):
            r = impl.r(ir,it)
            t = impl.t(it)
            File.writerow( [ t , r ] )
        File.writerow([])
        
# ------------------------------------
# Spit out shell radius
# ------------------------------------
def ShellTrajectory(impl):
    """Write the shell radius versus time."""
    File = csv.writer(open(os.path.join(OutputDir,'ShellR.csv'),'w'))
    File.writerow( [ "Time (ns)" , "Fuel Radius (cm)" ] )
    for it in range(impl.it_min(), impl.it_max() ):
        t = impl.t(it)
        r = impl.r( impl.ir_fuel() , it )
        File.writerow( [ t , r ] )

def run(impl):
    # input sanity check:
    if not isinstance(impl,Implosion):
        print("WARNING: invalid input.")
        return
        
    LagrangePlots(impl)
    ShellTrajectory(impl)