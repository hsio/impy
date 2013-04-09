# Run speed tests on an implosion
# Author: A Zylstra
# Date: 8/14/12

from Implosion import *
from datetime import *
import random

def run(impl):
    # input sanity check:
    if not isinstance(impl,Implosion):
        print("WARNING: invalid input.")
        return
        
    #seen RNG
    random.seed()
    
    print('Speed tests')
    
    # getting min, max times
    print("Get min, max times (1E5)x")
    t1 = datetime.now()
    for i in range(100000):
        temp = impl.tmin()
        temp = impl.tmax()
    t2 = datetime.now()
    print("    " + '{:.1f}'.format((t2-t1).total_seconds()) + "s")
    
    # getting min, max radii
    print("Get min, max radii (1E5)x at min, max times")
    t1 = datetime.now()
    tmin = impl.tmin()
    tmax = impl.tmax()
    for i in range(50000):
        temp = impl.rmin(tmin)
        temp = impl.rmin(tmax)
        temp = impl.rmax(tmin)
        temp = impl.rmax(tmax)
    t2 = datetime.now()
    print("    " + '{:.1f}'.format((t2-t1).total_seconds()) + "s")
    
    #generate random r, random t
    print("Generating random (r,t) coordinates.")
    RandR = []
    RandT = []
    for i in range(100000):
        RandT.append( random.uniform(tmin, tmax))
        rmin = impl.rmin(RandT[len(RandT)-1])
        rmax = impl.rmax(RandT[len(RandT)-1])
        RandR.append( random.uniform(rmin, rmax))
        
    #query ion density
    print("Get ni (1E5)x at random coordinates.")
    t1 = datetime.now()
    for i in range(len(RandT)):
        temp = impl.ni( RandR[i] , RandT[i] )
    t2 = datetime.now()
    print("    " + '{:.1f}'.format((t2-t1).total_seconds()) + "s")
    print("Get Ti (1E5)x at random coordinates.")
    t1 = datetime.now()
    for i in range(len(RandT)):
        temp = impl.Ti( RandR[i] , RandT[i] )
    t2 = datetime.now()
    print("    " + '{:.1f}'.format((t2-t1).total_seconds()) + "s")
    print("Get A,Z,F (1E5)x at random coordinates.")
    t1 = datetime.now()
    for i in range(len(RandT)):
        temp = impl.IonA( RandR[i] , RandT[i] )
        temp = impl.IonZ( RandR[i] , RandT[i] )
        temp = impl.IonF( RandR[i] , RandT[i] )
    t2 = datetime.now()
    print("    " + '{:.1f}'.format((t2-t1).total_seconds()) + "s")