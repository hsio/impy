# Run speed tests on an implosion
# Author: A Zylstra
# Date: 2013/04/13

from Implosion import *
from datetime import *
import numpy.random as nprand

def run(impl):
    # input sanity check:
    if not isinstance(impl,Implosion):
        print("WARNING: invalid input.")
        return
    
    print('Speed tests')
    
    # getting min, max times
    print("Get min, max times (1E6)x")
    t1 = datetime.now()
    for i in range(1000000):
        temp = impl.it_min()
        temp = impl.it_max()
    t2 = datetime.now()
    print("    " + '{:.1f}'.format((t2-t1).total_seconds()) + "s")
    
    # getting min, max radii
    print("Get min, max radii (1E6)x at min, max times")
    t1 = datetime.now()
    tmin = impl.it_min()
    tmax = impl.it_max()
    for i in range(1000000):
        temp = impl.ir_min()
        temp = impl.ir_max()
    t2 = datetime.now()
    print("    " + '{:.1f}'.format((t2-t1).total_seconds()) + "s")
    
    #generate random r, random t
    print("Generating random (r,t) coordinates.")
    RandR = []
    RandT = []
    rmin = impl.ir_min()
    rmax = impl.ir_max()
    for i in range(1000000):
        RandT.append( nprand.randint(tmin, tmax))
        RandR.append( nprand.randint(rmin, rmax))
        
    #query ion density
    print("Get ni (1E6)x at random coordinates.")
    t1 = datetime.now()
    for i in range(len(RandT)):
        temp = impl.ni( RandR[i] , RandT[i] )
    t2 = datetime.now()
    print("    " + '{:.1f}'.format((t2-t1).total_seconds()) + "s")
    print("Get Ti (1E6)x at random coordinates.")
    t1 = datetime.now()
    for i in range(len(RandT)):
        temp = impl.Ti( RandR[i] , RandT[i] )
    t2 = datetime.now()
    print("    " + '{:.1f}'.format((t2-t1).total_seconds()) + "s")
    print("Get A,Z,F (1E6)x at random coordinates.")
    t1 = datetime.now()
    for i in range(len(RandT)):
        temp = impl.IonA( RandR[i] , RandT[i] )
        temp = impl.IonZ( RandR[i] , RandT[i] )
        temp = impl.IonF( RandR[i] , RandT[i] )
    t2 = datetime.now()
    print("    " + '{:.1f}'.format((t2-t1).total_seconds()) + "s")