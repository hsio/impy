# Python-based Guderley imploding shock calculations
# A. Zylstra 2012/02/12

from Guderley import *
from LILAC import *
from Postproc import *
from datetime import *

print("----------------------------------------")
print("pyImplosion")
print("Author: Alex Zylstra")
print("Date: May 24, 2012")
print("----------------------------------------")
print("0 = Guderley")
print("1 = LILAC")
mode = float(input("mode = "))
print("----------------------------------------")
t1 = datetime.now()

if mode == 0:
    print("Solving Guderley ...")
    impl = Guderley('Guderley')

if mode == 1:
    print("Reading LILAC...")
    filename = input("LILAC file: ")
    impl = LILAC(filename)
    
print("done!")
t2 = datetime.now()
print( '{:.1f}'.format((t2-t1).total_seconds()) + "s elapsed")

print("----------------------------------------")
print("Initializing post-processor ...")
t1 = datetime.now()
p = Postproc('Postproc',impl)
t2 = datetime.now()
print("done!")
print( '{:.1f}'.format((t2-t1).total_seconds()) + "s elapsed")

t1 = datetime.now()
print("postproc calculations...")
p.write(0.5e-9)
p.run()
p.PrintTrajectories()
p.EnergyEntropy()
p.LagrangePlots()
t2 = datetime.now()
print("done!")
print( '{:.1f}'.format((t2-t1).total_seconds()) + "s elapsed")
