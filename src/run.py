# Python-based Guderley imploding shock calculations
# A. Zylstra 2012/02/06

from Guderley import *
from Postproc import *
from datetime import *

print("----------------------------------------")
print("Guderley simulation")
print("Author: Alex Zylstra")
print("Date: Feb 6, 2012")

print("----------------------------------------")
t1 = datetime.now()
print("Solving Guderley ...")
g = Guderley('Guderley')
t2 = datetime.now()
print("done!")
print( '{:.1f}'.format((t2-t1).total_seconds()) + "s elapsed")

print("----------------------------------------")
print("Initializing post-processor ...")
t1 = datetime.now()
p = Postproc('Postproc',g)
t2 = datetime.now()
print("done!")
print( '{:.1f}'.format((t2-t1).total_seconds()) + "s elapsed")

t1 = datetime.now()
print("postproc calculations...")
p.write(1.05e-9)
p.run()
p.PrintTrajectories()
p.PrintEnergy()
t2 = datetime.now()
print("done!")
print( '{:.1f}'.format((t2-t1).total_seconds()) + "s elapsed")
