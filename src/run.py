# Python-based Guderley imploding shock calculations
# A. Zylstra 2012/01/27

from Guderley import *
from Postproc import *

print("----------------------------------------")
print("Guderley simulation")
print("Author: Alex Zylstra")
print("Date: Jan 27, 2012")

print("----------------------------------------")
print("Solving Guderley ...")
g = Guderley('Guderley')
print("done!")
print("----------------------------------------")
print("Running post-processor ...")
p = Postproc('Postproc',g)
p.write(1.05e-9)
p.run()
p.PrintTrajectories()
