# Python-based implosion analyzer
# A. Zylstra 2012/08/28

# Command line arguments
# python3.2 run.py (y/n) [(mode) (mode opts)]
# (y/n): run all (or none) of the modules
# (mode): 0, 1, 2 for Guderley, LILAC, HYADES
# (mode opts): input parameters. Mode specific.

from datetime import *
import os, sys, inspect

print("----------------------------------------")
print("pyImplosion")
print("Author: Alex Zylstra")
print("Date: Aug 28, 2012")
print("v0.6.5")
print("----------------------------------------")

#path setup
# add Modules folder to path
cmd_subfolder = os.path.realpath(os.path.abspath(os.path.join(os.path.split(inspect.getfile( inspect.currentframe() ))[0],"Modules")))
if cmd_subfolder not in sys.path:
    sys.path.insert(0, cmd_subfolder)
# add Implosions folder to path
cmd_subfolder2 = os.path.realpath(os.path.abspath(os.path.join(os.path.split(inspect.getfile( inspect.currentframe() ))[0],"Implosions")))
if cmd_subfolder2 not in sys.path:
    sys.path.insert(0, cmd_subfolder2)

# Implosion setup. Can be done via CLI or argv
t1 = datetime.now()
from Config import *
impl = 0
if len(sys.argv) > 2:
    #automated setup
    print("Auto implosion generation...")
    impl = implAuto()
else:
    # Have the user select and set up and implosion type
    impl = implSelector()
    
t2 = datetime.now()
print( '{:.1f}'.format((t2-t1).total_seconds()) + "s elapsed")

print("----------------------------------------")

# Make a list of all modules in the 'Modules' folder
rawList = os.listdir("Modules")
moduleNames = []
for i in rawList:
    if ".py" in i:
        moduleNames.append( i[:-3] )

# import the modules to a map
modules = map(__import__, moduleNames)

# run all modules
runIndex = 0
runModule = 'n'
if len(sys.argv) < 2:
    sys.argv.append('n')
for i in modules:
    if not (sys.argv[1] == 'y' or sys.argv[1] == 'Y'):
        runModule = input("Run " + moduleNames[runIndex] + "? (y/n): ")
    else:
        print("Running " + moduleNames[runIndex] + "...")
    runIndex += 1
    t1 = datetime.now()
    
    if runModule == 'y' or runModule == 'Y' or sys.argv[1] == 'y' or sys.argv[1] == 'Y':
        i.run(impl)
        t2 = datetime.now()
        print( '{:.1f}'.format((t2-t1).total_seconds()) + "s elapsed")
