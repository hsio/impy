# Python-based implosion analyzer
# A. Zylstra 2012/08/07

from datetime import *
import os, sys, inspect


print("----------------------------------------")
print("pyImplosion")
print("Author: Alex Zylstra")
print("Date: Aug 8, 2012")
print("v0.5.0")
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

# Have the user select and set up and implosion type
from Config import *
impl = implSelector()
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
for i in modules:
    t1 = datetime.now()
    print( moduleNames[runIndex] )
    runIndex += 1
    
    i.run(impl)
    t2 = datetime.now()
    print( '{:.1f}'.format((t2-t1).total_seconds()) + "s elapsed")
