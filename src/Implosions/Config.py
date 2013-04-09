# Add implosion types to this file so they can be selected
# To add options, follow instructions where comments are
# marked [MODIFY]
# A. Zylstra 2012/08/15

import sys
from Implosion import * #abstract class

# [MODIFY] add line:
# from <filename> import *
# where <filename> is the name of your file
from Guderley import *
from LILAC import *
from HYADES import *
# [MODIFY] add line:
# Implosion.register(<class name>)
Implosion.register(Guderley)
Implosion.register(LILAC)
Implosion.register(HYADES)

def implSelector():
    """Let thte user select implosion type from options."""
    # [MODIFY] add the name of your file's class to this list.
    implTypes = ['Guderley', 'LILAC', 'HYADES']

    # Do mode selection
    index = 0
    for i in implTypes:
        print(str(index) + " = " + i)
        index += 1
    mode = int(input("mode = "))
    
    print("Generate implosion...")
    impl = 0
    # [MODIFY] add your class here, in same order as implTypes array
    if mode == 0:
        impl = Guderley()
    if mode == 1:
        impl = LILAC()
    if mode == 2:
        impl = HYADES()
        
    return impl
    
def implAuto():
    """Automatic implosion generation from command line arguments."""
    #sanity check
    if len(sys.argv) < 3:
        print("ERROR: not enough input parameters")
        sys.exit()
    
    mode = float(sys.argv[2])
    impl = 0
    # [MODIFY] add your class here, in same order as implTypes array
    if mode == 0:
        impl = Guderley(sys.argv[3:])
    if mode == 1:
        impl = LILAC(sys.argv[3])
    if mode == 2:
        impl = HYADES(sys.argv[3])
    
    return impl