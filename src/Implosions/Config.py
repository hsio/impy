# Add implosion types to this file so they can be selected
# To add options, follow instructions where comments are
# marked [MODIFY]
# A. Zylstra 2012/08/07

from datetime import *
from Implosion import * #abstract class

# [MODIFY] add line:
# from <filename> import *
# where <filename> is the name of your file
from Guderley import *
from LILAC import *
# [MODIFY] add line:
# Implosion.register(<class name>)
Implosion.register(Guderley)
Implosion.register(LILAC)

def implSelector():
    """Let thte user select implosion type from options."""
    # [MODIFY] add the name of your file's class to this list.
    implTypes = ['Guderley', 'LILAC']

    # Do mode selection
    index = 0
    for i in implTypes:
        print(str(index) + " = " + i)
        index += 1
    mode = int(input("mode = "))
    
    t1 = datetime.now()
    print("Generate implosion...")
    impl = 0
    # [MODIFY] add your class here, in same order as implTypes array
    if mode == 0:
        impl = Guderley()
    if mode == 1:
        impl = LILAC()
    
    t2 = datetime.now()
    print( '{:.1f}'.format((t2-t1).total_seconds()) + "s elapsed")
    
    return impl