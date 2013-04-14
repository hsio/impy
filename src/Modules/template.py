# module template
# Author:
# Date:

from Implosion import *

def run(impl):
    # input sanity check:
    if not isinstance(impl,Implosion):
        print("WARNING: invalid input.")
        return
        
    print('foo')