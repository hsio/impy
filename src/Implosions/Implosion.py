# Python-based abstract representation of an implosion
# All implosion types must implement these methods
# A. Zylstra 2013/04/13

from abc import ABCMeta

## Implosion meta class. Define methods that must be implemented
# in any type of implosion.
class Implosion(metaclass=ABCMeta):
    """Abstract class for implementing implosions."""
    name = ""
    
    ## Constructor
    # @param inName name of the Implosion object
    def __init__(self, inName):
        self.name = inName

    # required hydro variables
    def ni(self, ir, it):
        """Ion number density ni(r,t) with ir and it the radial and temporal indices. Returns ni in 1/cm^3."""
        return 0
    def ne(self, ir, it):
        """Electron number density ne(r,t) with ir and it the radial and temporal indices. Returns ne in 1/cm^3."""
        return 0
    def Ti(self, ir, it):
        """Ion temperature with ir and it the radial and temporal indices. Returns Ti in keV."""
        return 0
    def Te(self, ir, it):
        """Electron temperature with ir and it the radial and temporal indices. Returns Te in keV."""
        return 0
    def u(self, ir, it):
        """Fluid velocity u(r,t) with ir and it the radial and temporal indices. Returns u in um/ns."""
        return 0
    def c(self, ir, it):
        """Sound speed c(r,t) with ir and it the radial and temporal indices. Returns c in um/ns."""
        return 0
    def rho(self, ir, it):
        """Density rho(r,t) with ir and it the radial and temporal indices. Returns rho in g/cm3."""
        return 0
    def P(self, ir, it):
        """Pressure P(r,t) with ir and it the radial and temporal indices. Returns P in GBar."""
        return 0
    def vol(self, ir, it):
        """Volume of zone index ir at time index it."""
        return 0
    def rhoR(self, ir, it):
        """Areal density of zone index ir at time index it. Returns rhoR in g/cm2"""
        return 0

    # required time limits in the problem
    def it_min(self):
        """Minimum time for post-proc calculations."""
        return 0
    def it_tc(self):
        """Time index corresponding to right before shock coalescence."""
        return 0
    def it_max(self):
        """Maximum time for post-proc calculations."""
        return 0
    
    # required length limits in the problem
    def ir_min(self):
        """Minimum radius for post-proc calculations."""
        return 0
    def ir_fuel(self):
        """Maximum radius of fuel material."""
        return 0
    def ir_max(self):
        """Maximum radius for post-proc calculations."""
        return 0

    # Convert spatial and temporal indices to real values:
    def t(self, it):
        """Real time in s corresponding to index value it."""
        return 0
    def dt(self):
        """Get the post-processor time step in s."""
        return 0
    def r(self ,ir, it):
        """Real radius in cm corresponding to index value ir at time given by index it."""
        return 0

        
    # required material composition info
    def IonA(self, ir, it):
        """List of AMU masses for fuel ions. At zone indices ir and it."""
        return []
    def IonZ(self, ir, it):
        """List of ion Z for fuel ions. At zone indices ir and it."""
        return []
    def IonF(self, ir, it):
        """List of fuel ion relative populations. At zone indices ir and it."""
        return []
    def Abar(self, ir, it):
        """Average ion A, at zone indices ir and it."""
        return 0
    def Zbar(self, ir, it):
        """Average ion Z, at zone indices ir and it."""
        return 0
    def f(self, ir, it, A, Z):
        """Get composition fraction for fuel ion with A and Z at radius/time indices ir and it."""
        return 0