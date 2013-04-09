# Python-based abstract representation of an implosion
# All implosion types must implement these methods
# A. Zylstra 2012/08/07

from abc import ABCMeta

class Implosion(metaclass=ABCMeta):
    """Abstract class for implementing implosions."""
    name = ""
    
    # initialization
    def __init__(self, inName):
        self.name = inName

    # required hydro variables
    def ni(self, r, t):
        """Ion number density ni(r,t) with r in cm and t in s Returns ni in 1/cm^3."""
        return 0
    def ne(self, r, t):
        """Electron number density ne(r,t) with r in cm and t in s Returns ne in 1/cm^3."""
        return 0
    def Ti(self, r, t):
        """Ion temperature with r in cm and t in s Returns Ti in keV."""
        return 0
    def Te(self, r, t):
        """Electron temperature with r in cm and t in s Returns Te in keV."""
        return 0
    def u(self, r, t):
        """Fluid velocity u(r,t) with r in cm and t in s. Returns u in um/ns."""
        return 0
    def c(self, r, t):
        """Sound speed c(r,t) with r in cm and t in s Returns c in um/ns."""
        return 0
    def rho(self, r, t):
        """Density rho(r,t) with r in cm and t in s Returns rho in g/cm3."""
        return 0
    def P(self, r, t):
        """Pressure P(r,t) with r in cm and t in s Returns P in GBar."""
        return 0

    # required time limits in the problem
    def tmin(self):
        """Minimum time for post-proc calculations."""
        return 0
    def tmax(self):
        """Maximum time for post-proc calculations."""
        return 0
    
    # required length limits in the problem
    def rmin(self, t):
        """Minimum radius for post-proc calculations at time t in s."""
        return 0
    def rmax(self, t):
        """Maximum radius for post-proc calculations at time t in s."""
        return 0
        
    # required material composition info
    def IonA(self, r, t):
        """List of AMU masses for fuel ions. At r in cm and t in s."""
        return []
    def IonZ(self, r, t):
        """List of ion Z for fuel ions. At r in cm and t in s."""
        return []
    def IonF(self, r, t):
        """List of fuel ion relative populations. At r in cm and t in s."""
        return []