
from impy.implosions.Implosion import Implosion
from impy.gui.GenerateGuderley import Generate_Guderley
from impy.resources.constants import Erg2eV, AMU

import math
import scipy
import scipy.optimize
import scipy.integrate
import scipy.interpolate
import numpy as np

class Guderley(Implosion):
    """Python-based abstract representation of an implosion. All implosion types must implement these methods.
    Each implosion must support a few options for constructing.

    :param type: The type of constructor to be used. Available options are:

    'GUI': The class should interact directly with the user to get required options (default)

    'File': Create directly from a file. If the function `getFileTypes` returns a 0-length list, this will not be called

    'CLI': Interact with user via CLI, or take info from args.

    :param args: Additional information, which depends on the type of constructor:

    type='GUI': unused

    type='File': A full path to the file to open.

    type='CLI': A full list of arguments passed to the executable, to be interpreted as this Implosion pleases.

    :param wm: (optional) Window manager to use for displaying windows

    The behavior of all functions taking indices (time and radius, `it` and `ir` respectively) is as follows.

    Both indices can be single integers, in which case the return type will be a scalar number unless otherwise noted::

        >>> # foo is an Implosion
        >>> bar = foo.ni(5,5)
        >>> print(bar)
        5.

    Both indices can be length-2 `tuples` containing a range [min,max) of desired indices to sample at.
    In this case a 2-D :py:class:`numpy.ndarray` is returned where the first index corresponds to the time indices and
    the second axis corresponds to the radial indices::

        >>> # foo is an Implosion
        >>> it = (1,5)
        >>> ir = (8,10)
        >>> data = foo.ni(it, ir)
        >>> type(data)
        <class 'numpy.ndarray'>
        >>> data.shape
        (4,2)
        >>> data[0,0] == foo.ni(1,8)
        True

    One note is that there are several functions for material composition (`IonA`, `IonF`, `IonZ`) which return
    :py:class:`numpy.ndarray` for a single zone. These can also be called with tuples, in which case they return
    a 3-D array where the third axis has length equal to the maximum number of ions in a zone. Some zones may have
    padding zeros in this array. The other material composition functions behave like the other functions above,
    i.e. they can be called with either scalar or tuple arguments, and return either
    a scalar or a 2-D :py:class:`numpy.ndarray` (respectively).

    :author: Alex Zylstra
    :date: 2014-05-15
    """
    __author__ = 'Alex Zylstra'
    __date__ = '2014-05-15'
    __version__ = '0.2.0'

    # ----------------------------------------
    #           Generic methods
    # ----------------------------------------
    def __init__(self, type='GUI', args='', wm=None):
        """Construct a new implosion."""
        super(Guderley, self).__init__()
        self.__ready__ = False

        # mass densities in g/cc for one atm of gas
        pD = 0.08988*2/1000
        pT = 0.08988*3/1000
        p3He = 0.1786*(3/4)/1000

        # Some Guderley constants:
        self.gamma = 5./3. #ideal monoatomic gas
        self.n = 3. #spherical geometry
        self.alpha = 0.68838

        # Get file to open based on type:
        if type is 'GUI':
            dialog = Generate_Guderley(None)
            if dialog.cancelled or dialog.result is None or len(dialog.result) < 12:
                return

            r0 = dialog.result[0]
            tc = dialog.result[1]
            xish = dialog.result[2]
            ion1 = dialog.result[3]
            p1 = dialog.result[4]
            ion2 = dialog.result[5]
            p2 = dialog.result[6]
            eiCoup = dialog.result[7]
            dr = dialog.result[8]
            t0 = dialog.result[9]
            t1 = dialog.result[10]
            dt = dialog.result[11]

        elif type is 'CLI':
            r0 = float(input("Radius (um): "))
            tc = float(input("Collapse time (ns): "))
            xish = float(input("Shock strength (um/ns^a): "))
            ions = ['D','3He','T']
            ion1 = ''
            while ion1 not in ions:
                ion1 = input("Ion 1? (D/3He/T) ")
            ion2 = ''
            while ion2 not in ions:
                ion2 = input("Ion 2? (D/3He/T) ")

            p1 = float(input(ion1+" fill (atm): "))
            p2 = float(input(ion2+" fill (atm): "))
            eiCoup = float(input("Rygg coupling parameter: "))
            dr = float(input("Radial step (um): "))
            t0 = float(input("t_min (ns): "))
            t1 = float(input("t_max (ns): "))
            dt = float(input("dt (ns): "))

        else:
            raise ValueError('Type passed to Guderley constructor is invalid: ' + type)

        # Set instance variables:
        self.r0 = pow(10,-4.)*r0
        self.tc = pow(10,-9)*tc
        self.xish = pow(10,-4)*pow(10,9*self.alpha)*xish
        self.xirsh = 0.740 * self.xish
        self.t0 = pow(self.r0/self.xish,1./self.alpha)
        self.t_min = t0 * 1e-9
        self.t_max = t1 * 1e-9
        self.dr = dr * 1e-4
        self.tstep = dt * 1e-9

        # this is helpful numerically. Basically, it keeps tc from being too close to time steps, since
        # (t-tc) occurs in some denominators
        if math.fabs((self.tc-self.t_min) % self.tstep - self.tstep) <1e-12:
            self.tc += 1e-12

        #fuel info
        self.f1 = p1
        if ion1 == 'D':
            self.rho0 = self.f1*pD
            self.Ion1A = 2
            self.Ion1Z = 1
        elif ion1 == 'T':
            self.rho0 = self.f1*pT
            self.Ion1A = 3
            self.Ion1Z = 1
        elif ion1 == '3He':
            self.rho0 = self.f1*p3He
            self.Ion1A = 3
            self.Ion1Z = 2
        self.f2 = p2
        if ion2 == 'D':
            self.rho0 = self.f2*pD
            self.Ion2A = 2
            self.Ion2Z = 1
        elif ion2 == 'T':
            self.rho0 = self.f2*pT
            self.Ion2A = 3
            self.Ion2Z = 1
        elif ion2 == '3He':
            self.rho0 = self.f2*p3He
            self.Ion2A = 3
            self.Ion2Z = 2
        # figure out total 'pressure' including molecular dissociation
        temp = self.f1 + self.f2
        if ion1 == 'D' or ion1 == 'T':
            temp += self.f1
        if ion2 == 'D' or ion2 == 'T':
            temp += self.f2
        # calculate normalized fuel fractions:
        if ion1 == 'D' or ion1 == 'T':
            self.f1 = self.f1*2 / temp
        else:
            self.f1 = self.f1 / temp
        if ion2 == 'D' or ion2 == 'T':
            self.f2 = self.f2*2 / temp
        else:
            self.f2 = self.f2 / temp
        # Calculate average A and Z:
        self.FuelA = 0
        self.FuelZ = 0
        if ion1 == 'D':
            self.FuelA += 2*self.f1
            self.FuelZ += 1*self.f1
        elif ion1 == 'T':
            self.FuelA += 3*self.f1
            self.FuelZ += 1*self.f1
        elif ion1 == '3He':
            self.FuelA += 3*self.f1
            self.FuelZ += 2*self.f1
        if ion2 == 'D':
            self.FuelA += 2*self.f2
            self.FuelZ += 1*self.f2
        elif ion2 == 'T':
            self.FuelA += 3*self.f2
            self.FuelZ += 1*self.f2
        elif ion2 == '3He':
            self.FuelA += 3*self.f2
            self.FuelZ += 2*self.f2

        self.eiCoup = eiCoup

        # Some Guderley constants
        self.gKappa = 0.
        self.gLambda = 1./self.alpha - 1.
        self.gEpsilon = self.gKappa*(self.gamma-1)+2*self.gLambda
        self.gMu = 2/(self.gamma-1)
        self.gBeta = self.n - self.gMu*self.gLambda
        self.gNu = self.n*self.gamma + self.gKappa - 2*self.gLambda

        #Singular points in UC plane
        self.P1=[0,1]
        self.P2=[0,0]
        self.P3=[0,0]
        self.P4=[0,0]
        self.P5=[0,0]
        self.P6=[0,0]
        #Other points
        self.A = [math.sqrt(2*self.gamma*(self.gamma-1))/(self.gamma+1) , 2/(self.gamma+1)]
        self.S1=[10,0]
        self.S2=[0,0]

        #Trajectories in UC plane
        self.UincList = []
        self.UincInt = 0
        self.UouterInt = 0
        self.UcentralInt = 0

        #Trajectory
        self.rShellInt = 0 #for shell trajectory
        self.rLagrangeInt = 0 #for Lagrange plots

        #List of points for C vs xi
        self.CvsXiShock = []
        self.CvsXiShockInt = 0
        self.CvsXiCentral = []
        self.CvsXiCentralInt = 0
        self.CvsXiOuter = []
        self.CvsXiOuterInt = 0
        self.GA = (self.gamma+1)/(self.gamma-1)
        self.K3 = self.GA*pow(self.alpha*self.A[0],-1.*self.gMu*(self.n+self.gKappa)/self.gBeta)*pow(1-self.A[1],-1.*((self.gKappa+self.gMu*self.gLambda)/self.gBeta))

        self.__ready__ = True

    @classmethod
    def getFileTypes(cls):
        """Get a list containing extensions of file types supported by this implosion.
        Must be an array of dicts ready to pass to file dialog, e.g.::

            [('description', '*.ext')]
        """
        return []

    @classmethod
    def name(cls):
        """Get a string containing a name for this type of implosion."""
        return 'Guderley'

    def info(self):
        """Get a string of information about this specific implosion."""
        #TODO: implement more interesting info
        return 'Guderley'

    def ready(self):
        """Returns true if implosion object creation went OK and this object is ready for `generate` to be called."""
        return self.__ready__

    def generate(self):
        """Run the calculation to generate the implosion data."""
        # Top level construction, uses several helper functions:
        self.runProgress = 0
        self.__run__()

    def progress(self):
        """Get the implosion generation's progress estimate.

        :returns: Scalar number between 0 and 1.
        """
        return self.runProgress

    # ----------------------------------------
    #       Hydrodynamic Quantities
    # ----------------------------------------
    def ni(self, it, ir):
        """Ion number density ni(t,r).

         :param it: the temporal index
         :param ir: the spatial index
         :returns: ni [1/cm^3]
         """
        it, ir = self.__internalIndex__(it, ir)

        if np.isscalar(it) and np.isscalar(ir):
            return self.ni_raw[it,ir]
        return self.ni_raw[it[0]:it[1], ir[0]:ir[1]]

    def ne(self, it, ir):
        """Electron number density ne(t,r).

         :param it: the temporal index
         :param ir: the spatial index
         :returns: ne [1/cm^3]
         """
        it, ir = self.__internalIndex__(it, ir)

        if np.isscalar(it) and np.isscalar(ir):
            return self.ne_raw[it,ir]
        return self.ne_raw[it[0]:it[1], ir[0]:ir[1]]

    def Ti(self, it, ir):
        """Ion temperature.

         :param it: the temporal index
         :param ir: the spatial index
         :returns: Ti [keV]
         """
        it, ir = self.__internalIndex__(it, ir)

        if np.isscalar(it) and np.isscalar(ir):
            return self.Ti_raw[it,ir]
        return self.Ti_raw[it[0]:it[1], ir[0]:ir[1]]

    def Te(self, it, ir):
        """Electron temperature.

         :param it: the temporal index
         :param ir: the spatial index
         :returns: Te [keV]
         """
        it, ir = self.__internalIndex__(it, ir)

        if np.isscalar(it) and np.isscalar(ir):
            return self.Te_raw[it,ir]
        return self.Te_raw[it[0]:it[1], ir[0]:ir[1]]

    def u(self, it, ir):
        """Fluid velocity u(t,r).

         :param it: the temporal index
         :param ir: the spatial index
         :returns: u [um/ns]
         """
        it, ir = self.__internalIndex__(it, ir)

        if np.isscalar(it) and np.isscalar(ir):
            return self.u_raw[it,ir]
        return self.u_raw[it[0]:it[1], ir[0]:ir[1]]

    def c(self, it, ir):
        """Sound speed c(t,r).

         :param it: the temporal index
         :param ir: the spatial index
         :returns: c [um/ns]
         """
        # Using helper functions, statement below works for either scalar or tuple index calls:
        # Have to convert P from GBar to CGS (bayre)
        return np.sqrt( (5/3) * self.P(it,ir) * 1e15 / self.rho(it,ir) )

    def rho(self, it, ir):
        """Mass density rho(t,r).

         :param it: the temporal index
         :param ir: the spatial index
         :returns: rho [g/cm3]
         """
        it, ir = self.__internalIndex__(it, ir)

        if np.isscalar(it) and np.isscalar(ir):
            return self.rho_raw[it,ir]
        return self.rho_raw[it[0]:it[1], ir[0]:ir[1]]

    def P(self, it, ir):
        """Pressure P(t,r).

         :param it: the temporal index
         :param ir: the spatial index
         :returns: P [GBar]
         """
        it, ir = self.__internalIndex__(it, ir)

        if np.isscalar(it) and np.isscalar(ir):
            return self.P_raw[it,ir]
        return self.P_raw[it[0]:it[1], ir[0]:ir[1]]

    def vol(self, it, ir):
        """Volume of zone index ir at time index it.

         :param it: the temporal index
         :param ir: the spatial index
         :returns: Zone volume [cm^3]
         """
        it, ir = self.__internalIndex__(it, ir)

        if np.isscalar(it) and np.isscalar(ir):
            return self.vol_raw[it,ir]
        return self.vol_raw[it[0]:it[1], ir[0]:ir[1]]

    def rhoR(self, it, ir):
        """Areal density.

         :param it: the temporal index
         :param ir: the spatial index
         :returns: rhoR [g/cm2]
         """
        it, ir = self.__internalIndex__(it, ir)

        if np.isscalar(it) and np.isscalar(ir):
            return self.rhoR_raw[it,ir]
        return self.rhoR_raw[it[0]:it[1], ir[0]:ir[1]]

    # ----------------------------------------
    #       Time and spatial 'scales'
    # ----------------------------------------
    def it_min(self):
        """Minimum time (inclusive) in this implosion.

        :returns: The lowest acceptable value in `it` (typically will be 0)
        """
        return 0

    def it_tc(self):
        """Time right before shock coalescence.

        :returns: The time index for tc
        """
        return int((self.tc - self.t_min)/self.tstep)

    def it_tstag(self):
        """Time corresponding to stagnation, defined as the minimum shell radius.

        :returns: The time index for tstag
        """
        return np.argmin(self.r_raw[:,-1])

    def it_max(self):
        """Maximum time index (exclusive) in this implosion. Chosen to be exclusive so one can::

            for it in range(foo.it_min(), foo.it_max()):
                bar(it)

        :returns: The time index
        """
        return len(self.t_raw)

    def ir_min(self):
        """Minimum radius (inclusive) in this implosion.

        :returns: The lowest acceptable value in `ir` (typically will be 0)
        """
        return 0

    def ir_fuel(self):
        """Outer radius of fuel material.

        :returns: A radial index
        """
        return len(self.r_raw[0])

    def ir_max(self):
        """Maximum radial index (exclusive) in this implosion. Chosen to be exclusive so one can::

            for ir in range(foo.ir_min(), foo.ir_max()):
                bar(ir)

        :returns: The radial index
        """
        return len(self.r_raw[0])

    def t(self, it):
        """Convert indices to real time.

        :param it: The temporal index
        :returns: real time [s]
        """
        if np.isscalar(it):
            return self.t_min + self.tstep*it
        return np.arange(it[0], it[1])*self.tstep + self.t_min


    def dt(self, it, ir=None):
        """Get the post-processor time step.

        :param it: Time index (may be ignored if implosions have a constant time step)
        :param ir: Optional since dt is the same for all spatial zones. However, this gives the option of getting a
            2-D array as the return value, which is convenient for some calculations.
        :returns: The delta in time between it and it+1 [s]
        """
        if ir is None:
            if np.isscalar(it):
                return self.tstep
            return self.tstep * np.ones((it[1]-it[0]), dtype=np.float)

        # two axis:
        it, ir = self.__internalIndex__(it, ir)

        if np.isscalar(it) and np.isscalar(ir):
            return self.tstep
        return self.tstep * np.ones((it[1]-it[0], ir[1]-ir[0]), dtype=np.float)

    def r(self, it, ir):
        """Get physical radius for a zone.

        :param it: The temporal index
        :param ir: The spatial index
        :returns: r [cm]
        """
        it, ir = self.__internalIndex__(it, ir)

        if np.isscalar(it) and np.isscalar(ir):
            return self.r_raw[it,ir]
        return self.r_raw[it[0]:it[1], ir[0]:ir[1]]

    # ----------------------------------------
    #           Material info
    # ----------------------------------------
    def IonA(self, it, ir):
        """Masses for all ions in a zone.

        :param it: The temporal index
        :param ir: The spatial index
        :returns: A :py:class:`numpy.ndarray` of ion masses [AMU]
        """
        it, ir = self.__internalIndex__(it, ir)

        if np.isscalar(it) and np.isscalar(ir):
            return [self.Ion1A, self.Ion2A]
        ret = np.ones((it[1]-it[0],ir[1]-ir[0],2), dtype=np.float)
        ret[:,:,0] = self.Ion1A
        ret[:,:,1] = self.Ion2A
        return ret

    def IonZ(self, it, ir):
        """Ion atomic numbers for ions in a zone.

        :param it: The temporal index
        :param ir: The spatial index
        :returns: A :py:class:`numpy.ndarray` of ion atomic numbers [e]
        """
        it, ir = self.__internalIndex__(it, ir)

        if np.isscalar(it) and np.isscalar(ir):
            return [self.Ion1Z, self.Ion2Z]
        ret = np.ones((it[1]-it[0],ir[1]-ir[0],2), dtype=np.float)
        ret[:,:,0] = self.Ion1Z
        ret[:,:,1] = self.Ion2Z
        return ret

    def IonF(self, it, ir):
        """Ion fractions in a zone. Each fraction is between 0 and 1 (inclusive). The sum of all fractions is 1.

        :param it: The temporal index
        :param ir: The spatial index
        :returns: A :py:class:`numpy.ndarray` of ion fractions
        """
        it, ir = self.__internalIndex__(it, ir)

        if np.isscalar(it) and np.isscalar(ir):
            return [self.f1, self.f2]
        ret = np.ones((it[1]-it[0],ir[1]-ir[0],2), dtype=np.float)
        ret[:,:,0] = self.f1
        ret[:,:,1] = self.f2
        return ret

    def Abar(self, it, ir):
        """Average ion mass in a zone.

        :param it: The temporal index
        :param ir: The spatial index
        :returns: Average ion mass [AMU]
        """
        it, ir = self.__internalIndex__(it, ir)

        if np.isscalar(it) and np.isscalar(ir):
            return self.FuelA
        return np.ones((it[1]-it[0],ir[1]-ir[0]), dtype=np.float) * self.FuelA

    def Zbar(self, it, ir):
        """Average ion Z.

        :param it: The temporal index
        :param ir: The spatial index
        :returns: Average ion atomic number [e]
        """
        it, ir = self.__internalIndex__(it, ir)

        if np.isscalar(it) and np.isscalar(ir):
            return self.FuelZ
        return np.ones((it[1]-it[0], ir[1]-ir[0]), dtype=np.float) * self.FuelZ

    def f(self, it, ir, A, Z):
        """Get fraction of a specified ion in a zone.

        :param it: The temporal index
        :param ir: The spatial index
        :param A: The ion mass you're interested in (scalar)
        :param Z: The ion atomic number you're interested in (scalar)
        :returns: Fraction of that ion in a zone
        """
        it, ir = self.__internalIndex__(it, ir)

        # scalars:
        if np.isscalar(it) and np.isscalar(ir):
            if A == self.Ion1A and Z == self.Ion1Z:
                return self.f1
            if A == self.Ion2A and Z == self.Ion2Z:
                return self.f2
            return 0.

        # matrix
        dim = (it[1]-it[0], ir[1]-ir[0])
        if A == self.Ion1A and Z == self.Ion1Z:
            return np.ones(dim, dtype=np.float) * self.f1
        if A == self.Ion2A and Z == self.Ion2Z:
            return np.ones(dim, dtype=np.float) * self.f2
        return np.zeros(dim, dtype=np.float)

    # ----------------------------------------
    #           Helper functions
    # ----------------------------------------
    def __run__(self):
        """Top-level method to do the Guderley calculation."""
        self.calcSingPoints()
        self.runProgress += 0.2
        self.UCTrajectories()
        self.runProgress += 0.2
        self.CalcUCG()
        self.runProgress += 0.2
        self.rShellInit()
        self.runProgress += 0.2
        self.tFFcalc = self.tFF()

        # Finally, pre-compute relevant hydro variables:
        self.__precompute_hydro__()

    # ------------------------------------
    # Calculate various trajectories and velocities
    # ------------------------------------
    def rs(self,t):
        """Calculate the shock position at time t (s). Returns rs in cm."""
        strength = self.xish
        if t > self.tc:
            strength = self.xirsh
        return strength * pow( math.fabs(t-self.tc) , self.alpha )
    def us(self,t):
        """Calculate the shock velocity at time t (s). Returns us in cm/s."""
        strength = self.xish
        if t > self.tc:
            strength = self.xirsh
        return strength * self.alpha * pow( math.fabs(t-self.tc) , self.alpha-1 )
    def rFF(self,t):
        """Calculate the 'free fall' radius at time t (s). Returns rFF in cm."""
        uFF = self.uShock(self.r0 , (self.tc - self.t0 + 1e-12) )
        ret = max( 0. , (self.r0 - (t-(self.tc-self.t0))*uFF) )
        return ret
    def tFF(self):
        """Calculate when the 'free fall' mass hits the origin. Returns tFF in s."""
        uFF = self.uShock(self.r0 , (self.tc - self.t0 + 1e-12) )
        return ( self.r0/uFF + (self.tc-self.t0) )
    def rShellInit(self):
        """Do initial numerical integration for shell position vs time calculations."""
        times = list(np.arange(self.tc-self.t0,self.tFF(),0.5e-12))
        rList = scipy.integrate.odeint(self.uLab, self.r0, times)
        #fix list formatting
        rList2 = []
        for i in rList:
            rList2.append(i[0])
        #interpolate
        self.rShellInt = scipy.interpolate.InterpolatedUnivariateSpline(times, rList2)
    def rShell(self,t):
        """Shell position at time t (s). Returns rShell in cm."""
        return self.rShellInt(t)[()]
    def rLagrangeInit(self,ri):
        """Integrate a Lagrangian trajectory starting at r0. Must be run betfore calling rLagrange(t). takes r in cm."""
        times = list(np.arange(self.tc-self.t0,self.tFF(),0.5e-12))
        rList = scipy.integrate.odeint(self.uLab, ri, times)
        #fix list formatting
        rList2 = []
        for i in rList:
            rList2.append(i[0])
        #interpolate
        #self.rLagrangeInt = scipy.interpolate.interp1d(times, rList2, kind='cubic')
        self.rLagrangeInt = scipy.interpolate.InterpolatedUnivariateSpline(times, rList2)
        return
    def rLagrange(self,t):
        """Lagrangian trajectory position at time t (s). Returns r in cm."""
        return self.rLagrangeInt(t)[()]

    # ------------------------------------
    # Self-similarity coordinates
    # ------------------------------------
    def xi(self, r, t):
        """Calculate the self-similarity coordinate as a function of r (cm) and t (s)."""
        if t == self.tc:
            return 0
        return (r/self.r0) * np.power( np.fabs( (t-self.tc)/self.t0 ) , -1.*self.alpha )

    def U23(self,x):
        return (self.n-1)*self.gamma*np.power(x,2)+(self.gKappa-2*self.gLambda-self.gamma*(self.n-1-self.gLambda))*x-(self.gKappa-2*self.gLambda)

    def calcSingPoints(self):
        """Calculate the singular points in the UC plane."""
        U2 = scipy.optimize.fsolve(self.U23,0.4)[0]
        U3 = scipy.optimize.fsolve(self.U23,0.6)[0]
        self.P2 = [1-U2,U2]
        self.P3 = [1-U3,U3]
        self.P5 = [(math.sqrt(self.n)/(self.alpha*(self.n+self.gMu)))/math.sqrt(1.+(self.gKappa+self.gMu*self.gLambda)*(self.n+self.gMu)/((self.gMu+2)*(self.n-self.gMu*self.gLambda))) , self.gMu/(self.alpha*(self.n+self.gMu))]
        self.P6 = [10. , -(self.gKappa-2*self.gLambda)/(self.n*self.gamma)]


    # ------------------------------------
    # Some misc helper functions
    # ------------------------------------
    def Delta(self, x1, x2):
        return np.power(x2,2.) - np.power(1.-x1,2.)
    def Delta1(self, x1, x2):
        return x1*(1.-x1)*(1./self.alpha-x1)-pow(x2,2.)*(self.n*x1+(self.gKappa-2*self.gLambda)/self.gamma)
    def Delta2(self, x1, x2):
        return x2*((1.-x1)*(1/self.alpha-x1)+(x1/self.gMu)*(self.gLambda+(self.n-1)*(x1-1))-pow(x2,2.)+(self.gEpsilon/(2.*self.gamma))*pow(x2,2.)/(x1-1.))
    def dUdC(self, x1, x2):
        return self.Delta1(x1,x2) / self.Delta2(x1,x2)
    def ShockJump(self, x):
        """Shock jump in the UC plane. x is a 2D point."""
        ret1 = 2*self.gamma*(self.gamma-1)*pow((1-x[1]),2)/pow(self.gamma+1,2)+pow(x[0],2)*(1-2*pow(((self.gamma-1)/(self.gamma+1)),2)-2*(self.gamma-1)*pow((x[0]/(1-x[1])),2)/pow(self.gamma+1,2))
        ret2 = 1-(1-x[1])*((self.gamma-1)/(self.gamma+1)+2*pow((x[0]/(1-x[1])),2)/(self.gamma+1))
        return [ret1 , ret2]
    def ShockJumpErr1(self, x):
        return math.fabs(self.ShockJump([float(x),float(self.UCcentral(x))])[1] - self.UCouter(x))
    def ShockJumpErr2(self, x):
        return math.fabs(self.ShockJump([float(x),float(self.UCouter(x))])[1] - self.UCcentral(x))

    # ------------------------------------
    # Calculate trajectories in the UC plane
    # ------------------------------------
    def UCTrajectories(self):
        """Solve the Guderley solution for UC trajectories"""
        dInt = 0.0005

        #Incoming shock
        #Find the UC trajectory for shock between A and P3
        intRange = list(np.arange( self.A[0], self.P3[0], -1.0*dInt))
        temp = scipy.integrate.odeint( self.dUdC , self.A[1] , intRange )
        Uinc1 = []
        for i in range(0,len(temp)):
            Uinc1.append([ intRange[i] , temp[i,0] ])
        Uinc1.reverse() #ascending order

        #Find the UC trajectory for the shock between P4 and P3
        intRange = list(np.arange(self.P4[0]+dInt, self.P3[0]-0.02, dInt))
        temp = scipy.integrate.odeint( self.dUdC , self.P4[1]+dInt , intRange )
        Uinc2 = []
        for i in range(0,len(temp)):
            Uinc2.append([ intRange[i] , temp[i,0] ])

        #Stitch the whole trajectory together
        self.UincList = [self.P4] + Uinc2 + [self.P3] + Uinc1
        # set up interpolation
        UincListx = []
        UincListy = []
        for i in self.UincList:
            UincListx.append(i[0])
            UincListy.append(i[1])
        #self.UincInt = scipy.interpolate.interp1d(UincListx, UincListy, kind='linear')
        self.UincInt = scipy.interpolate.InterpolatedUnivariateSpline(UincListx, UincListy)

        #Outer flow region
        #numerical integration
        intRange = list(np.arange( dInt , 1.1 , dInt ))
        temp = scipy.integrate.odeint( self.dUdC , self.P4[0]-dInt , intRange )
        Uouterflow1 = []
        for i in range(0,len(temp)):
            Uouterflow1.append([ intRange[i] , temp[i,0] ])
        # set up interpolation
        x = []
        y = []
        for i in Uouterflow1:
            x.append(i[0])
            y.append(i[1])
        #self.UouterInt = scipy.interpolate.interp1d(x, y, kind='linear')
        self.UouterInt = scipy.interpolate.InterpolatedUnivariateSpline(x, y)

        #Central flow region
        intRange = list(np.arange( self.P6[0] , 1-dInt , -1.*dInt ))
        temp = scipy.integrate.odeint( self.dUdC , self.P6[1] , intRange )
        Ucentralflow1 = []
        for i in range(0,len(temp)):
            Ucentralflow1.append([ intRange[i] , temp[i,0] ])
        # set up interpolation
        x = []
        y = []
        for i in Ucentralflow1:
            x.append(i[0])
            y.append(i[1])
        x.reverse()
        y.reverse()
        #self.UcentralInt = scipy.interpolate.interp1d(x, y, kind='linear')
        self.UcentralInt = scipy.interpolate.InterpolatedUnivariateSpline(x, y)


        #Find the shock jump points
        self.S1 = [10,0]
        self.S2 = [0,0]
        S2C = scipy.optimize.fmin(self.ShockJumpErr1 , 1.089, disp=0)[0]
        S1C = self.ShockJump([S2C,self.UCcentral(S2C)])[0]
        self.S1 = [S1C, self.UCouter(S1C)]
        self.S2 = [S2C, self.UCcentral(S2C)]


    # ------------------------------------
    # Interpolating functions for UC trajectories
    # ------------------------------------
    def UCinc(self, x):
        """U(C) for the incoming shock."""
        if type(x) is list:
            x = x[0]
        return self.UincInt(x)[()]
    def UCouter(self, x):
        """U(C) for the outer flow."""
        if type(x) is list:
            x = x[0]
        if x < self.S1[0] and x > 0 and x < 1.1:
            return self.UouterInt(x)[()]
        return 0
    def UCcentral(self, x):
        """U(C) for the central flow."""
        if type(x) is list:
            x = x[0]
        if x > self.S2[0] and x < self.P6[0] and x > 1:
            return self.UcentralInt(x)[()]
        return 0

    # ------------------------------------
    # Initial calculation of U,C,G vs self-similarity index
    # ------------------------------------
    def CalcUCG(self):
        """Calculate U,C,G vs xi."""
        #Incoming shock
        dXi=0.005
        XiMin=1
        XiMax=100
        self.CvsXiShock.append([XiMin,self.A[0]])
        for i in range(1, int(1 + (XiMax-XiMin)/dXi)):
            Xii = self.CvsXiShock[i-1][0]
            Ci = self.CvsXiShock[i-1][1]
            self.CvsXiShock.append([Xii+dXi,Ci+self.Delta2(self.UCinc(Ci),Ci)/self.Delta(self.UCinc(Ci),Ci)*(math.log(Xii+dXi) - math.log(Xii))])
        #interpolation
        x = []
        y = []
        for i in self.CvsXiShock:
            x.append(i[0])
            y.append(i[1])
        #self.CvsXiShockInt = scipy.interpolate.interp1d(x, y, kind='linear')
        self.CvsXiShockInt = scipy.interpolate.InterpolatedUnivariateSpline(x, y)

        #Central flow region
        dXi = 0.001
        XiMin = dXi
        XiMax = 0.74
        self.CvsXiCentral.append([XiMax, self.S2[0]])
        for i in range(1, int(1+(XiMax-XiMin)/dXi)):
            Xii = self.CvsXiCentral[i-1][0]
            Ci = self.CvsXiCentral[i-1][1]
            self.CvsXiCentral.append([Xii-dXi,Ci-self.Delta2(self.UCcentral(Ci),Ci)/self.Delta(self.UCcentral(Ci),Ci)*(math.log(Xii+dXi) - math.log(Xii))])
        self.CvsXiCentral.reverse()
        #interpolation
        x = []
        y = []
        for i in self.CvsXiCentral:
            x.append(i[0])
            y.append(i[1])
        #self.CvsXiCentralInt = scipy.interpolate.interp1d(x, y, kind='linear')
        self.CvsXiCentralInt = scipy.interpolate.InterpolatedUnivariateSpline(x, y)

        #Outer flow region
        dXi = 0.005
        XiMin = 0.74
        XiMax = 100
        self.CvsXiOuter.append([XiMin,self.S1[0]])
        for i in range(1, int(1+(XiMax-XiMin)/dXi)):
            Xii = self.CvsXiOuter[i-1][0]
            Ci = self.CvsXiOuter[i-1][1]
            self.CvsXiOuter.append([Xii+dXi,Ci+self.Delta2(self.UCouter(Ci),Ci)/self.Delta(self.UCouter(Ci),Ci)*(math.log(Xii+dXi) - math.log(Xii))])
        #interpolation
        x = []
        y = []
        for i in self.CvsXiOuter:
            x.append(i[0])
            y.append(i[1])
        #self.CvsXiOuterInt = scipy.interpolate.interp1d(x, y, kind='linear')
        self.CvsXiOuterInt = scipy.interpolate.InterpolatedUnivariateSpline(x, y)

        return

    # ------------------------------------
    # U,C,G vs self-similarity index
    # ------------------------------------
    #incoming shock
    def Cshock(self, x):
        """C(xi) for the incoming shock."""
        if x > 1 and x < 100:
            return self.CvsXiShockInt(x)[()]
        return 0
    def Ushock(self, x):
        """U(xi) for the incoming shock."""
        if x > 1 and x < 100:
            return self.UCinc(self.CvsXiShockInt(x)[()])
        return 0
    def Gshock(self, x):
        """G(xi) for the incoming shock."""
        if x > 1 and x < 100:
            return self.K3*pow(self.alpha*pow(x,1/self.alpha)*self.Cshock(x),self.gMu*(self.n+self.gKappa)/self.gBeta)*pow(1-self.Ushock(x),(self.gKappa+self.gMu*self.gLambda)/self.gBeta)
        return 0
    #central flow region
    def Ccentral(self, x):
        """C(xi) for the central flow region."""
        if x > 0 and x <= 0.74:
            return self.CvsXiCentralInt(x)[()]
        return 0
    def Ucentral(self, x):
        """U(xi) for the central flow region."""
        if x > 0 and x <= 0.74:
            return self.UCcentral(self.CvsXiCentralInt(x)[()])
        return 0
    def Gcentral(self, x):
        """G(xi) for the central flow region."""
        if x > 0 and x <= 0.74:
            return self.K3*pow(self.alpha*pow(x,1/self.alpha)*self.Ccentral(x),self.gMu*(self.n+self.gKappa)/self.gBeta)*pow(1-self.Ucentral(x),(self.gKappa+self.gMu*self.gLambda)/self.gBeta)
        return 0
    #outer flow region
    def Couter(self, x):
        """C(xi) for the outer flow region."""
        if x > 0.74 and x < 100:
            return self.CvsXiOuterInt(x)[()]
        return 0
    def Uouter(self, x):
        """U(xi) for the outer flow region."""
        if x > 0.74 and x < 100:
            return self.UCouter(self.CvsXiOuterInt(x)[()])
        return 0
    def Gouter(self, x):
        """G(xi) for the outer flow region."""
        if x > 0.74 and x < 100:
            return self.K3*pow(self.alpha*pow(x,1/self.alpha)*self.Couter(x),self.gMu*(self.n+self.gKappa)/self.gBeta)*pow(1-self.Uouter(x),(self.gKappa+self.gMu*self.gLambda)/self.gBeta)
        return 0
    #Overall U,C,G vs r and t
    def C(self, r, t):
        """Calculate C(r,t) with r in cm and t in s."""
        if t < self.tc: #incoming shock
            return self.Cshock( self.xi(r,t) )
        if r < self.rs(t): #central flow
            return self.Ccentral( self.xi(r,t) )
        if r >= self.rs(t): #outer flow
            return self.Couter (self.xi(r,t) )
    def U(self, r, t):
        """Calculate U(r,t) with r in cm and t in s"""
        if t < self.tc: #incoming shock
            return self.Ushock( self.xi(r,t) )
        if r < self.rs(t): #central flow
            return self.Ucentral( self.xi(r,t) )
        if r >= self.rs(t): #outer flow
            return self.Uouter (self.xi(r,t) )
    def G(self, r, t):
        """Calculate G(r,t) with r in cm and t in s"""
        if t < self.tc: #incoming shock
            return self.Gshock( self.xi(r,t) )
        if r < self.rs(t): #central flow
            return self.Gcentral( self.xi(r,t) )
        if r >= self.rs(t): #outer flow
            return self.Gouter (self.xi(r,t) )

    # convenient to have velocity in two frames:
    def uShock(self, r, t):
        """Fluid velocity u(r,t) with r in cm and t in s. Returns u in um/ns, relative to shock direction."""
        if t != self.tc:
            return self.alpha*r/math.fabs(t-self.tc)*self.U(r,t)
        return 0
    def uLab(self, r, t):
        """Fluid velocity u(r,t) with r in cm and t in s. Returns u in um/ns, in lab frame."""
        ret = self.uShock(r,t)
        if t != self.tc:
            if t < self.tc: #incoming shock
                return -1.*ret
            if r < self.rs(t): #central flow
                return ret
            if r >= self.rs(t): #outer flow
                return ret
        return 0

    def __precompute_hydro__(self):
        """Precompute the hydro variables for faster use later."""
        import math
        dim = ( math.ceil((self.t_max-self.t_min)/self.tstep) , int(self.r0/self.dr) )

        self.t_raw = np.arange(self.t_min, self.t_max, self.tstep)
        self.r_raw = np.ndarray(dim, dtype=np.float)
        self.u_raw = np.ndarray(dim, dtype=np.float)
        self.c_raw = np.ndarray(dim, dtype=np.float)
        self.T_raw = np.ndarray(dim, dtype=np.float)
        self.rho_raw = np.ndarray(dim, dtype=np.float)
        self.P_raw = np.ndarray(dim, dtype=np.float)

        # Set initial conditions
        self.r_raw[0,:] = np.arange(self.dr, self.r0+self.dr, self.dr)
        for j in range(dim[1]):
            self.u_raw[0,j] = self.uLab(self.r_raw[0,j], self.t_raw[0])
            self.c_raw[0,j] = self.alpha * self.r_raw[0,j] / math.fabs(self.t_raw[0]-self.tc) * self.C(self.r_raw[0,j], self.t_raw[0])
            self.T_raw[0,j] = pow(10,-3) * self.FuelA * AMU * pow(self.c_raw[0,j],2) / Erg2eV
            self.rho_raw[0,j] = self.rho0 * np.power((self.r_raw[0,j]/self.r0), self.gKappa) * self.G(self.r_raw[0,j], self.t_raw[0])
            self.P_raw[0,j] = (1 + self.FuelZ) * self.rho_raw[0,j] * np.power(self.c_raw[0,j],2) * pow(10,-15)

        # Loop to set the rest
        for i in range(1,dim[0]):
            for j in range(dim[1]):
                self.r_raw[i,:] = np.arange(self.dr, self.r0+self.dr, self.dr)
                #self.r_raw[i,j] = self.r_raw[i-1,j] + self.u_raw[i-1,j]*self.tstep
                self.u_raw[i,j] = self.uLab(self.r_raw[i,j], self.t_raw[i])
                self.c_raw[i,j] = self.alpha * self.r_raw[i,j] / math.fabs(self.t_raw[i]-self.tc) * self.C(self.r_raw[i,j], self.t_raw[i])
                self.T_raw[i,j] = pow(10,-3) * self.FuelA * AMU * np.power(self.c_raw[i,j],2) / Erg2eV
                if self.T_raw[i,j] == 0:
                    self.T_raw[i,j] = 0.025
                self.rho_raw[i,j] = self.rho0 * pow((self.r_raw[i,j]/self.r0), self.gKappa) * self.G(self.r_raw[i,j], self.t_raw[i])
                if self.rho_raw[i,j] == 0:
                    self.rho_raw[i,j] = self.rho0
                self.P_raw[i,j] = (1 + self.FuelZ) * self.rho_raw[i,j] * np.power(self.c_raw[i,j],2) * pow(10,-15)

        # Some can be calculated more elegantly:
        self.Ti_raw = np.ndarray(dim, dtype=np.float)
        self.Te_raw = np.ndarray(dim, dtype=np.float)
        self.ni_raw = np.ndarray(dim, dtype=np.float)
        self.ne_raw = np.ndarray(dim, dtype=np.float)

        self.ni_raw = self.rho_raw / (self.FuelA * AMU)
        self.ne_raw = self.ni_raw * self.FuelZ
        self.Ti_raw = (0.001 * self.FuelA * AMU * np.power(self.c_raw,2) / Erg2eV) * (1.0 - self.eiCoup * self.FuelZ / (1.+self.FuelZ))
        self.Te_raw = (self.T_raw - self.Ti_raw) / self.FuelZ
        # make sure temperature ends up non-zero:
        self.Ti_raw = np.maximum(self.Ti_raw, 0.025*np.ones_like(self.Ti_raw))
        self.Te_raw = np.maximum(self.Te_raw, 0.025*np.ones_like(self.Te_raw))

        # Calculate the volume and rhoR in each zone at each time:
        self.vol_raw = np.ndarray(dim, dtype=np.float)
        self.rhoR_raw = np.ndarray(dim, dtype=np.float)
        for i in range(dim[0]):
            for j in range(dim[1]):
                if j == 0:
                    self.vol_raw[i,j] = (4*np.pi/3) * np.power(self.r_raw[i,j],3)
                    self.rhoR_raw[i,j] = self.rho_raw[i,j] * self.r_raw[i,j]
                else:
                    r0 = self.r_raw[i,j-1]
                    r1 = self.r_raw[i,j]
                    self.vol_raw[i,j] = (4*np.pi/3) * (np.power(r1,3)-np.power(r0,3))
                    self.rhoR_raw[i,j] = self.rho_raw[i,j] * (r1-r0)


