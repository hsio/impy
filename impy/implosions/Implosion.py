""" Python-based abstract representation of an implosion. All implosion types must implement these methods.

:author: Alex Zylstra
:date: 2014-01-23
"""
__author__ = 'Alex Zylstra'
__date__ = '2014-01-23'
__version__ = '1.0.0'

from abc import ABCMeta, abstractmethod
import numpy as np
import inspect
from impy.resources import fusion

class Implosion(metaclass=ABCMeta):
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
    :date: 2014-01-23
    """
    __author__ = 'Alex Zylstra'
    __date__ = '2014-01-23'
    __version__ = '1.0.0'
    
    # ----------------------------------------
    #           Generic methods
    # ----------------------------------------
    @abstractmethod
    def __init__(self, type='GUI', args='', wm=None):
        """Construct a new implosion."""
        self.yieldData = dict()

    @classmethod
    @abstractmethod
    def getFileTypes(cls):
        """Get a list containing extensions of file types supported by this implosion.
        Must be an array of dicts ready to pass to file dialog, e.g.::

            [('description', '*.ext')]
        """
        pass

    @classmethod
    @abstractmethod
    def name(cls):
        """Get a string containing a name for this type of implosion."""
        pass

    @abstractmethod
    def info(self):
        """Get a string of information about this specific implosion."""
        pass

    @abstractmethod
    def generate(self):
        """Run the calculation to generate the implosion data."""
        pass

    @abstractmethod
    def progress(self):
        """Get the implosion generation's progress estimate.

        :returns: Scalar number between 0 and 1.
        """
        pass

    # ----------------------------------------
    #       Hydrodynamic Quantities
    # ----------------------------------------
    @abstractmethod
    def ni(self, it, ir):
        """Ion number density ni(t,r).

         :param it: the temporal index
         :param ir: the spatial index
         :returns: ni [1/cm^3]
         """
        pass
    
    @abstractmethod
    def ne(self, it, ir):
        """Electron number density ne(t,r).

         :param it: the temporal index
         :param ir: the spatial index
         :returns: ne [1/cm^3]
         """
        pass
    
    @abstractmethod
    def Ti(self, it, ir):
        """Ion temperature.

         :param it: the temporal index
         :param ir: the spatial index
         :returns: Ti [keV]
         """
        pass
    
    @abstractmethod
    def Te(self, it, ir):
        """Electron temperature.

         :param it: the temporal index
         :param ir: the spatial index
         :returns: Te [keV]
         """
        pass
    
    @abstractmethod
    def u(self, it, ir):
        """Fluid velocity u(t,r).

         :param it: the temporal index
         :param ir: the spatial index
         :returns: u [um/ns]
         """
        pass
    
    @abstractmethod
    def c(self, it, ir):
        """Sound speed c(t,r).

         :param it: the temporal index
         :param ir: the spatial index
         :returns: c [um/ns]
         """
        pass
    
    @abstractmethod
    def rho(self, it, ir):
        """Mass density rho(t,r).

         :param it: the temporal index
         :param ir: the spatial index
         :returns: rho [g/cm3]
         """
        pass
    
    @abstractmethod
    def P(self, it, ir):
        """Pressure P(t,r).

         :param it: the temporal index
         :param ir: the spatial index
         :returns: P [GBar]
         """
        pass
    
    @abstractmethod
    def vol(self, it, ir):
        """Volume of zone index ir at time index it.

         :param it: the temporal index
         :param ir: the spatial index
         :returns: Zone volume [cm^3]
         """
        pass
    
    @abstractmethod
    def rhoR(self, it, ir):
        """Areal density.

         :param it: the temporal index
         :param ir: the spatial index
         :returns: rhoR [g/cm2]
         """
        pass


    # ----------------------------------------
    #       Time and spatial 'scales'
    # ----------------------------------------
    @abstractmethod
    def it_min(self):
        """Minimum time (inclusive) in this implosion.

        :returns: The lowest acceptable value in `it` (typically will be 0)
        """
        pass

    @abstractmethod
    def it_tc(self):
        """Time right before shock coalescence.

        :returns: The time index for tc
        """
        pass

    @abstractmethod
    def it_tstag(self):
        """Time corresponding to stagnation, defined as the minimum shell radius.

        :returns: The time index for tstag
        """
        pass

    @abstractmethod
    def it_max(self):
        """Maximum time index (exclusive) in this implosion. Chosen to be exclusive so one can::

            for it in range(foo.it_min(), foo.it_max()):
                bar(it)

        :returns: The time index
        """
        pass
    
    @abstractmethod
    def ir_min(self):
        """Minimum radius (inclusive) in this implosion.

        :returns: The lowest acceptable value in `ir` (typically will be 0)
        """
        pass

    @abstractmethod
    def ir_fuel(self):
        """Outer radius of fuel material.

        :returns: A radial index
        """
        pass

    @abstractmethod
    def ir_max(self):
        """Maximum radial index (exclusive) in this implosion. Chosen to be exclusive so one can::

            for ir in range(foo.ir_min(), foo.ir_max()):
                bar(ir)

        :returns: The radial index"""
        pass

    @abstractmethod
    def t(self, it):
        """Convert indices to real time.

        :param it: The temporal index
        :returns: real time [s]
        """
        pass

    @abstractmethod
    def dt(self, it, ir=None):
        """Get the post-processor time step.

        :param it: Time index (may be ignored if implosions have a constant time step)
        :param ir: Optional since dt is the same for all spatial zones. However, this gives the option of getting a
            2-D array as the return value, which is convenient for some calculations.
        :returns: The delta in time between it and it+1 [s]
        """
        pass

    @abstractmethod
    def r(self, it, ir):
        """Get physical radius for a zone.

        :param it: The temporal index
        :param ir: The spatial index
        :returns: r [cm]
        """
        pass

    def checkIndices(self, it, ir):
        """Validate it and ir.

        :param it: A temporal index or tuple
        :param ir: A spatial index or tuple
        :returns: True if it and ir are acceptable, False otherwise.
        """
        assert np.isscalar(it) or (isinstance(it, tuple) and len(it)==2)
        assert np.isscalar(ir) or (isinstance(ir, tuple) and len(ir)==2)

        # Handle scalar arguments:
        if np.isscalar(it) and np.isscalar(ir):
            return (it >= self.it_min and it < self.it_max() and ir >= self.ir_min() and ir < self.ir_max())

        # Handle ranges:
        # First, if one is scalar while the other is a tuple, fix it:
        if np.isscalar(it):
            it = np.asarray([it,it+1])
        if np.isscalar(ir):
            ir = np.asarray([ir,ir+1])

        # loop over all it, ir:
        for i in it:
            if i < self.it_min() or i >= self.it_max():
                return False
        for j in ir:
            if j < self.ir_min() or j >= self.ir_max():
                return False
        # If both loops complete, all indices are OK:
        return True

    def __internalIndex__(self, it, ir):
        """Internal validation of indices. Differs from above in return.
        Fails an assertion if either index is unacceptable.

        :param it: A temporal index or tuple
        :param ir: A spatial index or tuple
        :returns: modified it and ir (fixes cases where only one is a tuple).
        """
        assert np.isscalar(it) or (isinstance(it, tuple) and len(it)==2)
        assert np.isscalar(ir) or (isinstance(ir, tuple) and len(ir)==2)

        # Handle scalar arguments:
        if np.isscalar(it) and np.isscalar(ir):
            return it, ir

        # Handle ranges:
        # First, if one is scalar while the other is a tuple, fix it:
        if np.isscalar(it):
            it = np.asarray([it,it+1])
        if np.isscalar(ir):
            ir = np.asarray([ir,ir+1])
        return it, ir

    # ----------------------------------------
    #           Material info
    # ----------------------------------------
    @abstractmethod
    def IonA(self, it, ir):
        """Masses for all ions in a zone.

        :param it: The temporal index
        :param ir: The spatial index
        :returns: A :py:class:`numpy.ndarray` of ion masses [AMU]
        """
        pass

    @abstractmethod
    def IonZ(self, it, ir):
        """Ion atomic numbers for ions in a zone.

        :param it: The temporal index
        :param ir: The spatial index
        :returns: A :py:class:`numpy.ndarray` of ion atomic numbers [e]
        """
        pass

    @abstractmethod
    def IonF(self, it, ir):
        """Ion fractions in a zone. Each fraction is between 0 and 1 (inclusive). The sum of all fractions is 1.

        :param it: The temporal index
        :param ir: The spatial index
        :returns: A :py:class:`numpy.ndarray` of ion fractions
        """
        pass

    @abstractmethod
    def Abar(self, it, ir):
        """Average ion mass in a zone.

        :param it: The temporal index
        :param ir: The spatial index
        :returns: Average ion mass [AMU]
        """
        pass

    @abstractmethod
    def Zbar(self, it, ir):
        """Average ion Z.

        :param it: The temporal index
        :param ir: The spatial index
        :returns: Average ion atomic number [e]
        """
        pass

    @abstractmethod
    def f(self, it, ir, A, Z):
        """Get fraction of a specified ion in a zone.

        :param it: The temporal index
        :param ir: The spatial index
        :param A: The ion mass you're interested in (scalar)
        :param Z: The ion atomic number you're interested in (scalar)
        :returns: Fraction of that ion in a zone
        """
        pass

    # ----------------------------------------
    #           Fusion info
    # ----------------------------------------
    def calcYield(self, it, ir, rxn):
        """Calculate the yield produced in a zone for a given fusion reaction.

        :param it: The temporal index
        :param ir: The spatial index
        :param rxn: The `Reaction` to use for fusion burn. Pass the class.
        :return: Number of reactions per zone per time step
        """
        assert rxn in fusion.allReactions()
        it, ir = self.__internalIndex__(it, ir)

        # Check to see if the yield is pre-computed:
        if not rxn.name() in self.yieldData:
            self.__calcYield__(rxn)

        if np.isscalar(it) and np.isscalar(ir):
            return self.yieldData[it,ir]
        return self.yieldData[rxn.name()][it[0]:it[1], ir[0]:ir[1]]

    def __calcYield__(self, rxn):
        """Helper function for calculating fusion yields."""
        dblcount = 1
        if (rxn.reactantA[0] == rxn.reactantA[1]) and (rxn.reactantZ[0] == rxn.reactantZ[1]):
            dblcount = 2 #account for factor of 2 if reactants are identical

        it = (self.it_min(), self.it_max())
        ir = (self.ir_min(), self.ir_max())
        f1 = self.f(it, ir, rxn.reactantA[0], rxn.reactantZ[0])
        f2 = self.f(it, ir, rxn.reactantA[1], rxn.reactantZ[1])
        sigmav = rxn.reactivity(self.Ti(it,ir))

        Y = sigmav * np.power(self.ni(it,ir), 2) * (f1*f2/dblcount) * self.vol(it,ir) * self.dt(it,ir)

        self.yieldData[rxn.name()] = Y

    # ----------------------------------------
    #      Stuff needed for pickling
    # ----------------------------------------
    def __getstate__(self):
        """Get the current state of this object as a `dict`.
        If there are members of the object that cannot be pickled, such as
        open files or GUI windows, they must be removed from the returned `dict`.
        """
        state = self.__dict__.copy()
        return state

    def __setstate__(self, state):
        """Set the state of this object from a given `dict`"""
        # Restore instance attributes (i.e., filename and lineno).
        self.__dict__.update(state)


def allImplosions():
    """Get a list containing all implemented implosions. """
    temp = Implosion.__subclasses__() + [g for s in Implosion.__subclasses__()
                                   for g in s.__subclasses__()]
    temp2 = []
    for t in temp:
        if not inspect.isabstract(t):
            temp2.append(t)
    return temp2