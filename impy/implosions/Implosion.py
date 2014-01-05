""" Python-based abstract representation of an implosion. All implosion types must implement these methods.

:author: Alex Zylstra
:date: 2014-01-05
"""

__author__ = 'Alex Zylstra'
__date__ = '2014-01-05'
__version__ = '1.0.0'

from abc import ABCMeta, abstractmethod
import numpy as np
import inspect

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

    The behavior of all functions taking indices (time and radius, `it` and `ir` respectively) is as follows.

    Both indices can be single integers, in which case the return type will be a scalar number unless otherwise noted::

        >>> # foo is an Implosion
        >>> bar = foo.ni(5,5)
        >>> print(bar)
        5.

    Both indices can be numpy arrays (:py:class:`numpy.ndarray`) containing several desired indices to sample at.
    In this case a 2-D :py:class:`numpy.ndarray` is returned where the first index corresponds to the time indices and
    the second axis corresponds to the radial indices::

        >>> # foo is an Implosion
        >>> it = np.asarray([1,2,3,4])
        >>> ir = np.asarray([8,9])
        >>> data = foo.ni(it, ir)
        >>> type(data)
        <class 'numpy.ndarray'>
        >>> data.shape
        (4,2)
        >>> data[0,0] == foo.ni(1,8)
        True

    One note is that there are several functions for material composition (`IonA`, `IonF`, `IonZ`) which return
    :py:class:`numpy.ndarray` for a single zone. If called with :py:class:`numpy.ndarray` for `it` and `ir` then
    these functions return a 3-D :py:class:`numpy.ndarray`.

    :author: Alex Zylstra
    :date: 2014-01-05
    """
    
    # ----------------------------------------
    #           Generic methods
    # ----------------------------------------
    @abstractmethod
    def __init__(self, type='GUI', args=''):
        """Construct a new implosion."""
        pass

    @classmethod
    @abstractmethod
    def getFileTypes(cls):
        """Get a list containing extensions of file types supported by this implosion."""
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
    def abort(self):
        """Signal that the implosion generation should be interrupted."""
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
    def dt(self, it):
        """Get the post-processor time step.

        :param it: Time index (may be ignored if implosions have a constant time step)
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

        :returns: True if it and ir are acceptable, False otherwise.
        """
        # If multiple indices:
        if isinstance(it, np.ndarray) and isinstance(ir, np.ndarray):
            # loop over all it, ir:
            for i in it:
                if i < self.it_min() or i >= self.it_max():
                    return False
            for j in ir:
                if j < self.ir_min() or j >= self.ir_max():
                    return False
            # If both loops complete, all indices are OK:
            return True

        # Scalars:
        assert isinstance(it, int)
        assert isinstance(ir, int)
        return (it >= self.it_min and it < self.it_max() and ir >= self.ir_min() and ir < self.ir_max())

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


def allImplosions():
    """Get a list containing all implemented implosions. """
    temp = Implosion.__subclasses__() + [g for s in Implosion.__subclasses__()
                                   for g in s.__subclasses__()]
    temp2 = []
    for t in temp:
        if not inspect.isabstract(t):
            temp2.append(t)
    return temp2