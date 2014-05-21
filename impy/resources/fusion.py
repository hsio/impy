""" Fusion reaction cross-sections and reactivities.

:author: Alex Zylstra
:date: 2014-01-03
"""

__author__ = 'Alex Zylstra'
__date__ = '2014-01-03'
__version__ = '0.1.0'

import numpy as np
from impy.resources.constants import *
import scipy
import numpy as np
import os, csv
from abc import ABCMeta, abstractmethod

def T9(Ti):
    """Convert ion temperature from keV to T9.

    :param Ti: Ion temperature [keV]
    """
    return Ti*1.16045e-2

def Eg(Ti,Z1,Z2,A1,A2):
    """Gamow energy [keV].

    :param Ti: Ion temperature [keV]
    :param Z1: Reactant 1 charge [e]
    :param Z2: Reactant 2 charge [e]
    :param A1: Reactant 1 mass [AMU]
    :param A2: Reactant 2 mass [AMU]
    """
    xi = 6.2696*pow(Z1*Z2,2/3)*pow(A1*A2/(A1+A2),1/3)
    return xi * pow(Ti,2/3)

class Reaction(metaclass=ABCMeta):
    """Abstract representation of a reaction. All fusion reactions derive from this class.

    :author: Alex Zylstra
    :date: 2014-01-03
    """

    #: Tuple containing reactant masses [AMU]
    reactantA = (0,0)

    #: Tuple containing reactant atomic numbers [e]
    reactantZ = (0,0)

    @classmethod
    @abstractmethod
    def name(cls):
        """String containing a name for this reaction, e.g. 'DT'"""
        pass

    @classmethod
    @abstractmethod
    def reactivity(cls, Ti):
        """Calculate reactivity [cm^3/s].

        :param Ti: Temperature [keV]
        :raises: :py:exc:`ValueError`
        """
        pass

    @classmethod
    @abstractmethod
    def crossSection(cls, En):
        """Calculate cross section [cm^2].

        :param En: CM energy [keV]
        :raises: :py:exc:`ValueError`
        """
        pass

class BoschHale(Reaction):
    """Implement reactions based off of Bosch and Hale. Subclasses must replace class variables.
    This class remains abstract.
    Reference: H.-S. Bosch, G.M. Hale, Nuclear Fusion Vol 32, No. 4, 611 (1992)

    :author: Alex Zylstra
    :date: 2014-01-03
    """

    # Values for the cross sections:
    BG = 0
    A1 = 0
    A2 = 0
    A3 = 0
    A4 = 0
    A5 = 0
    B1 = 0
    B2 = 0
    B3 = 0
    B4 = 0
    energyMin = 0
    energyMax = 0

    # Values for the reactivity fits:
    mc2 = 0
    C1 = 0
    C2 = 0
    C3 = 0
    C4 = 0
    C5 = 0
    C6 = 0
    C7 = 0
    TiMin = 0
    TiMax = 0

    @classmethod
    def reactivity(cls, Ti):
        """Calculate reactivity [cm^3/s].

        :param Ti: Temperature [keV]
        :raises: :py:exc:`ValueError`
        """
        #TODO: general sanity checks for scalars and arrays
        #if Ti < cls.TiMin or Ti > cls.TiMax:
        #    raise ValueError('Temperature outside of valid range for ' + cls.name() + ' at ' + str(Ti))

        theta = Ti / (1-Ti*(cls.C2+Ti*(cls.C4+Ti*cls.C6))/(1+Ti*(cls.C3+Ti*(cls.C5+Ti*cls.C7))))
        xi = np.power(np.power(cls.BG,2)/(4*theta), 1/3)
        return cls.C1*theta*np.sqrt(xi/(cls.mc2*np.power(Ti,3)))*np.exp(-3*xi)

    @classmethod
    def crossSection(cls, En):
        """Calculate cross section [cm^2].

        :param En: CM energy [keV]
        :raises: :py:exc:`ValueError`
        """
        #TODO: general sanity checks for scalars and arrays
        #if En < cls.energyMin or En > cls.energyMax:
        #    raise ValueError('Energy outside of valid range for ' + cls.name() + ' at ' + str(En))

        SE = (cls.A1 + En*(cls.A2 + En*(cls.A3 + En*cls.A4))) / (1 + En*(cls.B1+En*(cls.B2+En*(cls.B3+En*cls.B4))))
        return 1e-27 * SE / ( En*np.exp(cls.BG/np.sqrt(En)) )

class DT(BoschHale):
    """D+T fusion, i.e. T(D,n)4He."""
    reactantA = (2,3)
    reactantZ = (1,1)
    @classmethod
    def name(cls):
        return 'DT'
    C1 = 1.17302e-9
    C2 = 1.51361e-2
    C3 = 7.51886e-2
    C4 = 4.60643e-3
    C5 = 1.35e-2
    C6 = -1.0675e-4
    C7 = 1.366e-5
    BG = 34.3827
    mc2 = 1124656
    A1 = 6.927e4
    A2 = 7.454e8
    A3 = 2.050e6
    A4 = 5.2002e4
    B1 = 6.38e1
    B2 = -9.95e-1
    B3 = 6.981e-5
    B4 = 1.728e-4
    BG = 34.3827

    energyMin = 0.5
    energyMax = 550
    TiMin = 0.2
    TiMax = 100

class DDn(BoschHale):
    """DDn fusion, i.e. D(D,n)3He."""
    reactantA = (2,2)
    reactantZ = (1,1)
    @classmethod
    def name(cls):
        return 'DDn'
    C1 = 5.43360e-12
    C2 = 5.85778e-3
    C3 = 7.68222e-3
    C4 = 0
    C5 = -2.964e-6
    C6 = 0
    C7 = 0
    BG = 31.3970
    mc2 = 937814
    A1 = 5.3701e4
    A2 = 3.3027e2
    A3 = -1.2706e-1
    A4 = 2.9327e-5
    A5 = -2.5151e-9

    energyMin = 0.5
    energyMax = 4900
    TiMin = 0.2
    TiMax = 100

class DDp(BoschHale):
    """DDp fusion, i.e. D(D,p)T."""
    reactantA = (2,2)
    reactantZ = (1,1)
    @classmethod
    def name(cls):
        return 'DDp'
    C1 = 5.65718e-12
    C2 = 3.41267e-3
    C3 = 1.99167e-3
    C4 = 0.
    C5 = 1.05060e-5
    C6 = 0.
    C7 = 0.
    BG = 31.3970
    mc2 = 937814
    A1 = 5.5576e4
    A2 = 2.1054e2
    A3 = -3.2638e-2
    A4 = 1.4987e-6
    A5 = 1.8181e-10

    energyMin = 0.5
    energyMax = 5000
    TiMin = 0.2
    TiMax = 100


class D3He(BoschHale):
    """D3He fusion, i.e. 3He(D,p)4He."""
    reactantA = (2,3)
    reactantZ = (1,2)
    @classmethod
    def name(cls):
        return 'D3He'
    C1 = 5.51036e-10
    C2 = 6.41918e-3
    C3 = -2.02896e-3
    C4 = -1.91080e-5
    C5 = 1.35776e-4
    C6 = 0
    C7 = 0
    BG = 68.7508
    mc2 = 1124572
    A1 = 5.7501e6
    A2 = 2.5226e3
    A3 = 4.5566e1
    B1 = -3.1995e-3
    B2 = -8.5530e-6
    B3 = 5.9014e-8

    energyMin = 0.3
    energyMax = 900
    TiMin = 0.5
    TiMax = 190

#TODO: T3He
#TODO: TT
#TODO: 3He3He
#TODO: HD
#TODO: p11B
#TODO: p15N


# For ease in iterating:
import inspect
def allReactions():
    """Get a list containing all classes which fully implement Reaction. For example::

        for r in allReactions():
            print(r.crossSection(50))

    prints the cross section value at 50 keV Ecm for all reactions defined.
    """
    temp = Reaction.__subclasses__() + [g for s in Reaction.__subclasses__()
                                   for g in s.__subclasses__()]
    temp2 = []
    for t in temp:
        if not inspect.isabstract(t):
            temp2.append(t)
    return temp2


def fuel(A, Z):
    """Check to see if this definition of ion composition contains something that could be fuel,
    i.e. it will produce yield for some defined reaction.

    :param A: A list of ion atomic mass [AMU]
    :param Z: A list of ion atomic numbers [e]

    A and Z must be either python lists or :py:class:`numpy.ndarray`
    """
    # sanity:
    assert isinstance(A, list) or isinstance(A, np.ndarray)
    assert isinstance(Z, list) or isinstance(Z, np.ndarray)

    # Come up with lists of fuel ion pairs:
    FuelA = []
    FuelZ = []
    for r in allReactions():
        FuelA.append(r.reactantA)
        FuelZ.append(r.reactantZ)

    # Now compare what we were given to the lists of stuff that's fuel:
    for i in range(len(FuelA)):
        for j in range(len(A)):
            if FuelA[i][0] == A[j] and FuelA[i][1] == A[j] \
                and FuelZ[i][0] == Z[j] and FuelZ[i][1] == Z[j]:
                return True
    return False