"""Define useful physical constants, using CGS units.

:author: Alex Zylstra
:date: 2014-01-03
"""

__author__ = 'Alex Zylstra'
__date__ = '2014-01-03'
__version__ = '0.1.0'


#: fundamental charge [statC]
e = 4.803e-10
#: speed of light [cm/s]
c = 3.0e10
#: electron mass [g]
me = 9.109e-28
#: proton mass [g]
mp = 1.672e-24
#: Boltzmann constant [erg/K]
kB = 1.381e-16
#: Planck's constant [erg*s]
h = 6.626e-27
#: Planck's constant [erg*s]
hbar = 1.0546e-27
#: Avagadro's number
Na = 6.022e23
#: Atomic Mass Unit [g]
AMU = 1.66053873e-24
#: Convert from erg to eV
Erg2eV = 1./6.2415e11
#: Convert pressure from CGS (Bayre) to GBar
P2GBAR = 1.0e-10/1.01325e5