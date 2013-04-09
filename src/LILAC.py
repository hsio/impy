# LILAC wrapper class
# Reads in a .lilac file and does interpolation
# to support the Python post-processor
# This code generally uses CGS units
#
# Author: Alex Zylstra
# Date: 2012/02/12

import math
import numpy
import scipy.optimize
import scipy.integrate
import scipy.interpolate
import csv
import os
import string

class LILAC:
    """A wrapper class for a LILAC simulation."""
    # -----------------------------------------------------------
    # Global variables within this class
    # -----------------------------------------------------------
    # user specified parameters (default values here)
    name = "LILAC"
    filename = ""
    file = 0
    #Gas info
    gamma = 5./3.
    #D info, fraction
    Ion1 = "D"
    Ion1A = 2.
    Ion1Z = 1.
    f1 = 0.5
    #3He info, fraction
    Ion2 = "3He"
    Ion2A = 3.
    Ion2Z = 2.
    f2 = 0.5
    FuelA = Ion1A*f1 + Ion2A*f2
    FuelZ = Ion1Z*f1 + Ion2Z*f2
    
    #Raw read from LILAC
    fuelRegions = []
    time = []
    TiRaw = []
    niRaw = []
    TeRaw = []
    PRaw = []
    uRaw = []
    #Indices for various fields in the raw data
    rIndex = 0
    TiIndex = 0
    niIndex = 0
    TeIndex = 0
    PIndex = 0
    uIndex = 0
    #convert from Dy/cm2 to GBar
    PUnitConv = 1e-15
    #convert from eV to keV
    TUnitConv = 1e-3
    #interpolating functions
    TiInt = 0
    TeInt = 0
    niInt = 0
    PInt = 0
    uInt = 0
    #some time scales in the problem
    t0 = 0
    tc = 0
    #Some physical constants
    mp = 1.672e-24 #g
    
    # ------------------------------------
    # Initialization
    # ------------------------------------
    def __init__(self,fname):
        """Initialization. 'file' is the path to the LILAC output."""
        self.filename = fname
        self.file = open(self.filename,'r')
        self.setup()
        self.readLILAC()
        self.find_tc()
        
    def setup(self):
        """Some setup routines, i.e. read in info from user."""
        self.f1 = float(input("D fraction = "))
        self.f2 = float(input("3He fraction = "))
        #normalization
        self.f1 = self.f1 / (self.f1 + self.f2)
        self.f2 = self.f2 / (self.f1 + self.f2)
        #Calc avg A and Z:
        self.FuelA = self.Ion1A*self.f1 + self.Ion2A*self.f2
        self.FuelZ = self.Ion1Z*self.f1 + self.Ion2Z*self.f2
    
        
    # ------------------------------------
    # Read in the LILAC file
    # ------------------------------------
    def readLILAC(self):
        """Read in the LILAC file."""
        if not self.file.readable():
            print("Error reading LILAC data!")
            return
            
        #read in the name
        self.name = self.file.readline()
        
        # Read in the fuel region info
        dataReader = csv.reader(self.file , delimiter=' ')
        line = self.cleanRead(dataReader)
        while line[0].isdigit():
            #add to fuel zones
            if 'D' in line[2] or 'He' in line[2]:
                self.fuelRegions.append( [line[3], line[4]] )
            line = self.cleanRead(dataReader)
        
        #Read in the header until we get to the start of the data
        #switch to comma delimited
        dataReader = csv.reader(self.file , delimiter=',')
        prevLine = line
        while not 'time' in line[0]:
            prevLine = line
            line = self.cleanRead(dataReader)
        header = prevLine
        #iterate through header, identify columns that we want
        for i in range( len(header) ):
            if "Distance" in header[i]:
                self.rIndex = i
            if ("Ion" or "ion") in header[i] and "eV" in header[i]:
                self.TiIndex = i
            if ("Ion" or "ion") in header[i] and "Density" in header[i]:
                self.niIndex = i
            if ("elec" or "Elec") in header[i] and "eV" in header[i]:
                self.TeIndex = i
            if "Velocity" in header[i]:
                self.uIndex = i
            if "Pressure" in header[i]:
                self.PIndex = i
        
        #Read in the actual data
        dataReader = csv.reader(self.file , delimiter=' ')
        #first time snap was previously read
        line = [ "time=" , line[0].split(' ')[1] ]
        while len(line) > 0 and 'time=' in line[0]:
            #read in one time snapshot:
            t = float( line[1] )
            self.time.append( t )
            line = self.cleanRead(dataReader)
            zone = 0
            while len(line) > 0 and line[0][0].isdigit():
                if self.fuel(zone):
                    r = float( line[self.rIndex] )
                    self.TiRaw.append([ r, t, float( line[self.TiIndex] )*self.TUnitConv ])
                    self.TeRaw.append([ r, t, float( line[self.TeIndex] )*self.TUnitConv ])
                    self.niRaw.append([ r, t, float( line[self.niIndex] ) ])
                    self.uRaw.append([ r, t, float( line[self.uIndex] ) ])
                    self.PRaw.append([ r, t, float( line[self.PIndex] )*self.PUnitConv ])
                line = self.cleanRead(dataReader)
                zone += 1
        
        #Set up interpolation
        #Ti
        rt = []
        z = numpy.array( [] , float)
        for i in self.TiRaw:
            rt.append( [i[0],i[1]] )
            z = numpy.append(z, i[2])
        self.TiInt = scipy.interpolate.CloughTocher2DInterpolator(rt, z, fill_value=0.)
        #Te
        rt = []
        z = numpy.array( [] , float)
        for i in self.TeRaw:
            rt.append( [i[0],i[1]] )
            z = numpy.append(z, i[2])
        self.TeInt = scipy.interpolate.NearestNDInterpolator(rt, z)
        #ni
        rt = []
        z = numpy.array( [] , float)
        for i in self.niRaw:
            rt.append( [i[0],i[1]] )
            z = numpy.append(z, i[2])
        self.niInt = scipy.interpolate.NearestNDInterpolator(rt, z)
        #u
        rt = []
        z = numpy.array( [] , float)
        for i in self.uRaw:
            rt.append( [i[0],i[1]] )
            z = numpy.append(z, i[2])
        self.uInt = scipy.interpolate.NearestNDInterpolator(rt, z)
        #P
        rt = []
        z = numpy.array( [] , float)
        for i in self.PRaw:
            rt.append( [i[0],i[1]] )
            z = numpy.append(z, i[2])
        self.PInt = scipy.interpolate.NearestNDInterpolator(rt, z)
        
          
            
    # ------------------------------------
    # Hydro variables (interpolation)
    # ------------------------------------
    def u(self, r, t):
        """Fluid velocity u(r,t) with r in cm and t in s. Returns u in um/ns."""
        t = t * 1e9 #convert from s to ns
        return uInt([[r,t]])[0]
    def uLab(self, r, t):
        """Fluid velocity u(r,t) with r in cm and t in s. Returns u in um/ns."""
        t = t * 1e9 #convert from s to ns
        return uInt([[r,t]])[0]
    def c(self, r, t):
        """Sound speed c(r,t) with r in cm and t in s Returns c in um/ns."""
        return math.sqrt( self.gamma*self.P(r,t) / self.rho(r,t) )
    def rho(self, r, t):
        """Density rho(r,t) with r in cm and t in s Returns rho in g/cm3."""
        return self.mp*self.FuelZ*self.ni(r,t)
    def T(self, r, t):
        """One-fluid hydro temperature T(r,t) with r in cm and t in s Returns T in keV."""
        return self.Ti(r,t)
    def P(self, r, t):
        """One-fluid hydro pressure P(r,t) with r in cm and t in s Returns P in GBar."""
        t = t * 1e9 #convert from s to ns
        return self.PInt([[r,t]])[0]
    def Ti(self, r, t):
        """Rygg-style ion temperature with defined e-i coupling Ti(r,t) with r in cm and t in s Returns Ti in keV."""
        t = t * 1e9 #convert from s to ns
        return self.TiInt([[r,t]])[0]
    def Te(self, r, t):
        """Rygg-style electron temperature with r in cm and t in s. Returns Te in keV."""
        t = t * 1e9 #convert from s to ns
        return self.TeInt([[r,t]])[0]
    def ni(self, r, t):
        """Ion number density ni(r,t) with r in cm and t in s Returns ni in 1/cm^3."""
        t = t * 1e9 #convert from s to ns
        return self.niInt([[r,t]])[0]
    def ne(self, r, t):
        """Electron number density ne(r,t) with r in cm and t in s Returns ne in 1/cm^3."""
        return self.ni(r,t) * self.FuelZ
    
        
    # ------------------------------------
    # Find timescales in the problem
    # ------------------------------------
    def find_tc(self):
        """Find the shock collapse time."""
        tc = self.time[0]*1e-9
        prevTi = self.Ti(0, self.time[0]*1e-9)
        Ti = self.Ti(5e-4, self.time[1]*1e-9)
        i = 1
        #iterate until the ion temp changes by more than 10%
        while (i+1) < len(self.time):
            i += 1
            prevTi = Ti
            Ti = self.Ti(5e-4, self.time[i]*1e-9)
            if Ti > 0. and (Ti-prevTi)/prevTi >= 1:
                tc = self.time[i]*1e-9
                print(tc)
                return tc
        return tc
            
    # ------------------------------------
    # Helper functions
    # ------------------------------------      
    def cleanRead(self, dataReader):
        """Helper function, reads a line and removes empty elements."""
        line = next(dataReader)
        while line.count('') > 0:
            line.remove('')
        while line.count(' ') > 0:
            line.remove(' ')
        while line.count(' ') > 0:
            line.remove(' ')
        return line
    def fuel(self, i):
        """Check to see if zone i is a fuel zone."""
        ret = 0
        for region in self.fuelRegions:
            ret = ret or (i >= float(region[0]) and i<= float(region[1]))
        return ret
        