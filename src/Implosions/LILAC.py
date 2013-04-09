# LILAC wrapper class
# Reads in a .lilac file and does interpolation
# to support the Python post-processor modules
# Implements all functions called in Implosion.py
# This code generally uses CGS units
#
# Author: Alex Zylstra
# Date: 2012/08/13

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
    #fuel info
    #names of valid fuel ions, if they are present,
    #and definitions of their A and Z
    FuelIonNames = ['H','D','T','3He']
    #Fuel composition info
    A = []
    Z = []
    F = []
    Abar = []
    Zbar = []
    
    #Raw read from LILAC
    NumRegions = 0
    Regions = []
    rRegionMax = []
    rRegionInt = []
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
    #some physical scales
    rMax = []
    rMaxInt = 0
    
    # ------------------------------------
    # Initialization
    # ------------------------------------
    def __init__(self):
        """Initialization. 'file' is the path to the LILAC output."""
        self.filename = input("LILAC file: ")
        self.file = open(self.filename,'r')
        self.readLILAC()
        self.tc = self.find_tc()
        self.rRegionInit()
        
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
            #add to regions
            self.Regions.append( [int(line[len(line)-2]), int(line[len(line)-1]) ] ) #region boundaries
            # update materials
            self.MatIdent(line[1])
            self.NumRegions += 1
            self.rRegionMax.append([])
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
                r = float( line[self.rIndex] ) #radius
                # read in all variables of interest
                self.TiRaw.append([ r, t, float( line[self.TiIndex] )*self.TUnitConv ])
                self.TeRaw.append([ r, t, float( line[self.TeIndex] )*self.TUnitConv ])
                self.niRaw.append([ r, t, float( line[self.niIndex] ) ])
                self.uRaw.append([ r, t, float( line[self.uIndex] ) ])
                self.PRaw.append([ r, t, float( line[self.PIndex] )*self.PUnitConv ])
                # keep track of region boundaries as a function of time
                for i in range(len(self.Regions)):
                    if zone+2 == self.Regions[i][1]:
                        self.rRegionMax[i].append( [r , t] )
                line = self.cleanRead(dataReader)
                zone += 1
        
        #Set up interpolation
        #Ti
        rt = []
        z = numpy.array( [] , float)
        for i in self.TiRaw:
            rt.append( [i[0],i[1]] )
            z = numpy.append(z, i[2])
        self.TiInt = scipy.interpolate.NearestNDInterpolator(rt, z)
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
    def uLab(self, r, t):
        """Fluid velocity u(r,t) with r in cm and t in s. Returns u in um/ns."""
        t = t * 1e9 #convert from s to ns
        return uInt([[r,t]])[0]
    def u(self, r, t):
        """Fluid velocity u(r,t) with r in cm and t in s. Returns u in um/ns."""
        return self.uLab(r,t)
    def c(self, r, t):
        """Sound speed c(r,t) with r in cm and t in s Returns c in um/ns."""
        return math.sqrt( self.gamma*self.P(r,t) / self.rho(r,t) )
    def rho(self, r, t):
        """Density rho(r,t) with r in cm and t in s Returns rho in g/cm3."""
        return self.mp*self.FuelA*self.ni(r,t)
    def T(self, r, t):
        """One-fluid hydro temperature T(r,t) with r in cm and t in s Returns T in keV."""
        return self.Ti(r,t)
    def P(self, r, t):
        """One-fluid hydro pressure P(r,t) with r in cm and t in s Returns P in GBar."""
        t = t * 1e9 #convert from s to ns
        return self.PInt([[r,t]])[0]
    def Ti(self, r, t):
        """with r in cm and t in s Returns Ti in keV."""
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
        #iterate until the ion temp changes by more than 1keV
        while (i+1) < len(self.time):
            i += 1
            prevTi = Ti
            Ti = self.Ti(5e-4, self.time[i]*1e-9)
            if Ti > 0. and (Ti-prevTi) >= 1:
                tc = self.time[i]*1e-9
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
    def rRegionInit(self):
        """Helper function to do initial region boundary position interpolation."""
        for ir in range(self.NumRegions):
            self.rRegionInt.append(0)
            times = []
            rList = []
            #populate lists from 2D data
            for i in self.rRegionMax[ir]:
                rList.append( i[0] )
                times.append( i[1] )
            #self.rRegionInt[ir] = scipy.interpolate.interp1d(times, rList, kind='cubic') #also works, slower
            self.rRegionInt[ir] = scipy.interpolate.UnivariateSpline(times, rList)
    def rRegion(self, t, ir):
        """Region boundary position at time t (s). Returns r in cm."""
        return self.rRegionInt[ir](t*1e9)[()]
    def RegionIdent(self, r, t):
        """For r in cm and t in s, returns material region number."""
        for ir in range(self.NumRegions):
            if r <= self.rRegion(t,ir):
                return ir
        return self.NumRegions-1
        
            
    # ------------------------------------
    # Implementation of remaining Implosion functions
    # ------------------------------------
    # required time limits in the problem
    def tmin(self):
        """Minimum time for post-proc calculations (s)."""
        return max(self.time[0]*1e-9, self.tc-1e-10)
    def tmax(self):
        """Maximum time for post-proc calculations (s)."""
        return self.time[ len(self.time) - 1 ]*1e-9
    # required length limits in the problem
    def rmin(self, t):
        """Minimum radius for post-proc calculations at time t in s."""
        return 0
    def rmax(self, t):
        """Maximum radius for post-proc calculations at time t in s."""
        iFuel = 0
        for j in range(self.NumRegions): #find out how many regions actually contain fuel
            if len(self.A[j]) > 0:
                iFuel = max(iFuel, j)
        return self.rRegion(t, iFuel)
        
    # required material composition info
    def IonA(self, r, t):
        """List of AMU masses for all ions at r in cm and t in s."""
        return self.A[self.RegionIdent(r,t)]
    def IonZ(self, r, t):
        """List of ion Z for all ions at r in cm and t in s."""
        return self.Z[self.RegionIdent(r,t)]
    def IonF(self, r, t):
        """List of ion relative populations."""
        return self.F[self.RegionIdent(r,t)]
        
    # ------------------------------------
    # Handle LILAC material info
    # ------------------------------------
    def MatIdent(self, Num):
        """Read in LILAC material definitions from CSV file, looking for ID # Num."""
        MatFile = open("Resources/LILAC_Materials.csv",'r')
        if not MatFile.readable(): #Sanity check
            print("Error reading LILAC material file!")
            return
            
        MatFile.readline() #discard header
        
        dataReader = csv.reader(MatFile , delimiter=',')
        
        A = []
        Z = []
        F = []
        Abar = 0
        Zbar = 0
        for row in dataReader:
            if int(row[0]) == int(Num):
                if float(row[1]) > 0: # H is present
                    A.append(1)
                    Z.append(1)
                    F.append(float(row[1]))
                if float(row[2]) > 0: # D is present
                    A.append(2)
                    Z.append(1)
                    F.append(float(row[2]))
                if float(row[3]) > 0: # T is present
                    A.append(3)
                    Z.append(1)
                    F.append(float(row[3]))
                if float(row[4]) > 0: # 3He is present
                    A.append(3)
                    Z.append(2)
                    F.append(float(row[4]))
                Abar = float(row[5])
                Zbar = float(row[6])
        
        #sanity check
        if Abar == 0 or Zbar == 0:
            print("ERROR: material not found!")
        
        self.A.append(A)
        self.Z.append(Z)
        self.F.append(F)
        self.Abar.append(Abar)
        self.Zbar.append(Zbar)
        return
