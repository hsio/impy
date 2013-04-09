# HYADES wrapper class
# Reads in a .nc file and does interpolation
# to support the Python post-processor modules
# Implements all functions called in Implosion.py
# This code generally uses CGS units
#
# Author: Alex Zylstra
# Date: 2012/08/14

#TODO: c,rho,P

import math
import numpy
import scipy.interpolate
import csv
import os
import string
import scipy.io.netcdf
from Resources.Constants import *

class HYADES:
    """A wrapper class for a HYADES simulation."""
    # -----------------------------------------------------------
    # Global variables within this class
    # -----------------------------------------------------------
    # user specified parameters (default values here)
    filename = ""
    file = 0
    
    #Raw read from LILAC
    NumRegions = 0
    Regions = []
    
    #raw data
    t = []
    r = []
    rcm = []
    rcm2 = []
    TiRaw = []
    TeRaw = []
    niRaw = []
    neRaw = []
    uRaw = []
    
    # raw material info
    RegNums = []
    NumMatsReg = []
    AtmFrc = []
    AtmNum = []
    AtmWgt = []
    NumRegions = 0 # number of defined regions
    # for figuring out region boundaries
    rRegionMax = []
    rRegionInt = []
    
    #interpolating functions
    TiInt = 0
    TeInt = 0
    niInt = 0
    neInt = 0
    uInt = 0
    
    #some time scales in the problem
    t0 = 0
    tc = 0
    
    # ------------------------------------
    # Initialization
    # ------------------------------------
    def __init__(self, fname=None):
        """Initialization. 'file' is the path to the HYADES file."""
        if fname == None:
            self.filename = input("HYADES file: ")
        else:
            self.filename = fname
        self.file = scipy.io.netcdf.netcdf_file(self.filename,'r')
        self.readHYADES()
        self.rRegionInit()
        self.tc = self.find_tc()
        
    # ------------------------------------
    # Read in the HYADES file
    # ------------------------------------
    def readHYADES(self):
        """Read in the HYADES file."""
        # raw hydro info
        self.t = self.file.variables['DumpTimes'][:]
        self.r = self.file.variables['R'][:]
        self.rcm = self.file.variables['Rcm'][:]
        # fix weirdness with rcm
        for i in range(len(self.rcm)):
            self.rcm2.append( self.rcm[i][1:len(self.rcm[i])-1] )
        self.TiRaw = self.file.variables['Ti'][:]
        self.TeRaw = self.file.variables['Te'][:]
        self.niRaw = self.file.variables['Deni'][:]
        self.neRaw = self.file.variables['Dene'][:]       
        self.uRaw = self.file.variables['Ucm'][:]
        
        # raw material info
        self.RegNums = self.file.variables['RegNums'][:]
        self.NumMatsReg = self.file.variables['NumMatsReg'][:]
        self.AtmFrc = self.file.variables['AtmFrc'][:]
        self.AtmNum = self.file.variables['AtmNum'][:]
        self.AtmWgt = self.file.variables['AtmWgt'][:]

        # INTERPOLATION
        # interpolation funcs in ns units to make scipy happy
        # for Ti
        rt = []
        z = numpy.array( [] , float )
        for i in range(len(self.t)): #time indices
            for j in range(len(self.rcm2[i])): #rcm2 indices
                rt.append( [ self.rcm2[i][j] , self.t[i]*1e9 ] )
                z = numpy.append( z, self.TiRaw[i][j] )
        self.TiInt = scipy.interpolate.NearestNDInterpolator(rt, z)
        # for Te
        rt = []
        z = numpy.array( [] , float )
        for i in range(len(self.t)): #time indices
            for j in range(len(self.rcm2[i])): #rcm2 indices
                rt.append( [ self.rcm2[i][j] , self.t[i]*1e9 ] )
                z = numpy.append( z, self.TeRaw[i][j] )
        self.TeInt = scipy.interpolate.NearestNDInterpolator(rt, z)
        # for ni
        rt = []
        z = numpy.array( [] , float )
        for i in range(len(self.t)): #time indices
            for j in range(len(self.rcm2[i])): #rcm2 indices
                rt.append( [ self.rcm2[i][j] , self.t[i]*1e9 ] )
                z = numpy.append( z, self.niRaw[i][j] )
        self.niInt = scipy.interpolate.NearestNDInterpolator(rt, z)
        # for ne
        rt = []
        z = numpy.array( [] , float )
        for i in range(len(self.t)): #time indices
            for j in range(len(self.rcm2[i])): #rcm2 indices
                rt.append( [ self.rcm2[i][j] , self.t[i]*1e9 ] )
                z = numpy.append( z, self.neRaw[i][j] )
        self.neInt = scipy.interpolate.NearestNDInterpolator(rt, z)
        # for u
        rt = []
        z = numpy.array( [] , float )
        for i in range(len(self.t)): #time indices
            for j in range(len(self.rcm2[i])): #rcm2 indices
                rt.append( [ self.rcm2[i][j] , self.t[i]*1e9 ] )
                z = numpy.append( z, self.uRaw[i][j] )
        self.uInt = scipy.interpolate.NearestNDInterpolator(rt, z)
        
        # Work on material definitions
        self.NumRegions = max(self.RegNums) # how many regions we have to deal with
            
    # ------------------------------------
    # Hydro variables (interpolation)
    # ------------------------------------
    def Ti(self, r, t):
        """With r in cm and t in s Returns Ti in keV."""
        return self.TiInt([[r,t*1e9]])[0]
    def Te(self, r, t):
        """Electron temperature with r in cm and t in s Returns Te in keV."""
        return self.TeInt([[r,t*1e9]])[0]
    def ni(self, r, t):
        """Ion number density ni(r,t) with r in cm and t in s Returns ni in 1/cm^3."""
        return self.niInt([[r,t*1e9]])[0]
    def ne(self, r, t):
        """Electron number density ne(r,t) with r in cm and t in s Returns ne in 1/cm^3."""
        return self.neInt([[r,t*1e9]])[0]
    def u(self, r, t):
        """Fluid velocity u(r,t) with r in cm and t in s. Returns u in um/ns."""
        return self.uInt([[r,t*1e9]])[0]
    def c(self, r, t):
        """Sound speed c(r,t) with r in cm and t in s Returns c in um/ns."""
        return math.sqrt( (5/3) * self.P(r,t) * pow(10,15) / self.rho(r,t) )
    def rho(self, r, t):
        """Density rho(r,t) with r in cm and t in s Returns rho in g/cm3."""
        return self.ni(r,t) * mp * numpy.dot( self.IonF(r,t) , self.IonA(r,t) )
    def P(self, r, t):
        """Pressure P(r,t) with r in cm and t in s Returns P in GBar."""
        return (kB * (self.ni(r,t)*self.Ti(r,t) + self.ne(r,t)*self.Te(r,t)) * 11600*1000) * pow(10,-15)
        
    # ------------------------------------
    # Find timescales in the problem
    # ------------------------------------
    def find_tc(self):
        """Find the shock collapse time."""
        tc = self.t[0]
        prevTi = self.Ti(5e-4, self.t[0])
        Ti = self.Ti(5e-4, self.t[1])
        i = 1
        #iterate until the ion temp changes by more than 1keV
        while (i+1) < len(self.t):
            i += 1
            prevTi = Ti
            Ti = self.Ti(5e-4, self.t[i])
            if Ti > 0. and (Ti-prevTi) >= 1:
                tc = self.t[i]
                return tc
        return tc
            
    # ------------------------------------
    # Helper functions
    # ------------------------------------      
    def rRegionInit(self):
        """Helper function to do initial region boundary position interpolation."""
        ZoneBoundaries = []
        for ir in range(1,self.NumRegions+1):
            iz = 0
            while iz < len(self.RegNums) and self.RegNums[iz] <= ir:
                iz += 1
            #subtract 2 b/c RegNums has 0s at start and end, 1 more for while loop overshoot
            ZoneBoundaries.append(iz-3) 
        
        for iz in ZoneBoundaries:
            times = []
            rList = []
            #populate lists from 2D data
            for i in range(len(self.t)):
                rList.append( self.rcm2[i][iz] )
                times.append( self.t[i]*1e9 )
            self.rRegionInt.append( scipy.interpolate.UnivariateSpline(times, rList) )
    def rRegion(self, t, ir):
        """Region boundary position at time t (s). Returns r in cm."""
        return self.rRegionInt[ir](t*1e9)[()]
    def RegionIdent(self, r, t):
        """For r in cm and t in s, returns material region number."""
        for ir in range(self.NumRegions):
            if r <= self.rRegion(t,ir):
                return ir
        return self.NumRegions-1
    def fuel(self, Z):
        """Check to see if this definition of Z contains something that could be fuel (Z=1,2)."""
        return (1 in Z) or (2 in Z)
        
            
    # ------------------------------------
    # Implementation of remaining Implosion functions
    # ------------------------------------
    # required time limits in the problem
    def tmin(self):
        """Minimum time for post-proc calculations (s)."""
        return max(self.t[0], self.tc-1e-10)
    def tmax(self):
        """Maximum time for post-proc calculations (s)."""
        return self.t[ len(self.t) - 1 ]
    # required length limits in the problem
    def rmin(self, t):
        """Minimum radius for post-proc calculations at time t in s."""
        return 0
    def rfuel(self, t):
        """Maximum radius of fuel material at time t in s."""
        iFuel = 0
        for j in range(self.NumRegions): #find out how many regions actually contain stuff
            if self.fuel(self.AtmNum[j]):
                iFuel = max(iFuel, j)
        return self.rRegion(t, iFuel)
    def rmax(self, t):
        """Maximum radius for post-proc calculations at time t in s."""
        return self.rRegion(t, self.NumRegions-1)
        
    # required material composition info
    def IonA(self, r, t):
        """List of AMU masses for all ions at r in cm and t in s."""
        return self.AtmWgt[self.RegionIdent(r,t)]
    def IonZ(self, r, t):
        """List of ion Z for all ions at r in cm and t in s."""
        return self.AtmNum[self.RegionIdent(r,t)]
    def IonF(self, r, t):
        """List of ion relative populations."""
        return self.AtmFrc[self.RegionIdent(r,t)]
    def Abar(self, r, t):
        """Average ion A, at r in cm and t in s."""
        return numpy.dot( self.IonA(r,t) , self.IonF(r,t) )
    def Zbar(self, r, t):
        """Average ion Z, at r in cm and t in s."""
        return numpy.dot( self.IonZ(r,t) , self.IonF(r,t) )
