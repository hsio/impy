# HYADES wrapper class
# Reads in a .nc file and does interpolation
# to support the Python post-processor modules
# Implements all functions called in Implosion.py
# This code generally uses CGS units
#
# Author: Alex Zylstra
# Date: 2013/11/05

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
    
    #Raw read from HYADES
    NumRegions = 0
    Regions = []
    
    #raw data
    t_raw = []
    r_raw = []
    rcm = []
    rcm2 = []
    TiRaw = []
    TeRaw = []
    niRaw = []
    neRaw = []
    uRaw = []
    VolRaw = []
    
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
    
    #some time scales in the problem
    t0 = 0
    tc = 0
    itc = 0
    # radial index corresponding to fuel limit:
    iFuel = 0
    
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
        self.find_tc()
        
    # ------------------------------------
    # Read in the HYADES file
    # ------------------------------------
    def readHYADES(self):
        """Read in the HYADES file."""
        # raw hydro info
        self.t_raw = numpy.asarray( self.file.variables['DumpTimes'][:] )
        self.r_raw = numpy.asarray( self.file.variables['R'][:] )
        self.rcm = numpy.asarray( self.file.variables['Rcm'][:] )
        
        # fix weirdness with rcm
        temp = []
        for i in range(len(self.rcm)):
            temp.append( self.rcm[i][1:len(self.rcm[i])-1] )
        self.rcm2 = numpy.asarray( temp )

        self.TiRaw = numpy.asarray( self.file.variables['Ti'][:] )
        self.TeRaw = numpy.asarray( self.file.variables['Te'][:] )
        self.niRaw = numpy.asarray( self.file.variables['Deni'][:] )
        self.neRaw = numpy.asarray( self.file.variables['Dene'][:] )       
        self.uRaw = numpy.asarray( self.file.variables['Ucm'][:] )
        self.VolRaw = numpy.asarray( self.file.variables['Vol'][:] )
        
        # raw material info
        self.RegNums = numpy.asarray( self.file.variables['RegNums'][:] )
        self.NumMatsReg = numpy.asarray( self.file.variables['NumMatsReg'][:] )
        self.AtmFrc = numpy.asarray( self.file.variables['AtmFrc'][:] )
        self.AtmNum = numpy.asarray( self.file.variables['AtmNum'][:] )
        self.AtmWgt = numpy.asarray( self.file.variables['AtmWgt'][:] )

        
        # Work on material definitions
        self.NumRegions = max(self.RegNums) # how many regions we have to deal with
        
    # ------------------------------------
    # Hydro variables (interpolation)
    # ------------------------------------
    def Ti(self, ir, it):
        """With ir and it radial and temporal indices, Returns Ti in keV."""
        return self.TiRaw[it][ir]
    def Te(self, ir, it):
        """Electron temperature with ir and it radial and temporal indices, Returns Te in keV."""
        return self.TeRaw[it][ir]
    def ni(self, ir, it):
        """Ion number density ni(r,t) with ir and it radial and temporal indices, Returns ni in 1/cm^3."""
        return self.niRaw[it][ir]
    def ne(self, ir, it):
        """Electron number density ne(r,t) with ir and it radial and temporal indices, Returns ne in 1/cm^3."""
        return self.neRaw[it][ir]
    def u(self, ir, it):
        """Fluid velocity u(r,t) with ir and it radial and temporal indices,. Returns u in um/ns."""
        return self.uRaw[it][ir]
    def c(self, ir, it):
        """Sound speed c(r,t) with ir and it radial and temporal indices, Returns c in um/ns."""
        return math.sqrt( (5/3) * self.P(ir,it) * pow(10,15) / self.rho(ir,it) )
    def rho(self, ir, it):
        """Density rho(r,t) with ir and it radial and temporal indices, Returns rho in g/cm3."""
        return self.ni(ir,it) * mp * numpy.dot( self.IonF(ir,it) , self.IonA(ir,it) )
    def P(self, ir, it):
        """Pressure P(r,t) with ir and it radial and temporal indices, Returns P in GBar."""
        return (kB * (self.ni(ir,it)*self.Ti(ir,it) + self.ne(ir,it)*self.Te(ir,it)) * 11600*1000) * pow(10,-15)
    def vol(self, ir, it):
        """Volume of zone index ir at time index it."""
        return self.VolRaw[it][ir]
    def rhoR(self, ir, it):
        """Areal density of zone index ir at time index it. Returns rhoR in g/cm2"""
        return self.rho(ir,it)*(self.r_raw[it][ir+1]-self.r_raw[it][ir])
        
    # ------------------------------------
    # Find timescales in the problem
    # ------------------------------------
    def find_tc(self):
        """Find the shock collapse time."""
        tc = self.t_raw[0]
        prevTi = self.Ti(0, 0)
        Ti = self.Ti(0, 1)
        i = 1
        #iterate until the ion temp changes by more than 1keV
        while (i+1) < len(self.t_raw):
            i += 1
            prevTi = Ti
            Ti = self.Ti(0, i)
            if Ti > 0. and (Ti-prevTi) >= 1:
                tc = self.t_raw[i]
                return tc
        self.tc = tc
        self.itc = max(0,i-1)
            
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
            for i in range(len(self.t_raw)):
                rList.append( self.rcm2[i][iz] )
                times.append( self.t_raw[i]*1e9 )
            self.rRegionInt.append( scipy.interpolate.UnivariateSpline(times, rList) )

        # find where the fuel is in the problem:
        iFuel = 0
        for j in range(self.ir_max()): #find out how many regions actually contain stuff
            if self.fuel(self.IonA(j,0),self.IonZ(j,0)):
                iFuel = max(iFuel, j)
        self.iFuel = iFuel

    def rRegion(self, it, ir):
        """Region boundary position at time index it. Returns r in cm."""
        return self.rRegionInt[ir](self.it)[()]
    def RegionIdent(self, ir):
        """For radial index ir, returns material region number."""
        return self.RegNums[ir]
    def fuel(self, A, Z):
        """Check to see if this definition of Z contains something that could be fuel (Z=1,2)."""
        FuelA = [2,3,3]
        FuelZ = [1,1,2]
        for i in range(len(A)):
            for j in range(len(FuelA)):
                if FuelA[j] == A[i] and FuelZ[j] == Z[i]:
                    return True
        return False
        
            
    # ------------------------------------
    # Implementation of remaining Implosion functions for radius/time
    # ------------------------------------
    # required time limits in the problem
    def it_min(self):
        """Minimum time index for post-proc calculations (s)."""
        return 0
    def it_tc(self):
        """Time index corresponding to right before shock coalescence."""
        return self.itc
    def it_max(self):
        """Maximum time index for post-proc calculations (s)."""
        return len(self.t_raw)
    # required length limits in the problem
    def ir_min(self):
        """Minimum radial index for post-proc calculations."""
        return 0
    def ir_fuel(self):
        """Maximum radial index of fuel material."""
        return self.iFuel
    def ir_max(self):
        """Maximum radial index for post-proc calculations."""
        return len(self.TiRaw[0])
    def t(self, it):
        """Real time in s corresponding to index value it."""
        return self.t_raw[it]
    def dt(self):
        """Get the post-processor time step in s."""
        return self.t_raw[1]-self.t_raw[0]
    def r(self ,ir, it):
        """Real radius in cm corresponding to index value ir."""
        return self.rcm2[it][ir]

    # ------------------------------------
    # Implementation of material info
    # ------------------------------------
    # required material composition info
    def IonA(self, ir, it):
        """List of AMU masses for all ions at radial, temporal indices ir and it."""
        return self.AtmWgt[max(self.RegNums[ir]-1, 0)]  # fix issue with first zone
    def IonZ(self, ir, it):
        """List of ion Z for all ions at radial, temporal indices ir and it."""
        return self.AtmNum[max(self.RegNums[ir]-1, 0)]  # fix issue with first zone
    def IonF(self, ir, it):
        """List of ion relative populations."""
        return self.AtmFrc[max(self.RegNums[ir]-1, 0)]  # fix issue with first zone
    def Abar(self, ir, it):
        """Average ion A, at radial, temporal indices ir and it."""
        return numpy.dot( self.IonA(ir,it) , self.IonF(ir,it) )
    def Zbar(self, ir, it):
        """Average ion Z, at radial, temporal indices ir and it."""
        return numpy.dot( self.IonZ(ir,it) , self.IonF(ir,it) )
    def f(self, ir, it, A, Z):
        """Get composition fraction for fuel ion with A and Z at radius/time indices ir and it."""
        IonA = self.IonA(ir,it)
        IonF = self.IonF(ir,it)
        IonZ = self.IonZ(ir,it)
        region = self.RegNums[ir]
        for i in range( len(IonF) ):
            if IonA[i] == A and IonZ[i] == Z:
                return IonF[i]
        return 0
