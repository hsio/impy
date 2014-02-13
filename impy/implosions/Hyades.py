
from impy.implosions.Implosion import Implosion
import os
import scipy
import scipy.io
import numpy as np
from impy.resources import fusion
from impy.resources import constants

class Hyades(Implosion):
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
    def __init__(self, type='GUI', args='', wm=None):
        """Construct a new implosion."""
        super(Hyades, self).__init__()

        # Get file to open based on type:
        if type is 'GUI':
            from tkinter.filedialog import askopenfilename
            FILEOPENOPTIONS = dict(defaultextension='.db',
                           filetypes=[('Hyades netCDF','*.nc')],
                           multiple=False)
            filename = askopenfilename(**FILEOPENOPTIONS)
        elif type is 'File':
            filename = args
        elif type is 'CLI':
            # TODO: CLI args
            filename = input("HYADES file: ")
        else:
            raise ValueError('Type passed to Hyades constructor is invalid: ' + type)

        # Check that file exists, and that extension is correct:
        if not os.path.isfile(filename):
            raise FileNotFoundError('Could not find HYADES file: ' + filename)
        if not filename[-3:] == '.nc':
            raise FileNotFoundError('Invalid file type in HYADES: ' + filename)

        # Set a few instance variables:
        self.filename = filename
        self.runProgress = 0.

    @classmethod
    def getFileTypes(cls):
        """Get a list containing extensions of file types supported by this implosion.
        Must be an array of dicts ready to pass to file dialog, e.g.::

            [('description', '*.ext')]
        """
        return [('Hyades netCDF','*.nc')]

    @classmethod
    def name(cls):
        """Get a string containing a name for this type of implosion."""
        return 'HYADES'

    def info(self):
        """Get a string of information about this specific implosion."""
        #TODO: implement more interesting info
        return 'HYADES file: ' + self.filename

    def generate(self):
        """Run the calculation to generate the implosion data."""
        # Top level construction, uses several helper functions:
        self.runProgress = 0.
        self.__readHYADES__()
        self.__precompute_ions__()
        self.runProgress += 0.05
        self.__rRegionInit__()
        self.runProgress += 0.05
        self.__find_tc__()
        self.runProgress += 0.05

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
            return self.niRaw[it][ir]
        return self.niRaw[it[0]:it[1], ir[0]:ir[1]]

    def ne(self, it, ir):
        """Electron number density ne(t,r).

         :param it: the temporal index
         :param ir: the spatial index
         :returns: ne [1/cm^3]
         """
        it, ir = self.__internalIndex__(it, ir)

        if np.isscalar(it) and np.isscalar(ir):
            return self.neRaw[it][ir]
        return self.neRaw[it[0]:it[1], ir[0]:ir[1]]

    def Ti(self, it, ir):
        """Ion temperature.

         :param it: the temporal index
         :param ir: the spatial index
         :returns: Ti [keV]
         """
        it, ir = self.__internalIndex__(it, ir)

        if np.isscalar(it) and np.isscalar(ir):
            return self.TiRaw[it][ir]
        return self.TiRaw[it[0]:it[1], ir[0]:ir[1]]

    def Te(self, it, ir):
        """Electron temperature.

         :param it: the temporal index
         :param ir: the spatial index
         :returns: Te [keV]
         """
        it, ir = self.__internalIndex__(it, ir)

        if np.isscalar(it) and np.isscalar(ir):
            return self.TeRaw[it][ir]
        return self.TeRaw[it[0]:it[1], ir[0]:ir[1]]

    def u(self, it, ir):
        """Fluid velocity u(t,r).

         :param it: the temporal index
         :param ir: the spatial index
         :returns: u [um/ns]
         """
        it, ir = self.__internalIndex__(it, ir)

        if np.isscalar(it) and np.isscalar(ir):
            return self.uRaw[it][ir]
        return self.uRaw[it[0]:it[1], ir[0]:ir[1]]

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
        # Using helper functions, statement below works for either scalar or tuple index calls:
        return self.ni(it,ir) * constants.mp * self.Abar(it,ir)

    def P(self, it, ir):
        """Pressure P(t,r).

         :param it: the temporal index
         :param ir: the spatial index
         :returns: P [GBar]
         """
        # Following statement should work for either type of specified index:
        return constants.kB * (self.ni(it,ir)*self.Ti(it,ir) + self.ne(it,ir)*self.Te(it,ir)) * 11600*1000*1e-15

    def vol(self, it, ir):
        """Volume of zone index ir at time index it.

         :param it: the temporal index
         :param ir: the spatial index
         :returns: Zone volume [cm^3]
         """
        it, ir = self.__internalIndex__(it, ir)

        if np.isscalar(it) and np.isscalar(ir):
            return self.VolRaw[it][ir]
        return self.VolRaw[it[0]:it[1], ir[0]:ir[1]]

    def rhoR(self, it, ir):
        """Areal density.

         :param it: the temporal index
         :param ir: the spatial index
         :returns: rhoR [g/cm2]
         """
        it, ir = self.__internalIndex__(it, ir)

        if np.isscalar(it) and np.isscalar(ir):
            return self.rho(it,ir)*(self.r_raw[it][ir+1]-self.r_raw[it][ir])

        ret = np.ndarray((it[1]-it[0],ir[1]-ir[0]), dtype=np.float)
        for i in np.arange(it[0], it[1], 1):
            for j in np.arange(ir[0], ir[1], 1):
                ret[i-it[0]][j-ir[0]] = self.rho(i,j)*(self.r_raw[i][j+1]-self.r_raw[i][j])
        return ret

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
        return self.itc

    def it_tstag(self):
        """Time corresponding to stagnation, defined as the minimum shell radius.

        :returns: The time index for tstag
        """
        return self.itstag

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
        return self.iFuel

    def ir_max(self):
        """Maximum radial index (exclusive) in this implosion. Chosen to be exclusive so one can::

            for ir in range(foo.ir_min(), foo.ir_max()):
                bar(ir)

        :returns: The radial index"""
        return len(self.TiRaw[0])

    def t(self, it):
        """Convert indices to real time.

        :param it: The temporal index
        :returns: real time [s]
        """
        assert np.isscalar(it) or (isinstance(it, tuple) and len(it)==2)

        if np.isscalar(it):
            return self.t_raw[it]
        return self.t_raw[it[0]:it[1]]


    def dt(self, it, ir=None):
        """Get the post-processor time step.

        :param it: Time index (may be ignored if implosions have a constant time step)
        :param ir: Optional since dt is the same for all spatial zones. However, this gives the option of getting a
            2-D array as the return value, which is convenient for some calculations.
        :returns: The delta in time between it and it+1 [s]
        """
        # Case where no ir is specified:
        if ir is None:
            assert np.isscalar(it) or (isinstance(it, tuple) and len(it)==2)
            if np.isscalar(it):
                return self.dtRaw[it,0]
            else:
                return self.dtRaw[it[0]:it[1], 0]

        # Both it and ir are specified:
        it, ir = self.__internalIndex__(it, ir)

        # Handle scalar it:
        if np.isscalar(it) and np.isscalar(ir):
            return self.dtRaw[it, ir]

        # Handle array:
        return self.dtRaw[it[0]:it[1], ir[0]:ir[1]]

    def r(self, it, ir):
        """Get physical radius for a zone.

        :param it: The temporal index
        :param ir: The spatial index
        :returns: r [cm]
        """
        it, ir = self.__internalIndex__(it, ir)

        # Handle scalar arguments:
        if np.isscalar(it) and np.isscalar(ir):
            return self.rcm2[it][ir]

        # Handle ranges:
        return self.rcm2[it[0]:it[1], ir[0]:ir[1]]

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
            return self.IonARaw[it,ir]
        return self.IonARaw[it[0]:it[1], ir[0]:ir[1]]

    def IonZ(self, it, ir):
        """Ion atomic numbers for ions in a zone.

        :param it: The temporal index
        :param ir: The spatial index
        :returns: A :py:class:`numpy.ndarray` of ion atomic numbers [e]
        """
        it, ir = self.__internalIndex__(it, ir)

        if np.isscalar(it) and np.isscalar(ir):
            return self.IonZRaw[it,ir]
        return self.IonZRaw[it[0]:it[1], ir[0]:ir[1]]

    def IonF(self, it, ir):
        """Ion fractions in a zone. Each fraction is between 0 and 1 (inclusive). The sum of all fractions is 1.

        :param it: The temporal index
        :param ir: The spatial index
        :returns: A :py:class:`numpy.ndarray` of ion fractions
        """
        it, ir = self.__internalIndex__(it, ir)

        if np.isscalar(it) and np.isscalar(ir):
            return self.IonFRaw[it,ir]
        return self.IonFRaw[it[0]:it[1], ir[0]:ir[1]]

    def Abar(self, it, ir):
        """Average ion mass in a zone.

        :param it: The temporal index
        :param ir: The spatial index
        :returns: Average ion mass [AMU]
        """
        it, ir = self.__internalIndex__(it, ir)

        if np.isscalar(it) and np.isscalar(ir):
            return self.IonAbarRaw[it,ir]
        return self.IonAbarRaw[it[0]:it[1], ir[0]:ir[1]]

    def Zbar(self, it, ir):
        """Average ion Z.

        :param it: The temporal index
        :param ir: The spatial index
        :returns: Average ion atomic number [e]
        """
        it, ir = self.__internalIndex__(it, ir)

        if np.isscalar(it) and np.isscalar(ir):
            return self.IonZbarRaw[it,ir]
        return self.IonZbarRaw[it[0]:it[1], ir[0]:ir[1]]

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
            Azone = self.IonA(it, ir)
            Zzone = self.IonZ(it, ir)
            for i in range(len(Azone)):
                if A == Azone[i] and Z == Zzone[i]:
                    return True
            return False

        # Truth values, for if a given entry corresponds to the ion we want
        shape = (it[1]-it[0], ir[1]-ir[0], len(self.IonARaw[0][0]))
        testA = np.ones(shape, dtype=np.float) * A
        testZ = np.ones(shape, dtype=np.float) * Z
        truthA = np.equal(testA, self.IonARaw)
        truthZ = np.equal(testZ, self.IonZRaw)
        return np.sum(self.IonFRaw*truthA*truthZ, axis=2)

    # ----------------------------------------
    #           Helper functions
    # ----------------------------------------
    def __readHYADES__(self):
        """Read hydro data from the file."""
        file = scipy.io.netcdf.netcdf_file(self.filename, 'r')
        self.runProgress += 0.05

        # raw hydro info
        self.t_raw = np.asarray( file.variables['DumpTimes'][:] )
        self.runProgress += 0.05
        self.r_raw = np.asarray( file.variables['R'][:] )
        self.runProgress += 0.05
        self.rcm = np.asarray( file.variables['Rcm'][:] )
        self.runProgress += 0.05

        # fix weirdness with rcm
        temp = []
        for i in range(len(self.rcm)):
            temp.append( self.rcm[i][1:len(self.rcm[i])-1] )
        self.rcm2 = np.asarray( temp )
        self.runProgress += 0.05

        self.TiRaw = np.asarray( file.variables['Ti'][:] )
        self.runProgress += 0.05
        self.TeRaw = np.asarray( file.variables['Te'][:] )
        self.runProgress += 0.05
        self.niRaw = np.asarray( file.variables['Deni'][:] )
        self.runProgress += 0.05
        self.neRaw = np.asarray( file.variables['Dene'][:] )
        self.runProgress += 0.05
        self.uRaw = np.asarray( file.variables['Ucm'][:] )
        self.runProgress += 0.05
        self.VolRaw = np.asarray( file.variables['Vol'][:] )
        self.runProgress += 0.05

        # raw material info
        self.RegNums = np.asarray( file.variables['RegNums'][:] )
        self.runProgress += 0.05
        self.NumMatsReg = np.asarray( file.variables['NumMatsReg'][:] )
        self.runProgress += 0.05
        self.AtmFrc = np.asarray( file.variables['AtmFrc'][:] )
        self.runProgress += 0.05
        self.AtmNum = np.asarray( file.variables['AtmNum'][:] )
        self.runProgress += 0.05
        self.AtmWgt = np.asarray( file.variables['AtmWgt'][:] )
        self.runProgress += 0.05

        # Work on material definitions
        self.NumRegions = max(self.RegNums) # how many regions we have to deal with

        # Precompute the dt matrix:
        self.dtRaw = np.ndarray(shape=(len(self.t_raw), len(self.r_raw)), dtype=np.float)
        for i in range(len(self.dtRaw)):
            # edge case:
            if i >= self.it_max()-1:
                self.dtRaw[i, self.ir_min():self.ir_max()] = self.t_raw[self.it_max()-1] - self.t_raw[self.it_max()-2]
            else:
                self.dtRaw[i, self.ir_min():self.ir_max()] = self.t_raw[i+1] - self.t_raw[i]

        self.runProgress += 0.05
        file.close()

    def __rRegionInit__(self):
        """Helper function to do initial region boundary position identification."""
        ZoneBoundaries = []
        for ir in range(1,self.NumRegions+1):
            iz = 0
            while iz < len(self.RegNums) and self.RegNums[iz] <= ir:
                iz += 1
            #subtract 2 b/c RegNums has 0s at start and end, 1 more for while loop overshoot
            ZoneBoundaries.append(iz-3)

        self.rRegionInt = []
        for iz in ZoneBoundaries:
            temp = []
            #populate lists from 2D data
            for i in range(len(self.t_raw)):
                temp.append( self.rcm2[i][iz] )
            self.rRegionInt.append( temp )

        # find where the fuel is in the problem:
        iFuel = 0
        for j in range(self.ir_max()): #find out how many regions actually contain stuff
            if fusion.fuel(self.IonA(0,j), self.IonZ(0,j)):
                iFuel = max(iFuel, j)
        self.iFuel = iFuel

    def __rRegion__(self, it, ir):
        """Region boundary position at time index it. Returns r in cm."""
        return self.rRegionInt[ir][it]

    def __RegionIdent__(self, ir):
        """For radial index ir, returns material region number."""
        return self.RegNums[ir]

    def __find_tc__(self):
        """Find the critical timescales in the problem"""
        # Find shock coalescence time by somewhat crude method:
        tc = self.t_raw[0]
        prevTi = self.Ti(0, 0)
        Ti = self.Ti(1, 0)
        i = 1
        #iterate until the ion temp changes by more than 1keV
        while (i+1) < len(self.t_raw):
            i += 1
            prevTi = Ti
            Ti = self.Ti(i, 0)
            if Ti > 0. and (Ti-prevTi) >= 1:
                tc = self.t_raw[i]
                return tc
        self.tc = tc
        self.itc = max(0,i-1)

        # Find stagnation time via minimum shell radius:
        min_R = self.__rRegion__(0, self.iFuel)
        stag_it = 0
        for it in range(len(self.t_raw)):
            R = self.__rRegion__(it, self.iFuel)
            if R < min_R:
                min_R = R
                stag_it = it
        self.tstag = self.t_raw[stag_it]
        self.itstag = stag_it

    def __precompute_ions__(self):
        """Create arrays for ion composition in a more convenient structure."""
        # Need to know the max number of ions per zone:
        maxIons = 0
        for i in range(len(self.AtmWgt)):
            maxIons = max(len(self.AtmWgt[i]), maxIons)

        # Size of created 3-D arrays for IonA, IonZ, IonF:
        shape = (self.it_max()-self.it_min(), self.ir_max()-self.ir_min(), maxIons)
        self.IonARaw = np.zeros(shape, dtype=np.float)
        self.IonZRaw = np.zeros(shape, dtype=np.float)
        self.IonFRaw = np.zeros(shape, dtype=np.float)
        # Create 2-D arrays for averages (Abar and Zbar):
        shape = (self.it_max()-self.it_min(), self.ir_max()-self.ir_min())
        self.IonAbarRaw = np.zeros(shape, dtype=np.float)
        self.IonZbarRaw = np.zeros(shape, dtype=np.float)

        # Now we unfortunately have to loop to populate these:
        for ir in np.arange(self.ir_min(), self.ir_max(), 1):
            for k in range(len(self.AtmWgt[max(self.RegNums[ir]-1, 0)])):
                self.IonARaw[self.it_min():self.it_max(),ir,k] = self.AtmWgt[max(self.RegNums[ir]-1, 0)][k]
                self.IonZRaw[self.it_min():self.it_max(),ir,k] = self.AtmNum[max(self.RegNums[ir]-1, 0)][k]
                self.IonFRaw[self.it_min():self.it_max(),ir,k] = self.AtmFrc[max(self.RegNums[ir]-1, 0)][k]

        self.IonAbarRaw = np.sum(self.IonARaw*self.IonFRaw, axis=2)
        self.IonZbarRaw = np.sum(self.IonZRaw*self.IonFRaw, axis=2)
