# impy - a post-processor for HYADES implosion simulations
# Copyright (c) Massachusetts Institute of Technology / Alex Zylstra
# Distributed under the MIT License

from impy.implosions.Implosion import Implosion
import os
import csv
import sys
import re
import numpy as np
from impy.resources import fusion
from impy.resources import constants


class LILAC(Implosion):
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
    """
    __author__ = 'Alex Zylstra & Hong Sio'
    __date__ = '2015-01-02'
    __version__ = '0.2.1'

    # ----------------------------------------
    #           Generic methods
    # ----------------------------------------
    def __init__(self, type='GUI', args='', wm=None):
        """Construct a new implosion."""
        super(LILAC, self).__init__()
        self.__ready__ = False
        self.__cancelled__ = False

        # Get file to open based on type:
        # we allow all extension, because apparently Darwin doesn't have pull-down menu
        # in open file that allow you to switch between extensions
        if type is 'GUI':
            from tkinter.filedialog import askopenfilename
            FILEOPENOPTIONS = dict(filetypes=[('LILAC','*.*')],
                           multiple=False)
            filename = askopenfilename(**FILEOPENOPTIONS)
        elif type is 'File':
            filename = args
        elif type is 'CLI':
            # TODO: CLI args
            filename = input("LILAC file: ")
        else:
            raise ValueError('Type passed to LILAC constructor is invalid: ' + type)

        # Check that file exists, and that extension is correct:
        if not os.path.isfile(filename):
            raise FileNotFoundError('Could not find LILAC file: ' + filename)
        if not (filename[-6:] == '.lilac' or filename[-5].isdigit()):
            raise FileNotFoundError('Invalid file type in LILAC: ' + filename)

        # Set a few instance variables:
        self.filename = filename
        self.runProgress = 0.
        if self.__cancelled__:
            print(self.ready())
            return

        self.__ready__ = True
        print("Ready to generate implosion?  ",self.__ready__)

    @classmethod
    def getFileTypes(cls):
        """Get a list containing extensions of file types supported by this implosion.
        Must be an array of dicts ready to pass to file dialog, e.g.::

            [('description', '*.ext')]
        """
        return [('reduced LILAC','*.lilac')]

    @classmethod
    def name(cls):
        """Get a string containing a name for this type of implosion."""
        return 'LILAC'

    def info(self):
        """Get a string of information about this specific implosion."""
        #TODO: implement more interesting info
        return 'LILAC file: ' + self.filename

    def ready(self):
        """Returns true if implosion object creation went OK and this object is ready for `generate` to be called."""
        return self.__ready__

    def generate(self):
        """Run the calculation to generate the implosion data."""
        # Top level construction, uses several helper functions:
        print("Beginning generating implosion:  ")
        self.runProgress = 0.
        # for old style, reduced lilac
        if ".lilac" == self.filename[-6:]:
            self.__read_materials_reduced_LILAC__()
            self.__read_reduced_LILAC__()
        # otherwise, treat it as the raw lilac dump file
        else:
            self.__read_LILAC__()

        self.__post_process_raw_data__()

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
            return self.niRaw[it][ir] * self.IonZbarRaw[it][ir]
        return self.niRaw[it[0]:it[1], ir[0]:ir[1]] * self.IonZbarRaw[it[0]:it[1], ir[0]:ir[1]]

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
        return np.sqrt((5 / 3) * self.P(it, ir) * 1e15 / self.rho(it, ir))

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

        ret = self.rho(it, (ir[0],ir[1]-1)) * np.diff(self.r(it,ir), n=1)
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
            return self.r_raw[it][ir]

        # Handle ranges:
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
                    return self.IonF(it,ir)[i]
            return 0.

        # Truth values, for if a given entry corresponds to the ion we want
        shape = (it[1]-it[0], ir[1]-ir[0], len(self.IonARaw[0][0]))
        testA = np.ones(shape, dtype=np.float) * A
        testZ = np.ones(shape, dtype=np.float) * Z
        truthA = np.equal(testA, self.IonARaw)
        truthZ = np.equal(testZ, self.IonZRaw)
        return np.sum(self.IonFRaw*truthA*truthZ, axis=2)



    # ----------------------------------------
    #           Helper functions for parsing reduced lilac file
    # ----------------------------------------
    def __read_materials_reduced_LILAC__(self):
        """Read material info from the top of the file."""
        # initialize a few lists for figuring out the material regions:
        self.NumRegions = 0
        self.Regions = []
        self.rRegionMax = []
        self.rRegionInt = []

        self.region_A = []
        self.region_Z = []
        self.region_F = []
        self.region_Aavg = []
        self.region_Zavg = []

        with open(self.filename, 'r') as file:
            # First line is the name
            name = file.readline()

            # Read in the fuel region info
            dataReader = csv.reader(file , delimiter=' ')
            line = self.cleanRead(dataReader)
            while line[0].isdigit() and not self.__cancelled__:
                #add to regions
                self.Regions.append( [int(line[len(line)-2]), int(line[len(line)-1]) ] ) #region boundaries
                # update materials
                self.MatIdent(line[1], self.Regions[-1])
                self.NumRegions += 1
                line = self.cleanRead(dataReader)

    def __read_reduced_LILAC__(self):
        """Read hydro data from the file."""

        with open(self.filename, 'r') as file:
            #read in the name
            self.sim_name = file.readline()
            #convert from Dy/cm2 to GBar
            PUnitConv = 1e-15
            #convert from eV to keV
            TUnitConv = 1e-3

            dataReader = csv.reader(file, delimiter=',')
            line = self.cleanRead(dataReader)

            #Read in the header until we get to the start of the data
            #switch to comma delimited
            dataReader = csv.reader(file , delimiter=',')
            prevLine = line
            while not 'time' in line[0]:
                prevLine = line
                line = self.cleanRead(dataReader)
            header = prevLine
            #iterate through header, identify columns that we want
            for i in range( len(header) ):
                if "Distance" in header[i]:
                    rIndex = i
                if ("Ion" or "ion") in header[i] and "eV" in header[i]:
                    TiIndex = i
                if ("Ion" or "ion") in header[i] and "Density" in header[i]:
                    niIndex = i
                if ("elec" or "Elec") in header[i] and "eV" in header[i]:
                    TeIndex = i
                if "Velocity" in header[i]:
                    uIndex = i
                if "Pressure" in header[i]:
                    PIndex = i

            self.runProgress += 0.1

            # 2-D arrays for the actual data
            self.r_raw = []
            self.t_raw = []
            self.TiRaw = []
            self.TeRaw = []
            self.niRaw = []
            self.uRaw = []
            self.PRaw = []

            #Read in the actual data
            dataReader = csv.reader(file , delimiter=' ')
            #first time snap was previously read
            line = [ "time=" , line[0].split(' ')[1] ]
            while len(line) > 0 and 'time=' in line[0]:
                #read in one time snapshot:
                t = float( line[1] )
                # for some reason LILAC spits out double sometimes (le sigh)
                if len(self.t_raw) == 0 or t != self.t_raw[-1]:
                    self.t_raw.append( t )
                    self.r_raw.append( [] )
                    self.TiRaw.append( [] )
                    self.TeRaw.append( [] )
                    self.niRaw.append( [] )
                    self.uRaw.append( [] )
                    self.PRaw.append( [] )
                    line = self.cleanRead(dataReader)
                    zone = 0
                    while len(line) > 0 and line[0][0].isdigit():
                        r = float( line[rIndex] ) #radius
                        self.r_raw[-1].append(r)
                        # read in all variables of interest
                        self.TiRaw[-1].append(float( line[TiIndex] )*TUnitConv)
                        self.TeRaw[-1].append(float( line[TeIndex] )*TUnitConv)
                        self.niRaw[-1].append(float( line[niIndex] ))
                        self.uRaw[-1].append(float( line[uIndex] ))
                        self.PRaw[-1].append(float( line[PIndex] )*PUnitConv)
                        line = self.cleanRead(dataReader)
                        zone += 1
                else:  # need to discard this output
                    line = self.cleanRead(dataReader)
                    while len(line) > 0 and line[0][0].isdigit():
                        line = self.cleanRead(dataReader)

            self.runProgress += 0.4

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

    def MatIdent(self, Num, Zones):
        """Read in LILAC material definitions from CSV file, looking for ID # Num."""
        try:
            if getattr(sys, 'frozen', False):
                # The application is frozen
                datadir = os.path.dirname(sys.executable)
            else:
                # The application is not frozen
                # Change this bit to match where you store your data files:
                datadir = os.path.dirname(__file__)
            #dir = os.path.split(__file__)
            fname = os.path.join(datadir, 'LILAC_Materials.csv')
            print(fname)
            MatFile = open(fname, 'r')
        except:
            print('could not open materials file')
            return

        # TODO: data file input needs to be robust
        # TODO: better material error handling
        dataReader = csv.reader(MatFile , delimiter=',')

        A = []
        Z = []
        F = []
        Abar = 0
        Zbar = 0
        header = []
        for row in dataReader:
            # first row of the file is the header:
            if len(header) == 0:
                header = row
            # found row corresponding to this material:
            elif int(row[0]) == int(Num):
                Abar = float(row[1])
                Zbar = float(row[2])
                # number of ions:
                n = (len(row) - 3)/3
                A = np.ndarray(n, dtype=np.float)
                Z = np.ndarray(n, dtype=np.float)
                F = np.ndarray(n, dtype=np.float)

                # add each individual ion to the arrays:
                for i in range(int(n)):
                    A[i] = float(row[3*(1+i)])
                    Z[i] = float(row[3*(1+i)+1])
                    F[i] = float(row[3*(1+i)+2])

                break

        #sanity check
        if Abar == 0 or Zbar == 0:
            print("ERROR: material not found!")
            # Prompt the user:
            from impy.gui.LILAC_Material import LILAC_Material
            dialog = LILAC_Material(None, Num, Zones)

            # If the user cancelled, abort:
            if dialog.cancelled:
                self.__cancelled__ = True
                return

            # Get the results:
            A = np.asarray(dialog.result[1])
            Z = np.asarray(dialog.result[2])
            F = np.asarray(dialog.result[3])
            Abar = np.dot(A,F)
            Zbar = np.dot(Z,F)

        self.region_A.append(A)
        self.region_Z.append(Z)
        self.region_F.append(F)
        self.region_Aavg.append(Abar)
        self.region_Zavg.append(Zbar)


    # ----------------------------------------
    #           Helper functions for full lilac dump
    # ----------------------------------------
    def find_next_num(self, L, ii=0):
        """"Helper function to find next number in a list given starting index.

        :param L: A list of str
        :param ii: starting index to search
        :return: the first str that can be interpreted as a number, cast as a number
        """
        for item in L[ii:]:
            try:
                return float(item)
            except:
                continue
        return None

    def float_(self, str):
        """exception safe version of float().

        :param str: a string, with the intention of being cast a number
        :return:  the str, cast as a float.  If not possible, return None
        """
        try:
            ans = float(str)
        except:
            ans = None
        return ans

    def split_str_into_words_min_two_spaces(self, str):
        """Split string into list, separated by more than one space.

        :param str:  a string, to be split into a list of words
        :return:  a list of strings, split if separated by at least two spaces
        """

        # This is a regular expression that split by more than or equal two spaces
        line = re.split('\s{2}', str)
        # clean up the white space
        while '' in line:
            line.remove('')
        headings = []
        for item in line[:]:
            if item[0] == " ":
                item = item[1:]
            else:
                item = item
            headings.append(item)
        return headings

    # ----------------------------------------
    #          Dictionaries for the full lilac dump
    # ----------------------------------------

    lilacKeys = {
        "L_TOKEN_ENDL":'\n',
        "L_TOKEN_WHITESPACE":" ",                    #white space
        "L_TOKEN_NUMERICAL":"0123456789eE+-. ",      #valid numbers

        # Keys used for the master parse function
        "KEY_LILAC":"LILAC",
        "ZONE_INDEX_SHIFT":-2,
        "KEY_LASER_HEADER":"0-----",
        "KEY_PELLET_HEADER":"*** pellet structure ***",
        "KEY_PELLET_FOOTER":" ****",
        "KEY_TIME_STEP_HEADER":"0********** conditions after time step:",
        "KEY_TIME_STEP_FOOTER":" ****",
        "KEY_TIME_STEP_FOOTER_COUNT":3,
        "KEY_RUN_SUMMARY_BREAK":"rs*****",
        "KEY_RUN_SUMMARY_HEADER":"********** run summary",
        "KEY_RUN_SUMMARY_FOOTER":"+ exec",
        "KEY_ZONING_HEADER":"***zoning information***",
        "KEY_ZONING_FOOTER":"0pellet total:",

        # Keys used for the pellet parsing function
        "KEY_PELLET_LAYER_DELIMITER":"-------------------------------------",
        "KEY_PELLET_INDEX_ELEM":1,
        "KEY_PELLET_INDEX_FRAC":2,
        "KEY_PELLET_INDEX_A":5,
        "KEY_PELLET_INDEX_ELEM_NEXT":0,      # pellet structure past first line of each layer
        "KEY_PELLET_INDEX_FRAC_NEXT":1,      # pellet structure past first line of each layer
        "KEY_PELLET_INDEX_NUM_CELLS":8,
        "KEY_PELLET_NUM_LAYER_INFO_NEXT":3,  # number of elements per row past first line of each layer

        # Tokens used for zoning info parsing
        "KEY_ZONE_INFO_INDEX":"j",
        "KEY_ZONE_INFO_MAT":"mat",
        "KEY_ZONE_INFO_R":"R center",
        "KEY_ZONE_INFO_DR":"dr",
        "KEY_ZONE_INFO_VOL":"volume",
        "KEY_ZONE_INFO_MASS_C":"mass(c)",
        "KEY_ZONE_INFO_MASS_B":"mass(b)",
        "KEY_ZONE_INFO_IONS":"# of ions",

        #  Tokens used for zone data parsing
        "KEY_ZONE_DATA_INDEX":"j",
        "KEY_ZONE_DATA_MAT":"m",
        "KEY_ZONE_DATA_R":"R",
        "KEY_ZONE_DATA_VEL":"vel",
        "KEY_ZONE_DATA_TI":"T (ion)",
        "KEY_ZONE_DATA_TE":"T (el.)",
        "KEY_ZONE_DATA_DENS":"density",
        "KEY_ZONE_DATA_EL_DENS":"el. dens",
        "KEY_ZONE_DATA_Z":"Z",
        "KEY_ZONE_DATA_PRESS":"press",
        "KEY_ZONE_DATA_NEUTRON":"neutron",
        "KEY_ZONE_DATA_E_DEP":"E. dep.",
        "KEY_ZONE_DATA_E_RAD":"E. rad.",
        "KEY_ZONE_DATA_VISC_P":"visc. p.",
        "KEY_ZONE_DATA_COMP":"comp.",

        # Tokens used by time-step parsing routines
        "KEY_ZONE_VEL":"vel",
        "KEY_ZONE_DENSITY":"density",
        "KEY_ZONE_END":"****************************",
        "KEY_SECTION_LASER_START":"********** laser conditions:",
        "KEY_SECTION_PLASMA_START":"********** plasma conditions:",
        "KEY_SECTION_RADIATION_START":"********** radiation:",
        "KEY_SECTION_THERMONUCLEAR_START":"********** thermonuclear reactions:",
        "KEY_SECTION_DYNAMIC_START":"********** dynamic summary:",
        "KEY_SECTION_SMOOTHING_START":"***** smoothing information *****",
        "KEY_SECTION_LAYERENERGY_START":"***** cumulative layer-averaged energies *****",
        "KEY_SECTION_TIMESTEP_START":"***** Time step control *****",
        "KEY_SECTION_LANGDON_START":"***** langdon factor *****",
        "KEY_TIME":"and time:",
        "KEY_LASER_POWER":"power:",
        "KEY_LASER_INTENSITY":"exp. intensity:",
        "KEY_LASER_TOTAL_ENERGY":"energy supplied:",
        "KEY_LASER_ABSORPTION_TOTAL":"absorption fraction  :",
        "KEY_LASER_ABSORPTION_STEP":"step abs. fraction   :",
        "KEY_LASER_CRITICAL_INDEX":"critical surface: j=",
        "KEY_LASER_QUARTER_CRITICAL_INDEX":"quarter critical: j=",
        "KEY_LASER_TENTH_CRITICAL_INDEX":"tenth critical  : j=",
        "KEY_LASER_CRITICAL_POS":"radius:",
        "KEY_LASER_QUARTER_CRITICAL_POS":"radius:",
        "KEY_LASER_TENTH_CRITICAL_POS":"radius:",
        "KEY_TOTAL_ENERGY":"total energy:",
        "KEY_THERMAL_ENERGY":"thermal mode:",
        "KEY_KINETIC_ENERGY":"kinetic mode:",
        "KEY_DYNAMICS_HEADER_LINE":"kin. energy",
        "KEY_DYNAMICS_UNIT_LINE":"Joules",
        "KEY_DYNAMICS_IMPLODED_LINE":"imploded:",
        "KEY_DYNAMICS_ABLATED_LINE":"ablated :",
        "KEY_DYNAMICS_ABLATION_RATE":"mass ablation rate:",
        "KEY_DYNAMICS_MASS":"mass",
        "KEY_DYNAMICS_ENERGY":"kin. energy",

        #Useful parameters
        "ELEC_FLUX_LIMITER":"electron flux limit coef.",
        "ARTIFICIAL_VISCOSITY_COEF":"art. visc. coef.",
        "SIMULATION_TIME_END":"tcut"
    }

    elemDict = {        # in the format of [A, Z]
        "H":[1,1],
        "D":[2,1],
        "T":[3,1],
        "He3":[3,2],
        "He":[4,2],
        "Li":[6.94, 3],
        "Be":[9.01,4],
        "C":[12.011,6],
        "CARBON":[12.011,6],
        "N":[14.007,7],
        "O":[15.999,8],
        "Ne":[20.180,10],
        "Na":[22.99,11],
        "Al":[26.982, 13],
        "Mg":[24.305,12],
        "Si":[28.005,14],
        "S":[32.06,16],
        "Ar":[39.948,18],
        "Ti":[47.867,22],
        "Fe":[55.845,26],
        "Ni":[58.693,28],
        "Cu":[63.546,29]
    }

    # ----------------------------------------
    #           Primary functions for full lilac dump
    # ----------------------------------------
    def __read_LILAC__(self):

        # 2-D arrays for the actual data
        self.r_raw = []
        self.t_raw = []
        self.TiRaw = []
        self.TeRaw = []
        self.niRaw = []
        self.uRaw = []
        self.PRaw = []
        self.neRaw = []
        self.ZRaw = []
        self.rhoRaw = []

        # parsing flags - not all are used
        self.is_valid_lilac_file   = False
        self.laser_info_processed  = False
        self.pellet_info_processed = False
        self.run_summary_processed = False
        self.time_steps_processed  = False
        self.zoning_info_processed = False

        #convert from MBar to GBar
        self.PUnitConv = 1e-3
        #convert from eV to keV
        self.TUnitConv = 1e-3
        #convert microns to cm
        self.RUnitConv = 1e-4

        with open(self.filename, 'r') as file:
            print(self.filename)

            dataReader = csv.reader(file)

            while True:
                try:
                    line = next(dataReader)
                    line = line[0]
                except StopIteration:           #finish parsing, when end of file is reached
                    print("EOF reached")
                    break
                except:                         #if something wrong with a line, go on to the next one
                    continue

                # If we arrived at the beginning of a time block, pass it to a parsing function
                if self.lilacKeys["KEY_TIME_STEP_HEADER"] in line:
                    self.lilac_full_parse_time_step(dataReader)

                # if we run into the pellet information, record it
                elif self.lilacKeys["KEY_PELLET_HEADER"] in line:
                    self.lilac_full_parse_pellet_info(dataReader)

                # if it runs into the time simulation ends, record it
                elif self.lilacKeys["SIMULATION_TIME_END"] in line:
                    temp = line[line.find(self.lilacKeys["SIMULATION_TIME_END"]):]
                    temp = temp.split()
                    self.t_end_ns = self.find_next_num(temp)
                    print("End of simulation time is ",self.t_end_ns," ns")

                # if it runs into the artificial viscosity multiplier, record it
                elif self.lilacKeys["ARTIFICIAL_VISCOSITY_COEF"] in line:
                    temp = line[:line.find(self.lilacKeys["ARTIFICIAL_VISCOSITY_COEF"])]
                    temp = temp.split()
                    self.artificial_visc_coeff =self.find_next_num(temp)
                    print("Artificial viscosity coefficient is ",self.artificial_visc_coeff)

                # if it runs into the electron flux limiter, record it
                elif self.lilacKeys["ELEC_FLUX_LIMITER"] in line:
                    temp = line[:line.find(self.lilacKeys["ELEC_FLUX_LIMITER"])]
                    temp = temp.split()
                    self.electron_flux_limiter = self.find_next_num(temp)
                    print("electron flux limiter is ",self.electron_flux_limiter)

                # At some point before time-step start the word "LILAC" has to appear
                # for this to be a valid file
                elif not self.is_valid_lilac_file and self.lilacKeys["KEY_LILAC"] in line:
                    is_valid_lilac_file = True

            time_steps_processed = True
            self.runProgress += 0.4

    def lilac_full_parse_time_step(self,dataReader):
        """Parse one time block in the full lilac dump file."""
        # if this function is called, assume it is at beginning of time block

        # this will store an ordered list of the quantity name
        headings = None
        # this will provide the next list of numbers to be inserted
        line = None

        def insert_data( dict_key):
            try:
                temp_index = headings.index(self.lilacKeys[dict_key])
                if "KEY_ZONE_DATA_R" == dict_key:
                    self.r_raw[-1].append(line[temp_index] * self.RUnitConv)
                elif "KEY_ZONE_DATA_TI" == dict_key:
                    self.TiRaw[-1].append(line[temp_index] * self.TUnitConv)
                elif "KEY_ZONE_DATA_TE" == dict_key:
                    self.TeRaw[-1].append(line[temp_index] * self.TUnitConv)
                elif "KEY_ZONE_DATA_PRESS" == dict_key:
                    self.PRaw[-1].append(line[temp_index] * self.PUnitConv)
                elif "KEY_ZONE_DATA_VEL" == dict_key:
                    self.uRaw[-1].append(line[temp_index])
                elif "KEY_ZONE_DATA_EL_DENS" == dict_key:
                    self.neRaw[-1].append(line[temp_index])
                elif "KEY_ZONE_DATA_Z" == dict_key:
                    self.ZRaw[-1].append(line[temp_index])
                elif "KEY_ZONE_DATA_DENS" == dict_key:
                    self.rhoRaw[-1].append(line[temp_index])
            except:
                print("some basic hydro information is missing")
                pass

        # add a new row to each 2D array for this new time step
        self.r_raw.append([])
        self.TiRaw.append([])
        self.TeRaw.append([])
        self.uRaw.append([])
        self.PRaw.append([])
        self.neRaw.append([])
        self.ZRaw.append([])
        self.niRaw.append([])
        self.rhoRaw.append([])

        #first, fetch the next line,  and time should be the first
        while True:
            try:
                line = next(dataReader)
                line = line[0]
            except:
                continue

            # This handle the case of a blank line sometime
            if line.split() == []:
                continue

            # get the time.  we want to store the time in ns
            elif self.lilacKeys["KEY_TIME"] in line:
                line = line.split()
                time_ns = self.find_next_num(line)
                if "fs" in line:
                    time_ns *= 1e-6
                elif "ps" in line:
                    time_ns *= 1e-3
                elif "ns" in line:
                    pass
                else:
                    print("no proper time unit found in lilac_full_parse_time_step")
                print( "time is ", time_ns, " ns")
                self.t_raw.append( time_ns )

            # if first word is index, we arrived at the headings
            elif self.lilacKeys["KEY_ZONE_DATA_INDEX"] == line.split()[0]:
                # now we need to find the order of the headings
                headings = self.split_str_into_words_min_two_spaces(line)
                if time_ns == 0.0:
                    print("\nheading for zone info")
                    print(headings)

            # hey, are you a zone information row?  yes you are.  yes you are!
            elif line.split()[0].isdigit():
                line = line.split()
                # if it starts with a digit but has only three lne, end of time step info
                if len(line) == self.lilacKeys["KEY_TIME_STEP_FOOTER_COUNT"]:
                    break
                line = [self.float_(item) for item in line]
                # attempt to insert data into respective array
                insert_data("KEY_ZONE_DATA_R")
                insert_data("KEY_ZONE_DATA_TI")
                insert_data("KEY_ZONE_DATA_TE")
                insert_data("KEY_ZONE_DATA_PRESS")
                insert_data("KEY_ZONE_DATA_VEL")
                insert_data("KEY_ZONE_DATA_EL_DENS")
                insert_data("KEY_ZONE_DATA_Z")
                insert_data("KEY_ZONE_DATA_DENS")

                # ion density is density divided by (avg. atomic number) * proton mass
                cell = line[0]
                for L in self.Regions:
                    if cell >= L[0] and cell <= L[1]:
                        region = L
                region_index = self.Regions.index(region)
                Aavg = self.region_Aavg[region_index]
                mass_per_ion = Aavg * constants.mp
                self.niRaw[-1].append( self.rhoRaw[-1][-1] / mass_per_ion )

    def lilac_full_parse_pellet_info(self,dataReader):
        """Parse pellet information in full lilac dump file."""
        # At this point, assume that the pellet header has been found

        # next line is not only useless, but can be accidentally confused with being the footer
        line = next(dataReader)

        # We will try to populate these variables
        self.NumRegions = 0
        self.Regions = []
        self.region_A = []
        self.region_Z = []
        self.region_F = []
        self.region_Aavg = []
        self.region_Zavg = []

        print("\nPellet info found!")

        #first, fetch the next line,  and time should be the first
        while True:
            try:
                line = next(dataReader)
                line = line[0]
            except:
                continue

            # This handle the case of a blank line sometime
            if line.split() == []:
                continue

            # all layers info start with an integer or some sort
            elif line.split()[0].isdigit():
                layer = self.split_str_into_words_min_two_spaces(line)
                # if this is the beginning of a new layer
                if len(layer) >self.lilacKeys["KEY_PELLET_NUM_LAYER_INFO_NEXT"]:
                    self.region_A.append([])
                    self.region_Z.append([])
                    self.region_F.append([])
                    self.NumRegions += 1
                    numCell = float(layer[self.lilacKeys["KEY_PELLET_INDEX_NUM_CELLS"]])
                    if len( self.Regions ) == 0:
                        self.Regions.append([2, numCell + 1])
                    else:
                       self.Regions.append([self.Regions[-1][1]+1,self.Regions[-1][1]+ numCell] )
                    # and then, fetch the element name
                    name_elem = layer[self.lilacKeys["KEY_PELLET_INDEX_ELEM"]]
                    # and then, fetch the atomic fraction
                    frac_elem = float(layer[self.lilacKeys["KEY_PELLET_INDEX_FRAC"]])
                # atomic fraction and name_elem has a different index if it is not on the first line
                else:
                    name_elem = layer[self.lilacKeys["KEY_PELLET_INDEX_ELEM_NEXT"]]
                    frac_elem = float(layer[self.lilacKeys["KEY_PELLET_INDEX_FRAC_NEXT"]])
                # the element name actually has a number before it we need to get ride of
                name_elem = name_elem.split()[1]
                try:
                    A_temp = self.elemDict[name_elem][0]
                    Z_temp = self.elemDict[name_elem][1]
                except KeyError:
                    print("element not found in dictionary:  ",name_elem)
                # attach these information to class variables
                self.region_A[-1].append(A_temp)
                self.region_Z[-1].append(Z_temp)
                self.region_F[-1].append(frac_elem)

            # if we finished processing one layer, wrap it up
            elif self.lilacKeys["KEY_PELLET_LAYER_DELIMITER"] in line and self.NumRegions > 0:
                self.region_A[-1] = np.asarray(self.region_A[-1])
                self.region_Z[-1] = np.asarray(self.region_Z[-1])
                self.region_F[-1] = np.asarray(self.region_F[-1])
                Abar = np.dot(self.region_A[-1],self.region_F[-1] )
                Zbar = np.dot(self.region_Z[-1],self.region_F[-1] )
                self.region_Aavg.append(Abar)
                self.region_Zavg.append(Zbar)

            elif self.lilacKeys["KEY_PELLET_FOOTER"] in line:
                print(self.region_A)
                print(self.region_Z)
                print(self.region_F)
                print(self.region_Aavg)
                break

    # ----------------------------------------
    #           Helper functions for parsing reduced lilac file
    # ----------------------------------------

    def __post_process_raw_data__(self):
        """Post process the raw data after __read_reduced_LILAC__ and __read_LILAC__."""

        # Convert data to ndarray
        self.r_raw = np.asarray(self.r_raw)
        self.t_raw = np.asarray(self.t_raw)
        self.TiRaw = np.asarray(self.TiRaw)
        self.TeRaw = np.asarray(self.TeRaw)
        self.niRaw = np.asarray(self.niRaw)
        self.uRaw = np.asarray(self.uRaw)
        self.PRaw = np.asarray(self.PRaw)

        # need to convert t to s:
        self.t_raw = self.t_raw * 1e-9

        # Precompute volume of each zone
        self.VolRaw = np.zeros_like(self.TiRaw)
        # First zone's volume is spherical:
        self.VolRaw[:,0] = (4*np.pi/3) * np.power(self.r_raw[:,0], 3)
        # subsequent zones are a difference in spherical volume
        self.VolRaw[:,1:] = (4*np.pi/3) * np.diff(np.power(self.r_raw, 3))

        self.runProgress += 0.25

        # Precompute time difference
        self.dtRaw = np.zeros_like(self.r_raw)
        # Most of this can be computed directly with diff:
        self.dtRaw[0:-1,:] = np.outer(np.diff(self.t_raw), np.ones(self.r_raw.shape[1]))
        # last element set manually
        self.dtRaw[-1,:] = self.t_raw[-1] - self.t_raw[-2]

        # Populate arrays for the material info in each zone/time
        # First, need to know the max # of ions in any zone
        maxIons = 0
        for i in range(len(self.region_A)):
            maxIons = max(len(self.region_A[i]), maxIons)
        shape = (len(self.t_raw), len(self.r_raw[0]), maxIons)  # shape for the material arrays
        self.IonARaw = np.zeros(shape, dtype=np.float)
        self.IonZRaw = np.zeros(shape, dtype=np.float)
        self.IonFRaw = np.zeros(shape, dtype=np.float)
        # Create 2-D arrays for averages (Abar and Zbar):
        shape = (len(self.t_raw), len(self.r_raw[0]))
        self.IonAbarRaw = np.zeros(shape, dtype=np.float)
        self.IonZbarRaw = np.zeros(shape, dtype=np.float)

        # populate the arrays just created
        # Do this by looping over regions. Each region defines material for set of zones for all time
        for i in range(len(self.Regions)):
            ir0 = self.Regions[i][0]
            ir1 = self.Regions[i][1]+1  # specified as inclusive limits by LILAC
            self.IonARaw[:,ir0:ir1] = self.region_A[i]
            self.IonZRaw[:,ir0:ir1] = self.region_Z[i]
            self.IonFRaw[:,ir0:ir1] = self.region_F[i]
            self.IonAbarRaw[:,ir0:ir1] = self.region_Aavg[i]
            self.IonZbarRaw[:,ir0:ir1] = self.region_Zavg[i]

        self.runProgress += 0.25

        self.__find_scales__()

    def __find_scales__(self):
        """Find radial and time scales in the problem."""
        # find where the fuel is in the problem:
        iFuel = 0
        for j in range(self.ir_max()): #find out how many regions actually contain stuff
            if fusion.fuel(self.IonA(0,j), self.IonZ(0,j)):
                iFuel = max(iFuel, j)
        self.iFuel = iFuel

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
        self.tc = tc
        self.itc = max(0,i-1)

        # Find stagnation time via minimum shell radius:
        min_R = self.r_raw[0, self.iFuel]
        stag_it = 0
        for it in range(self.it_max()):
            R = self.r_raw[it, self.iFuel]
            if R < min_R:
                min_R = R
                stag_it = it
        self.tstag = self.t_raw[stag_it]
        self.itstag = stag_it






