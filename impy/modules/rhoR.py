
from impy.modules.Module import Module
from impy.implosions.Implosion import Implosion
from impy.resources.fusion import *

import tkinter as tk
import tkinter.ttk as ttk
import matplotlib, matplotlib.pyplot
import pickle
import platform


class RhoR(Module, tk.Toplevel):
    """RhoR plotting module.

    :param type: The type of constructor to be used. Available options are:

    'GUI': The class should interact directly with the user to get required options (default)

    'CLI': Interact with user via CLI, or take info from args.

    :param args: Additional information, which depends on the type of constructor:

    type='GUI': unused

    type='CLI': A full list of arguments passed to the executable, to be interpreted as the module pleases.

    :param wm: A window manager to use for GUI windows during construction

    :author: Alex Zylstra
    :date: 2014-02-10
    """
    __author__ = 'Alex Zylstra'
    __date__ = '2014-02-12'
    __version__ = '0.1.0'

    # ----------------------------------------
    #           Generic methods
    # ----------------------------------------
    def __init__(self, type='GUI', args='', wm=None):
        """Construct a new instance of this module."""
        Module.__init__(self)

        # set a few instance variables:
        self.runProgress = 0

        # Keep track of all plots created:
        self.plots = dict()

        # results:
        self.t = []
        self.rhoR = []

    @classmethod
    def name(cls):
        """Get a string containing a name for this type of module."""
        return 'rhoR'

    @classmethod
    def info(cls):
        """Get a brief description of this specific module."""
        return 'rhoR plots'

    @classmethod
    def detailedInfo(cls):
        """Get any detailed description and explanation of this module."""
        return 'tbd'


    # ----------------------------------------
    #       Execution and GUI control
    # ----------------------------------------
    def run(self, imp):
        """Run the calculation.

        :param imp: An `Implosion` object.
        """
        assert isinstance(imp, Implosion)
        self.runProgress = 0.

        # limits for calculations (all space/time):
        it = (imp.it_min(), imp.it_max())
        ir = (imp.ir_min(), imp.ir_max())
        ir_fuel = (imp.ir_min(), imp.ir_fuel())
        ir_shell = (imp.ir_fuel(), imp.ir_max())

        # Get and save the info for plots:
        self.time = imp.t(it)
        self.rhoR = np.zeros((it[1]-it[0],3), dtype=np.float64)
        self.rhoR[:,0] = np.sum(imp.rhoR(it, ir_fuel), 1)
        self.runProgress = 0.33
        self.rhoR[:,1] = np.sum(imp.rhoR(it, ir_shell), 1)
        self.runProgress = 0.67
        self.rhoR[:,2] = self.rhoR[:,0] + self.rhoR[:,1]
        self.runProgress = 1.

    def display(self, type='GUI', refresh=False, wm=None):
        """Display the results.

        :param type: The type of interface to be used. Available options are 'GUI' and 'CLI'. Default is GUI.
        :param refresh: (optional) If set to True, only existing windows should be updated.
        :param wm: (optional) Window manager to use for displaying windows
        """

        if type is 'GUI' and not refresh:
            # Configure main window GUI:
            tk.Toplevel.__init__(self, None)
            self.grid()
            # stretch the column to fill all space:
            tk.Grid.columnconfigure(self, 0, weight=1)
            tk.Grid.columnconfigure(self, 1, weight=1)

            # Set window background
            if platform.system() == 'Darwin':
                self.configure(background='#E8E9E8')
            else:
                self.configure(background='#F1F1F1')

            self.__createWidgets__()
            self.title('rhoR')
            self.update_idletasks()
            if wm is not None:
                wm.addWindow(self)

            self.wm = wm
            self.__createWidgets__()

        elif type is 'CLI':
            pass

        elif refresh:
            self.__plot__(refresh=True)

    def progress(self):
        """Get the calculation's progress estimate.

        :returns: Scalar number between 0 and 1.
        """
        return self.runProgress

    def copy(self, other):
        """Copy the results from `other` to this module.

        :param other: Another instance of this class
        """
        assert isinstance(other, RhoR)
        try:
            self.time = np.copy(other.time)
            self.rhoR = np.copy(other.rhoR)
        except:
            pass

    def getPlots(self):
        """Return a list of all :py:class:`matplot.pyplot.Figure` instances created by this module."""
        return self.plots.values()

    # ----------------------------------------
    #           Results handling
    # ----------------------------------------
    def save(self, filename, type='CSV'):
        """Save the results to a file.

        :param filename: The file to save to
        :param type: The type of save file to create:

        'CSV': A CSV containing top-level results

        'pickle': Pickle the calculation results
        """
        if type == 'CSV':
            with open(filename, 'wb') as file:
                np.savetxt(file, self.time, delimiter=',', header='Time:\n')
            with open(filename, 'ab') as file:
                np.savetxt(file, self.rhoR, delimiter=',', header='rhoR:\n')

        elif type == 'pickle':
            with open(filename, 'wb') as file:
                pickle.dump(self.time, file)
                pickle.dump(self.rhoR, file)

    def savePlots(self, dir, prefix, type):
        """Save all generated plots.

        :param dir: The directory to save the plots to
        :param prefix: A prefix to append to any files generated by this module
        :param type: The type of plot to generate (e.g. 'eps')

        As a result, files generated in `dir` should have names like::

            prefix + my_name + '.' + type

        """
        for key in self.plots.keys():
            fname = os.path.join(dir, prefix+'_'+key+'.'+type)
            self.plots[key].savefig(fname, bbox_inches='tight')

    # ----------------------------------------
    #      Stuff needed for pickling
    # ----------------------------------------
    def __getstate__(self):
        """Get the current state of this object as a `dict`"""
        state = self.__dict__.copy()
        badKeys = ['wm','master','tk','_w','widgetName','plots','_tclCommands','_name','children','scalars','fuelVar','shellVar','totalVar']

        for key in badKeys:
            if key in state.keys():
                del state[key]
        return state

    # ----------------------------------------
    #           GUI stuff
    # ----------------------------------------
    def __createWidgets__(self):
        """Helper function which creates GUI elements."""
        plotLabel = ttk.Label(self, text='Plot Options')
        plotLabel.grid(row=2, column=0, columnspan=2, sticky='ns')

        label1 = ttk.Label(self, text='Fuel')
        label1.grid(row=3, column=0)
        self.fuelVar = tk.BooleanVar(value=True)
        fuelCheck = ttk.Checkbutton(self, variable=self.fuelVar)
        fuelCheck.grid(row=3, column=1)

        label2 = ttk.Label(self, text='Shell')
        label2.grid(row=4, column=0)
        self.shellVar = tk.BooleanVar(value=True)
        shellCheck = ttk.Checkbutton(self, variable=self.shellVar)
        shellCheck.grid(row=4, column=1)

        label3 = ttk.Label(self, text='Total')
        label3.grid(row=5, column=0)
        self.totalVar = tk.BooleanVar(value=True)
        totalCheck = ttk.Checkbutton(self, variable=self.totalVar)
        totalCheck.grid(row=5, column=1)

        burnRateButton = ttk.Button(self, text='Plot', command=self.__plot__)
        burnRateButton.grid(row=9, column=0, columnspan=2)

    def __plot__(self, refresh=False, *args):
        """Helper function to generate plot of burn rate"""
        # Check for a closed window:
        if 'rhoR' in self.plots.keys() and not matplotlib.pyplot.fignum_exists(self.plots['rhoR'].number):
            del self.plots['rhoR']
            refresh = False
        # Update the existing plot, if it exists
        refresh = refresh or 'rhoR' in self.plots.keys()
        if refresh:
            if 'rhoR' in self.plots.keys():
                fig = self.plots['rhoR']
                fig = matplotlib.pyplot.figure(fig.number)
                fig.clear()
            else:
                return
        # Make a new window:
        else:
            fig = matplotlib.pyplot.figure(figsize=(4,3))
            fig.canvas.set_window_title('rhoR')
        ax = fig.add_subplot(111)

        # Plot
        if self.fuelVar.get():
            ax.plot(1e9*self.time, 1e3*self.rhoR[:,0], 'r-', label='Fuel')
        if self.shellVar.get():
            ax.plot(1e9*self.time, 1e3*self.rhoR[:,1], 'b-', label='Shell')
        if self.totalVar.get():
            ax.plot(1e9*self.time, 1e3*self.rhoR[:,2], 'k-', label='Total')

        ax.set_xlabel('Time (ns)', fontsize=12)
        ax.set_ylabel('rhoR (mg/cm2)', fontsize=12)
        ax.legend(loc=2)

        matplotlib.pyplot.tight_layout()

        if not refresh:
            if self.wm is not None:
                self.wm.addWindow(matplotlib.pyplot.get_current_fig_manager().window)
            fig.show()
            fig.canvas.draw()
        self.plots['rhoR'] = fig
