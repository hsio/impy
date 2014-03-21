
from impy.modules.Module import Module
from impy.implosions.Implosion import Implosion
from impy.resources.fusion import *

import tkinter as tk
import tkinter.ttk as ttk
import matplotlib, matplotlib.pyplot
import pickle


class Lagrange(Module, tk.Toplevel):
    """Langrangian (mass element) plotting module.

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
    __date__ = '2014-02-10'
    __version__ = '1.0.0'

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
        self.Lagrange = []
        self.fuelShell = []

    @classmethod
    def name(cls):
        """Get a string containing a name for this type of module."""
        return 'Lagrange'

    @classmethod
    def info(cls):
        """Get a brief description of this specific module."""
        return 'Lagrangian mass element trajectories'

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
        # Get and save the info for plots:
        self.time = imp.t(it)
        self.fuelShell = imp.r(it, (imp.ir_fuel()))
        self.Lagrange = imp.r(it, ir)

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
            self.configure(background='#eeeeee')
            self.__createWidgets__()
            self.title('Lagrange')
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
        assert isinstance(other, Lagrange)
        try:
            self.time = np.copy(other.time)
            self.fuelShell = np.copy(other.fuelShell)
            self.Lagrange = np.copy(other.Lagrange)
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
                np.savetxt(file, self.fuelShell, delimiter=',', header='Fuel-shell interface radius:\n')
            with open(filename, 'ab') as file:
                np.savetxt(file, self.Lagrange, delimiter=',', header='radius points for each time step for each zone:\n')

        elif type == 'pickle':
            with open(filename, 'wb') as file:
                pickle.dump(self.time, file)
                pickle.dump(self.fuelShell, file)
                pickle.dump(self.Lagrange, file)

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
        badKeys = ['wm','master','tk','_w','widgetName','plots','_tclCommands','_name','children','scalars',
                   'incrementVar', 'minRVar', 'maxRVar', 'minTVar', 'maxTVar', 'fuelShellVar']

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

        label1 = ttk.Label(self, text='Plot every')
        label1.grid(row=3, column=0)
        self.incrementVar = tk.StringVar(value=10)
        incrementEntry = ttk.Entry(self, textvariable=self.incrementVar, width=8)
        incrementEntry.grid(row=3, column=1)

        label2 = ttk.Label(self, text='Min R (um)')
        label2.grid(row=4, column=0)
        self.minRVar = tk.StringVar(value=0)
        minREntry = ttk.Entry(self, textvariable=self.minRVar, width=8)
        minREntry.grid(row=4, column=1)

        label3 = ttk.Label(self, text='Max R (um)')
        label3.grid(row=5, column=0)
        self.maxRVar = tk.StringVar(value=int(1e4*1.5*np.max(self.Lagrange[0])))
        maxREntry = ttk.Entry(self, textvariable=self.maxRVar, width=8)
        maxREntry.grid(row=5, column=1)

        label4 = ttk.Label(self, text='Min t (ns)')
        label4.grid(row=6, column=0)
        self.minTVar = tk.StringVar(value=0)
        minTEntry = ttk.Entry(self, textvariable=self.minTVar, width=8)
        minTEntry.grid(row=6, column=1)

        label5 = ttk.Label(self, text='Max t (ns)')
        label5.grid(row=7, column=0)
        self.maxTVar = tk.StringVar(value=int(1e9*self.time[-1]))
        maxTEntry = ttk.Entry(self, textvariable=self.maxTVar, width=8)
        maxTEntry.grid(row=7, column=1)

        label6 = ttk.Label(self, text='Show Fuel-Shell')
        label6.grid(row=8, column=0)
        self.fuelShellVar = tk.BooleanVar(value=True)
        fuelShellCheck = ttk.Checkbutton(self, variable=self.fuelShellVar)
        fuelShellCheck.grid(row=8, column=1)

        burnRateButton = ttk.Button(self, text='Plot', command=self.__plot__)
        burnRateButton.grid(row=9, column=0, columnspan=2)

    def __plot__(self, refresh=False, *args):
        """Helper function to generate plot of burn rate"""
        # Check for a closed window:
        if 'Lagrange' in self.plots.keys() and not matplotlib.pyplot.fignum_exists(self.plots['Lagrange'].number):
            del self.plots['Lagrange']
            refresh = False
        # Update the existing plot, if it exists
        refresh = refresh or 'Lagrange' in self.plots.keys()
        if refresh:
            if 'Lagrange' in self.plots.keys():
                fig = self.plots['Lagrange']
                fig = matplotlib.pyplot.figure(fig.number)
                fig.clear()
            else:
                return
        # Make a new window:
        else:
            fig = matplotlib.pyplot.figure(figsize=(4,3))
            fig.canvas.set_window_title('Lagrange')
        ax = fig.add_subplot(111)

        # Plot mass elements:
        for i in range(0, len(self.Lagrange[0]), int(float(self.incrementVar.get()))):
            ax.plot(1e9*self.time, 1e4*self.Lagrange[:,i], 'k-')
        # Plot fuel-shell interface:
        if self.fuelShellVar.get():
            ax.plot(1e9*self.time, 1e4*self.fuelShell, 'r-', linewidth=2)

        ax.set_xlim(float(self.minTVar.get()), float(self.maxTVar.get()))
        ax.set_ylim(float(self.minRVar.get()), float(self.maxRVar.get()))

        ax.set_xlabel('Time (ns)', fontsize=12)
        ax.set_ylabel('Radius (um)', fontsize=12)

        matplotlib.pyplot.tight_layout()

        if not refresh:
            if self.wm is not None:
                self.wm.addWindow(matplotlib.pyplot.get_current_fig_manager().window)
            fig.show()
            fig.canvas.draw()
        self.plots['Lagrange'] = fig
