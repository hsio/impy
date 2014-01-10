__author__ = 'Alex Zylstra'
__date__ = '2014-01-08'
__version__ = '1.0.0'

from impy.modules.Module import Module
from impy.implosions.Implosion import Implosion
from impy.resources.fusion import *

import tkinter as tk
import tkinter.ttk as ttk
from impy.gui.widgets.Table_View import Table_Viewer
import matplotlib, matplotlib.pyplot


class Burn(Module, tk.Toplevel):
    """Fusion burn module.

    :param type: The type of constructor to be used. Available options are:

    'GUI': The class should interact directly with the user to get required options (default)

    'CLI': Interact with user via CLI, or take info from args.

    :param args: Additional information, which depends on the type of constructor:

    type='GUI': unused

    type='CLI': A full list of arguments passed to the executable, to be interpreted as the module pleases.

    :param wm: A window manager to use for GUI windows during construction

    :author: Alex Zylstra
    :date: 2014-01-08
    """

    # ----------------------------------------
    #           Generic methods
    # ----------------------------------------
    def __init__(self, type='GUI', args='', wm=None):
        """Construct a new instance of this module."""
        Module.__init__(self)

        #TODO: implement some options

        # set a few instance variables:
        self.progress = 0

        # Keep track of all plots created:
        self.plots = []

        # results:
        self.Y = dict()
        self.Ti = dict()
        self.burnRate = dict()

    def __createWidgets__(self):
        pass

    @classmethod
    def name(cls):
        """Get a string containing a name for this type of module."""
        return 'Burn'

    def info(self):
        """Get a string of information about this specific module."""
        #TODO: implement more interesting info string
        return 'Burn TBD'


    # ----------------------------------------
    #       Execution and GUI control
    # ----------------------------------------
    def run(self, imp):
        """Run the calculation.

        :param imp: An `Implosion` object.
        """
        #TODO: Bang time, burn width
        #TODO: Store data for plotting
        assert isinstance(imp, Implosion)

        # limits for calculations (all space/time):
        it = (imp.it_min(), imp.it_max())
        ir = (imp.ir_min(), imp.ir_max())
        # Ion temperature
        Ti = imp.Ti(it, ir)

        # Times (saved for plotting)
        self.time = imp.t(it)

        # Iterate over reactions:
        for rxn in allReactions():
            Y = imp.calcYield(it,ir,rxn)  # yield in each zone/time

            # save yield numbers:
            self.Y[rxn.name()] = np.sum(Y)

            # Calculate burn-weighted Ti if yield is non-zero:
            if self.Y[rxn.name()] > 0:
                self.Ti[rxn.name()] = np.sum((Y*Ti) / self.Y[rxn.name()])
            else:
                self.Ti[rxn.name()] = 0

            # Calculate burn rate:
            if self.Y[rxn.name()] > 0:
                self.burnRate[rxn.name()] = np.sum(Y, axis=1)

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
            self.__createWidgets__()
            self.title('Burn')
            self.update_idletasks()
            if wm is not None:
                wm.addWindow(self)

            self.wm = wm
            self.__createWidgets__()

        elif type is 'CLI':
            # iterate over all reactions with results calculated:
            for k in self.Y.keys():
                if self.Y[k] > 0:  # only print reactions with non-zero yield
                    print(k)
                    print('     Yield = ' + '{:.2e}'.format(self.Y[k]))
                    print('     Ti    = ' + '{:.2f}'.format(self.Ti[k]))

        elif refresh:
            pass
        #TODO: implement refresh options for display

    def progress(self):
        """Get the calculation's progress estimate.

        :returns: Scalar number between 0 and 1.
        """
        #TODO: implement progress
        return self.progress

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
        #TODO: implement save functionality
        pass

    def savePlots(self, dir, prefix, type):
        """Save all generated plots.

        :param dir: The directory to save the plots to
        :param prefix: A prefix to append to any files generated by this module
        :param type: The type of plot to generate (e.g. 'eps')

        As a result, files generated in `dir` should have names like::

            prefix + my_name + '.' + type

        """
        #TODO: save plots
        pass

    # ----------------------------------------
    #           GUI stuff
    # ----------------------------------------
    def __createWidgets__(self):
        """Helper function which creates GUI elements."""
        #self.label1 = ttk.Label(self, text="Show")
        #self.label1.grid(row=0, column=0)
        scalarButton = ttk.Button(self, text="Scalars", command=self.__scalarViewer__)
        scalarButton.grid(row=0, column=0)

        burnRateButton = ttk.Button(self, text="Burn Rate", command=self.__burnRate__)
        burnRateButton.grid(row=1, column=0)

    def __scalarViewer__(self, *args):
        """Helper function which is called to create the scalar viewer"""
        # Make a table to view the scalar results:
        columns = ('Reaction', 'Yield', 'Ti (keV)')
        data = []
        for k in self.Y.keys():
            if self.Y[k] > 0:  # only display reactions with non-zero yield
                data.append( (k, '{:.2e}'.format(self.Y[k]), '{:.2f}'.format(self.Ti[k])) )
        #TODO: implement main GUI!
        TV = Table_Viewer(columns=columns, data=data, parent=self)
        TV.geometry('230x150+100+100')
        TV.update_idletasks()
        if self.wm is not None:
            self.wm.addWindow(TV)

    def __burnRate__(self, *args):
        """Helper function to generate plot of burn rate"""
        fig = matplotlib.pyplot.figure(figsize=(4,3))
        ax = fig.add_subplot(111)

        for k in self.burnRate.keys():
            ax.plot(self.time, self.burnRate[k], label=k)

        ax.set_xlabel('Time (s)', fontsize=10)
        ax.set_ylabel('Burn Rate (1/s)', fontsize=10)
        ax.legend()

        matplotlib.pyplot.tight_layout()

        if self.wm is not None:
            self.wm.addWindow(matplotlib.pyplot.get_current_fig_manager().window)

        fig.show()