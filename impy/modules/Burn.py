
from impy.modules.Module import Module
from impy.implosions.Implosion import Implosion
from impy.resources.fusion import *

import tkinter as tk
import tkinter.ttk as ttk
from impy.gui.widgets.Table_View import Table_Viewer_Frame
import matplotlib, matplotlib.pyplot
import pickle


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
    :date: 2014-01-23
    """
    __author__ = 'Alex Zylstra'
    __date__ = '2014-01-23'
    __version__ = '1.0.0'

    # ----------------------------------------
    #           Generic methods
    # ----------------------------------------
    def __init__(self, type='GUI', args='', wm=None):
        """Construct a new instance of this module."""
        Module.__init__(self)

        #TODO: implement some options

        # set a few instance variables:
        self.runProgress = 0

        # Keep track of all plots created:
        self.plots = dict()

        # results:
        self.Y = dict()
        self.Ti = dict()
        self.burnRate = dict()
        self.bangTime = dict()
        self.burnRadius = dict()
        self.burnRadiusBins = dict()

    @classmethod
    def name(cls):
        """Get a string containing a name for this type of module."""
        return 'Burn'

    @classmethod
    def info(cls):
        """Get a brief description of this specific module."""
        return 'Fusion burn calculations'

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
        #TODO: burn width
        assert isinstance(imp, Implosion)
        self.runProgress = 0.

        # limits for calculations (all space/time):
        it = (imp.it_min(), imp.it_max())
        ir = (imp.ir_min(), imp.ir_max())
        # Ion temperature
        Ti = imp.Ti(it, ir)
        r = imp.r(it, ir)
        dt = imp.dt(it, ir)

        # Times (saved for plotting)
        self.time = imp.t(it)

        # Iterate over reactions:
        for rxn in allReactions():
            Y = imp.calcYield(it,ir,rxn)  # yield in each zone/time

            # save yield numbers:
            self.Y[rxn.name()] = np.sum(Y)

            # Do the rest of the calculations only for non-zero yield:
            if self.Y[rxn.name()] > 0:
                # Burn-weighted Ti:
                self.Ti[rxn.name()] = np.sum((Y*Ti) / self.Y[rxn.name()])

                # Calculate burn rate:
                self.burnRate[rxn.name()] = np.sum(Y/dt, axis=1)

                # Find bang time, defined as peak burn rate:
                itBT = np.argmax(self.burnRate[rxn.name()])
                self.bangTime[rxn.name()] = imp.t(itBT)*1e9  # also convert to ns

                # Calculate burn versus radius
                hist, bins = np.histogram(r.flatten(), bins=5000, weights=Y.flatten())
                self.burnRadius[rxn.name()] = hist
                dr = bins[1]-bins[0]
                self.burnRadiusBins[rxn.name()] = (bins[:-1]+dr/2)*1e4  # convert to um

            # Update progress:
            self.runProgress += 1./len(allReactions())

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
                    print('     Ti    = ' + '{:.2f}'.format(self.Ti[k]) + ' keV')
                    print('     BT    = ' + '{:.2f}'.format(self.bangTime[k]) + ' ns')

        elif refresh:
            self.__scalarViewer__(refresh=True)
            self.__burnRate__(refresh=True)
            self.__burnRadius__(refresh=True)

    def progress(self):
        """Get the calculation's progress estimate.

        :returns: Scalar number between 0 and 1.
        """
        return self.runProgress

    def copy(self, other):
        """Copy the results from `other` to this module.

        :param other: Another instance of this class
        """
        assert isinstance(other, Burn)
        try:
            self.time = np.copy(other.time)
            self.Y = other.Y.copy()
            self.Ti = other.Ti.copy()
            self.burnRate = other.burnRate.copy()
            self.bangTime = other.bangTime.copy()
            self.burnRadius = other.burnRadius.copy()
            self.burnRadiusBins = other.burnRadiusBins.copy()
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
            with open(filename, 'w') as file:
                csvwriter = csv.writer(file)
                file.write('Yields\n')
                for row in self.Y.items():
                    csvwriter.writerow(row)
                file.write('Ti (keV)\n')
                for row in self.Ti.items():
                    csvwriter.writerow(row)

            with open(filename, 'ab') as file:
                np.savetxt(file, self.time, header='Time:')
                for key in self.burnRate.keys():
                    np.savetxt(file, self.burnRate[key], header='Burn rate '+key)
                    np.savetxt(file, self.burnRadius[key], header='Burn radius '+key)
                    np.savetxt(file, self.burnRadiusBins[key], header='Burn radius bins '+key)

        elif type == 'pickle':
            with open(filename, 'wb') as file:
                pickle.dump(self.time, file)
                pickle.dump(self.Y, file)
                pickle.dump(self.Ti, file)
                pickle.dump(self.burnRate, file)
                pickle.dump(self.burnRadius, file)
                pickle.dump(self.burnRadiusBins, file)

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
        badKeys = ['wm','master','tk','_w','widgetName','plots','_tclCommands','_name','children','scalars']

        for key in badKeys:
            if key in state.keys():
                del state[key]
        return state

    # ----------------------------------------
    #           GUI stuff
    # ----------------------------------------
    def __createWidgets__(self):
        """Helper function which creates GUI elements."""
        scalarLabel = ttk.Label(self, text='Burn scalars')
        scalarLabel.grid(row=0, column=0, columnspan=2)

        self.__scalarViewer__()
        self.scalars.grid(row=1, column=0, columnspan=2, sticky='nsew')
        self.scalars.pack_propagate(0)
        self.scalars.grid_propagate(0)
        self.scalars.configure(width=250, height=100)

        plotLabel = ttk.Label(self, text='Plots')
        plotLabel.grid(row=2, column=0, columnspan=2, sticky='ns')

        burnRateButton = ttk.Button(self, text="Burn Rate", command=self.__burnRate__)
        burnRateButton.grid(row=3, column=0)

        burnRadiusButton = ttk.Button(self, text="Burn Radius", command=self.__burnRadius__)
        burnRadiusButton.grid(row=3, column=1)

    def __scalarViewer__(self, refresh=False, *args):
        """Helper function which is called to create the scalar viewer"""
        # Make a table to view the scalar results:
        columns = ('Reaction', 'Yield', 'Ti (keV)', 'BT (ns)')
        data = []
        for k in self.Y.keys():
            if self.Y[k] > 0:  # only display reactions with non-zero yield
                data.append( (k, '{:.2e}'.format(self.Y[k]), '{:.2f}'.format(self.Ti[k]), '{:.2f}'.format(self.bangTime[k])) )
        if not refresh:
            self.scalars = Table_Viewer_Frame(columns=columns, data=data, parent=self)
        else:
            self.scalars.setData(columns, data)

    def __burnRate__(self, refresh=False, *args):
        """Helper function to generate plot of burn rate"""
        # Update the existing plot, if it exists
        refresh = refresh or 'burnRate' in self.plots.keys()
        if refresh:
            if 'burnRate' in self.plots.keys():
                fig = self.plots['burnRate']
                fig = matplotlib.pyplot.figure(fig.number)
                fig.clear()
            else:
                return
        # Make a new window:
        else:
            fig = matplotlib.pyplot.figure(figsize=(4,3))
            fig.canvas.set_window_title('Burn Rate')
        ax = fig.add_subplot(111)

        for k in self.burnRate.keys():
            ax.plot(self.time, self.burnRate[k], label=k)

        ax.set_xlabel('Time (s)', fontsize=12)
        ax.set_ylabel('Burn Rate (1/s)', fontsize=12)
        ax.legend()

        matplotlib.pyplot.tight_layout()

        if not refresh:
            if self.wm is not None:
                self.wm.addWindow(matplotlib.pyplot.get_current_fig_manager().window)
            fig.show()
            fig.canvas.draw()
        self.plots['burnRate'] = fig

    def __burnRadius__(self, refresh=False, *args):
        """Helper function to generate plot of burn rate"""
        refresh = refresh or 'burnRadius' in self.plots.keys()
        if refresh:
            if 'burnRadius' in self.plots.keys():
                fig = self.plots['burnRadius']
                fig = matplotlib.pyplot.figure(fig.number)
                fig.clear()
            else:
                return
        # Make a new window:
        else:
            fig = matplotlib.pyplot.figure(figsize=(4,3))
            fig.canvas.set_window_title('Burn Radius')
        ax = fig.add_subplot(111)

        maxr = 0  # for plotting
        for k in self.burnRadius.keys():
            ax.plot(self.burnRadiusBins[k], self.burnRadius[k], label=k)
            maxr = max(maxr, self.__r1__(k))
        print(maxr)
        ax.set_xlim(0, maxr)
        ax.set_xlabel('Radius (um)', fontsize=12)
        ax.set_ylabel('Burn (1/um)', fontsize=12)
        ax.legend()

        matplotlib.pyplot.tight_layout()

        if not refresh:
            if self.wm is not None:
                self.wm.addWindow(matplotlib.pyplot.get_current_fig_manager().window)
            fig.show()
        self.plots['burnRadius'] = fig

    def __r1__(self, k):
        """Find the 10% burn radius for reaction k. Helper function for plotting"""
        maxVal = np.max(self.burnRadius[k])
        iMax = np.argmax(self.burnRadius[k])
        for i in range(iMax, len(self.burnRadius[k])):
            if self.burnRadius[k][i] < maxVal/100.:
                return self.burnRadiusBins[k][i]
        return self.burnRadiusBins[k][-1]