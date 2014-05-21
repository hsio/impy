
from impy.modules.Module import Module
from impy.implosions.Implosion import Implosion
from impy.resources.fusion import *

import tkinter as tk
import tkinter.ttk as ttk
import matplotlib, matplotlib.pyplot
import math
import platform


class Snapshot(Module, tk.Toplevel):
    """Snapshot (i.e. hydro variables at given time) plotting module.

    :param type: The type of constructor to be used. Available options are:

    'GUI': The class should interact directly with the user to get required options (default)

    'CLI': Interact with user via CLI, or take info from args.

    :param args: Additional information, which depends on the type of constructor:

    type='GUI': unused

    type='CLI': A full list of arguments passed to the executable, to be interpreted as the module pleases.

    :param wm: A window manager to use for GUI windows during construction

    :author: Alex Zylstra
    :date: 2014-02-12
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

    @classmethod
    def name(cls):
        """Get a string containing a name for this type of module."""
        return 'Snapshot'

    @classmethod
    def info(cls):
        """Get a brief description of this specific module."""
        return 'Snapshot (i.e. hydro variables at given time) plotting module.'

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
        self.imp = imp
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
            self.title('Snapshot')
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
        assert isinstance(other, Snapshot)
        self.imp = other.imp

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
        dt = self.imp.dt(0)
        it = int(self.imp.it_min() + ( float(self.timeVar.get())*1e-9 - self.imp.t(self.imp.it_min())) / dt)
        ir = (self.imp.ir_min(), self.imp.ir_max())

        results = np.zeros((ir[1]-ir[0],8), dtype=np.float64)
        results[:,0] = self.imp.r((it), ir)
        results[:,1] = self.imp.rho((it), ir)
        results[:,2] = self.imp.P((it), ir)
        results[:,3] = self.imp.u((it), ir)
        results[:,4] = self.imp.ne((it), ir)
        results[:,5] = self.imp.ni((it), ir)
        results[:,6] = self.imp.Te((it), ir)
        results[:,7] = self.imp.Ti((it), ir)
        header = 'time = '+str(self.imp.t(it)) + '\n' + 'r (cm), rho (g/cc), P (GBar), u (cm/s), ne (1/cc), ni (1/cc), Te (keV), Ti (keV)'
        if type == 'CSV':
            np.savetxt(filename, results, delimiter=',', header=header)

        elif type == 'pickle':
            results.dump(filename)

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
                   'plotRhoVar', 'plotPressureVar', 'plotVelocityVar', 'plotneVar', 'plotniVar', 'plotTeVar', 'plotTiVar',
                   'timeVar', 'logxVar', 'logyVar']

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
        plotLabel.grid(row=1, column=0, columnspan=2, sticky='ns')

        label1 = ttk.Label(self, text='ρ')
        label1.grid(row=2, column=0)
        self.plotRhoVar = tk.BooleanVar(value=True)
        plotRhoCheck = ttk.Checkbutton(self, variable=self.plotRhoVar)
        plotRhoCheck.grid(row=2, column=1)

        label2 = ttk.Label(self, text='P')
        label2.grid(row=3, column=0)
        self.plotPressureVar = tk.BooleanVar(value=True)
        plotPressureCheck = ttk.Checkbutton(self, variable=self.plotPressureVar)
        plotPressureCheck.grid(row=3, column=1)

        label3 = ttk.Label(self, text='u')
        label3.grid(row=4, column=0)
        self.plotVelocityVar = tk.BooleanVar(value=True)
        plotVelocityCheck = ttk.Checkbutton(self, variable=self.plotVelocityVar)
        plotVelocityCheck.grid(row=4, column=1)

        label4 = ttk.Label(self, text='ne')
        label4.grid(row=5, column=0)
        self.plotneVar = tk.BooleanVar(value=True)
        plotneCheck = ttk.Checkbutton(self, variable=self.plotneVar)
        plotneCheck.grid(row=5, column=1)

        label5 = ttk.Label(self, text='ni')
        label5.grid(row=6, column=0)
        self.plotniVar = tk.BooleanVar(value=True)
        plotniCheck = ttk.Checkbutton(self, variable=self.plotniVar)
        plotniCheck.grid(row=6, column=1)

        label6 = ttk.Label(self, text='Te')
        label6.grid(row=7, column=0)
        self.plotTeVar = tk.BooleanVar(value=True)
        plotTeCheck = ttk.Checkbutton(self, variable=self.plotTeVar)
        plotTeCheck.grid(row=7, column=1)

        label7 = ttk.Label(self, text='Ti')
        label7.grid(row=8, column=0)
        self.plotTiVar = tk.BooleanVar(value=True)
        plotTiCheck = ttk.Checkbutton(self, variable=self.plotTiVar)
        plotTiCheck.grid(row=8, column=1)

        label8 = ttk.Label(self, text='t (ns)')
        label8.grid(row=9, column=0)
        self.timeVar = tk.StringVar(value=0)
        timeEntry = ttk.Entry(self, textvariable=self.timeVar, width=8)
        timeEntry.grid(row=9, column=1)

        split1 = ttk.Separator(self)
        split1.grid(row=10, column=0, columnspan=2, sticky='nsew')

        label9 = ttk.Label(self, text='Log x')
        label9.grid(row=11, column=0)
        self.logxVar = tk.BooleanVar(value=False)
        logxCheck = ttk.Checkbutton(self, variable=self.logxVar)
        logxCheck.grid(row=11, column=1)

        label9 = ttk.Label(self, text='Log y')
        label9.grid(row=12, column=0)
        self.logyVar = tk.BooleanVar(value=False)
        logyCheck = ttk.Checkbutton(self, variable=self.logyVar)
        logyCheck.grid(row=12, column=1)

        split2 = ttk.Separator(self)
        split2.grid(row=13, column=0, columnspan=2, sticky='nsew')

        burnRateButton = ttk.Button(self, text='Plot', command=self.__plot__)
        burnRateButton.grid(row=14, column=0, columnspan=2)

    def __plot__(self, refresh=False, *args):
        # Set up time and spatial steps:
        dt = self.imp.dt(0)
        self.it = self.imp.it_min() + math.ceil(( float(self.timeVar.get())*1e-9 - self.imp.t(self.imp.it_min())) / dt)
        self.ir = (self.imp.ir_min(), self.imp.ir_max())

        # call all plot functions:
        self.__plot_rho__(refresh, *args)
        self.__plot_pres__(refresh, *args)
        self.__plot_u__(refresh, *args)
        self.__plot_n__(refresh, *args)
        self.__plot_T__(refresh, *args)

    def __plot_rho__(self, refresh=False, *args):
        """Helper function to generate plot of density"""
        # If plot is not requested, return:
        if not self.plotRhoVar.get():
            return

        # Check for a closed window:
        if 'rho' in self.plots.keys() and not matplotlib.pyplot.fignum_exists(self.plots['rho'].number):
            del self.plots['rho']
            refresh = False
        # Update the existing plot, if it exists
        refresh = refresh or 'rho' in self.plots.keys()
        if refresh:
            if 'rho' in self.plots.keys():
                fig = self.plots['rho']
                fig = matplotlib.pyplot.figure(fig.number)
                fig.clear()
            else:
                return
        # Make a new window:
        else:
            fig = matplotlib.pyplot.figure(figsize=(4,3))
        fig.canvas.set_window_title('rho, time = ' + '{:.3f}'.format(1e9*self.imp.t(self.it)))
        ax = fig.add_subplot(111)

        # Plot:
        ax.plot(1e4*self.imp.r((self.it), self.ir)[0], self.imp.rho((self.it), self.ir)[0], 'k-')

        ax.set_xlabel('r (um)', fontsize=12)
        ax.set_ylabel('Rho (g/cc)', fontsize=12)

        if self.logxVar.get():
            ax.set_xscale('log')
        if self.logyVar.get():
            ax.set_yscale('log')

        matplotlib.pyplot.tight_layout()

        if not refresh:
            fig.show()
            fig.canvas.draw()
            if self.wm is not None:
                self.wm.addWindow(matplotlib.pyplot.get_current_fig_manager().window)
        self.plots['rho'] = fig

    def __plot_pres__(self, refresh=False, *args):
        """Helper function to generate plot of pressure"""
        # If plot is not requested, return:
        if not self.plotPressureVar.get():
            return

        # Check for a closed window:
        if 'pressure' in self.plots.keys() and not matplotlib.pyplot.fignum_exists(self.plots['pressure'].number):
            del self.plots['pressure']
            refresh = False
        # Update the existing plot, if it exists
        refresh = refresh or 'pressure' in self.plots.keys()
        if refresh:
            if 'pressure' in self.plots.keys():
                fig = self.plots['pressure']
                fig = matplotlib.pyplot.figure(fig.number)
                fig.clear()
            else:
                return
        # Make a new window:
        else:
            fig = matplotlib.pyplot.figure(figsize=(4,3))
        fig.canvas.set_window_title('pressure, time = ' + '{:.3f}'.format(1e9*self.imp.t(self.it)))
        ax = fig.add_subplot(111)

        # Plot:
        ax.plot(1e4*self.imp.r((self.it), self.ir)[0], self.imp.P((self.it), self.ir)[0], 'k-')

        ax.set_xlabel('r (um)', fontsize=12)
        ax.set_ylabel('Pressure (GBar)', fontsize=12)

        if self.logxVar.get():
            ax.set_xscale('log')
        if self.logyVar.get():
            ax.set_yscale('log')

        matplotlib.pyplot.tight_layout()

        if not refresh:
            fig.show()
            fig.canvas.draw()
            if self.wm is not None:
                self.wm.addWindow(matplotlib.pyplot.get_current_fig_manager().window)
        self.plots['pressure'] = fig

    def __plot_u__(self, refresh=False, *args):
        """Helper function to generate plot of pressure"""
        # If plot is not requested, return:
        if not self.plotVelocityVar.get():
            return

        # Check for a closed window:
        if 'u' in self.plots.keys() and not matplotlib.pyplot.fignum_exists(self.plots['u'].number):
            del self.plots['u']
            refresh = False
        # Update the existing plot, if it exists
        refresh = refresh or 'u' in self.plots.keys()
        if refresh:
            if 'u' in self.plots.keys():
                fig = self.plots['u']
                fig = matplotlib.pyplot.figure(fig.number)
                fig.clear()
            else:
                return
        # Make a new window:
        else:
            fig = matplotlib.pyplot.figure(figsize=(4,3))
        fig.canvas.set_window_title('u, time = ' + '{:.3f}'.format(1e9*self.imp.t(self.it)))
        ax = fig.add_subplot(111)

        # Plot:
        ax.plot(1e4*self.imp.r((self.it), self.ir)[0], 1e-5*self.imp.u((self.it), self.ir)[0], 'k-')

        ax.set_xlabel('r (um)', fontsize=12)
        ax.set_ylabel('Velocity (km/s)', fontsize=12)

        if self.logxVar.get():
            ax.set_xscale('log')
        if self.logyVar.get():
            ax.set_yscale('log')

        matplotlib.pyplot.tight_layout()

        if not refresh:
            fig.show()
            fig.canvas.draw()
            if self.wm is not None:
                self.wm.addWindow(matplotlib.pyplot.get_current_fig_manager().window)
        self.plots['u'] = fig

    def __plot_n__(self, refresh=False, *args):
        """Helper function to generate plot of burn rate"""
        # If plot is not requested, return:
        if not self.plotneVar.get() or not self.plotniVar.get():
            return

        # Check for a closed window:
        if 'n' in self.plots.keys() and not matplotlib.pyplot.fignum_exists(self.plots['n'].number):
            del self.plots['n']
            refresh = False
        # Update the existing plot, if it exists
        refresh = refresh or 'n' in self.plots.keys()
        if refresh:
            if 'n' in self.plots.keys():
                fig = self.plots['n']
                fig = matplotlib.pyplot.figure(fig.number)
                fig.clear()
            else:
                return
        # Make a new window:
        else:
            fig = matplotlib.pyplot.figure(figsize=(4,3))
        fig.canvas.set_window_title('n, time = ' + '{:.3f}'.format(1e9*self.imp.t(self.it)))
        ax = fig.add_subplot(111)

        # Plot:
        if self.plotneVar.get():
            ax.plot(1e4*self.imp.r((self.it), self.ir)[0], self.imp.ne((self.it), self.ir)[0], 'r-', label='e')
        if self.plotniVar.get():
            ax.plot(1e4*self.imp.r((self.it), self.ir)[0], self.imp.ni((self.it), self.ir)[0], 'b-', label='i')

        ax.set_xlabel('r (um)', fontsize=12)
        ax.set_ylabel('n (1/cc)', fontsize=12)
        ax.legend()

        if self.logxVar.get():
            ax.set_xscale('log')
        if self.logyVar.get():
            ax.set_yscale('log')

        matplotlib.pyplot.tight_layout()

        if not refresh:
            fig.show()
            fig.canvas.draw()
            if self.wm is not None:
                self.wm.addWindow(matplotlib.pyplot.get_current_fig_manager().window)
        self.plots['n'] = fig

    def __plot_T__(self, refresh=False, *args):
        """Helper function to generate plot of temperature"""
        # If plot is not requested, return:
        if not self.plotTeVar.get() or not self.plotTiVar.get():
            return

        # Check for a closed window:
        if 'T' in self.plots.keys() and not matplotlib.pyplot.fignum_exists(self.plots['T'].number):
            del self.plots['T']
            refresh = False
        # Update the existing plot, if it exists
        refresh = refresh or 'T' in self.plots.keys()
        if refresh:
            if 'T' in self.plots.keys():
                fig = self.plots['T']
                fig = matplotlib.pyplot.figure(fig.number)
                fig.clear()
            else:
                return
        # Make a Tew window:
        else:
            fig = matplotlib.pyplot.figure(figsize=(4,3))
        fig.canvas.set_window_title('T, time = ' + '{:.3f}'.format(1e9*self.imp.t(self.it)))
        ax = fig.add_subplot(111)

        # Plot:
        if self.plotTeVar.get():
            ax.plot(1e4*self.imp.r((self.it), self.ir)[0], self.imp.Te((self.it), self.ir)[0], 'r-', label='e')
        if self.plotTiVar.get():
            ax.plot(1e4*self.imp.r((self.it), self.ir)[0], self.imp.Ti((self.it), self.ir)[0], 'b-', label='i')

        ax.set_xlabel('r (um)', fontsize=12)
        ax.set_ylabel('T (keV)', fontsize=12)
        ax.legend()

        if self.logxVar.get():
            ax.set_xscale('log')
        if self.logyVar.get():
            ax.set_yscale('log')

        matplotlib.pyplot.tight_layout()

        if not refresh:
            fig.show()
            fig.canvas.draw()
            if self.wm is not None:
                self.wm.addWindow(matplotlib.pyplot.get_current_fig_manager().window)
        self.plots['T'] = fig