# impy - a post-processor for HYADES implosion simulations
# Copyright (c) Massachusetts Institute of Technology / Alex Zylstra
# Distributed under the MIT License

import tkinter as tk
import tkinter.ttk as ttk
import platform

class Generate_Guderley(tk.Toplevel):
    """Implement a dialog window to prompt a user to configure a Guderley model

    :author: Alex Zylstra
    :date: 2014-05-15
    """
    __author__ = 'Alex Zylstra'
    __date__ = '2014-05-15'
    __version__ = '0.2.0'

    def __init__(self, parent, title='Configure Guderley'):
        """Initialize the dialog window"""
        super(Generate_Guderley, self).__init__(parent)
        self.transient(parent)
        self.parent = parent
        self.lift()
        self.grid()
        tk.Grid.columnconfigure(self, 0, weight=1)
        tk.Grid.columnconfigure(self, 1, weight=1)

        self.grab_set()

        self.result = None
        self.cancelled = False
        self.ions = ['', 'D', '3He', 'T']  # choices for type of fuel ion
        self.__create_widgets__(title)

        self.protocol("WM_DELETE_WINDOW", self.__cancel__)

        # a couple key bindings:
        self.bind('<Return>', self.__ok__)
        self.bind('<Escape>', self.__cancel__)

        # Set window background
        if platform.system() == 'Darwin':
            self.configure(background='#E8E9E8')
        else:
            self.configure(background='#F1F1F1')

        self.wait_window(self)

    def __create_widgets__(self, title):
        """Create the UI"""
        if title is not None:
            self.title(title)

        row = 0

        # input for initial radius
        label1 = ttk.Label(self, text='Radius (um):')
        label1.grid(row=row, column=0)
        self.radiusVar = tk.StringVar()
        self.radiusVar.set('430')
        radius = ttk.Entry(self, textvariable=self.radiusVar, width=8)
        radius.grid(row=row, column=1)
        row += 1

        # input for collapse time
        label2 = ttk.Label(self, text='tc (ns):')
        label2.grid(row=row, column=0)
        self.tcVar = tk.StringVar()
        self.tcVar.set('1.0')
        tc = ttk.Entry(self, textvariable=self.tcVar, width=8)
        tc.grid(row=row, column=1)
        row += 1

        # input for shock strength
        label2 = ttk.Label(self, text='xi (um/ns^α):')
        label2.grid(row=row, column=0)
        self.xiVar = tk.StringVar()
        self.xiVar.set('400')
        xi = ttk.Entry(self, textvariable=self.xiVar, width=8)
        xi.grid(row=row, column=1)
        row += 1

        # choose ion type 1
        label3 = ttk.Label(self, text='Ion 1:')
        label3.grid(row=row, column=0)
        self.Ion1Var = tk.StringVar()
        self.Ion1Var.set('D')
        Ion1 = ttk.OptionMenu(self, self.Ion1Var, *self.ions)
        Ion1.grid(row=row, column=1)
        row += 1

        # choose pressure of ion type #1
        label4 = ttk.Label(self, text='P (atm)')
        label4.grid(row=row, column=0)
        self.Ion1PVar = tk.StringVar()
        self.Ion1PVar.set('6')
        Ion1P = ttk.Entry(self, textvariable=self.Ion1PVar, width=8)
        Ion1P.grid(row=row, column=1)
        row += 1

        # choose ion type 2
        label4 = ttk.Label(self, text='Ion 2:')
        label4.grid(row=row, column=0)
        self.Ion2Var = tk.StringVar()
        self.Ion2Var.set('3He')
        Ion1 = ttk.OptionMenu(self, self.Ion2Var, *self.ions)
        Ion1.grid(row=row, column=1)
        row += 1

        # choose pressure of ion type #2
        label5 = ttk.Label(self, text='P (atm)')
        label5.grid(row=row, column=0)
        self.Ion2PVar = tk.StringVar()
        self.Ion2PVar.set('12')
        Ion2P = ttk.Entry(self, textvariable=self.Ion2PVar, width=8)
        Ion2P.grid(row=row, column=1)
        row += 1

        # choose e-i coupling
        label6 = ttk.Label(self, text='Rygg e-i coupling')
        label6.grid(row=row, column=0)
        self.CouplingVar = tk.StringVar()
        self.CouplingVar.set('0.5')
        Coupling = ttk.Entry(self, textvariable=self.CouplingVar, width=8)
        Coupling.grid(row=row, column=1)
        row += 1

        # choose radial step size
        label7 = ttk.Label(self, text='r step (um)')
        label7.grid(row=row, column=0)
        self.drVar = tk.StringVar()
        self.drVar.set('2.0')
        dr = ttk.Entry(self, textvariable=self.drVar, width=8)
        dr.grid(row=row, column=1)
        row += 1

        # choose time limits and step:
        label8 = ttk.Label(self, text='t_min (ns)')
        label8.grid(row=row, column=0)
        self.t0Var = tk.StringVar()
        self.t0Var.set('0.0')
        t0 = ttk.Entry(self, textvariable=self.t0Var, width=8)
        t0.grid(row=row, column=1)
        row += 1
        label9 = ttk.Label(self, text='t_max (ns)')
        label9.grid(row=row, column=0)
        self.t1Var = tk.StringVar()
        self.t1Var.set('2.0')
        t1 = ttk.Entry(self, textvariable=self.t1Var, width=8)
        t1.grid(row=row, column=1)
        row += 1
        label10 = ttk.Label(self, text='dt (ns)')
        label10.grid(row=row, column=0)
        self.dtVar = tk.StringVar()
        self.dtVar.set('0.01')
        dt = ttk.Entry(self, textvariable=self.dtVar, width=8)
        dt.grid(row=row, column=1)
        row += 1

        #Add the OK and cancel buttons
        box = ttk.Frame(self)
        w = ttk.Button(box, text="OK", width=10, command=self.__ok__)
        w.pack(side=tk.LEFT, padx=5, pady=5)
        w = ttk.Button(box, text="Cancel", width=10, command=self.__cancel__)
        w.pack(side=tk.LEFT, padx=5, pady=5)
        self.bind("<Return>", self.__ok__)
        self.bind("<Escape>", self.__cancel__)
        box.grid(row=row, columnspan=2)

    def __ok__(self, event=None):
        """Handle activation of the OK button."""
        if not self.__validate__():
            print('not valid')
            return

        self.__apply__()
        self.withdraw()
        self.update_idletasks()

        # put focus back to the parent window
        if self.parent is not None:
            self.parent.focus_set()
        self.destroy()

    def __cancel__(self, event=None):
        """Handle cancel button"""
        self.cancelled = True
        # put focus back to the parent window
        if self.parent is not None:
            self.parent.focus_set()
        self.destroy()

    def __validate__(self):
        """Validate the selection, returns true if it is OK"""
        try:
            r = float(self.radiusVar.get())
            tc = float(self.tcVar.get())
            xi = float(self.xiVar.get())
            p1 = float(self.Ion1PVar.get())
            p2 = float(self.Ion2PVar.get())
            c = float(self.CouplingVar.get())
            dr = float(self.drVar.get())
            t0 = float(self.t0Var.get())
            t1 = float(self.t1Var.get())
            dt = float(self.dtVar.get())

            # sanity checks:
            if r > 0 and xi > 0 and p1 > 0 and p2 > 0 and c >= 0 and c <= 1 \
                    and dr > 0 and dr < r and t0 < t1 and dt < (t1-t0):
                return True
        except:
            pass
        return False

    def __apply__(self):
        """Set the result"""
        r = float(self.radiusVar.get())
        tc = float(self.tcVar.get())
        xi = float(self.xiVar.get())
        ion1 = self.Ion1Var.get()
        p1 = float(self.Ion1PVar.get())
        ion2 = self.Ion2Var.get()
        p2 = float(self.Ion2PVar.get())
        c = float(self.CouplingVar.get())
        dr = float(self.drVar.get())
        t0 = float(self.t0Var.get())
        t1 = float(self.t1Var.get())
        dt = float(self.dtVar.get())
        self.result = [r,tc,xi,ion1,p1,ion2,p2,c,dr,t0,t1,dt]