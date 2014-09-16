# impy - a post-processor for HYADES implosion simulations
# Copyright (c) Massachusetts Institute of Technology / Alex Zylstra
# Distributed under the MIT License

import tkinter as tk
import tkinter.ttk as ttk
import math
import platform

class LILAC_Material(tk.Toplevel):
    """Implement a dialog window to prompt a user to configure LILAC materials.

    :author: Alex Zylstra
    :date: 2014-05-20
    """
    __author__ = 'Alex Zylstra'
    __date__ = '2014-05-20'
    __version__ = '0.2.0'

    def __init__(self, parent, material, zones, title='Configure LIlAC Material'):
        """Initialize the dialog window"""
        super(LILAC_Material, self).__init__(parent)
        self.transient(parent)
        self.parent = parent
        self.lift()
        self.grid()
        tk.Grid.columnconfigure(self, 0, weight=1)
        tk.Grid.columnconfigure(self, 1, weight=1)
        tk.Grid.columnconfigure(self, 2, weight=1)

        self.grab_set()

        self.result = None
        self.cancelled = False

        self.__create_widgets__(title, material, zones)

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

    def __create_widgets__(self, title, material, zones):
        """Create the UI"""
        if title is not None:
            self.title(title)

        row = 0

        # labels for top-level info
        label1 = ttk.Label(self, text='Material # ' + str(material))
        label1.grid(row=row, column=0, columnspan=3)
        row += 1

        label2 = ttk.Label(self, text='Zones ' + str(zones[0]) + '-' + str(zones[1]))
        label2.grid(row=row, column=0, columnspan=3)
        row += 1

        label3 = ttk.Label(self, text='# of ions')
        label3.grid(row=row, column=0, columnspan=2)
        self.numVar = tk.StringVar()
        self.numEntry = ttk.Entry(self, textvariable=self.numVar, width=6)
        self.numVar.set('1')
        self.numEntry.grid(row=row, column=2)
        row += 1

        label4 = ttk.Label(self, text='A')
        label4.grid(row=row, column=0)
        label5 = ttk.Label(self, text='Z')
        label5.grid(row=row, column=1)
        label6 = ttk.Label(self, text='F')
        label6.grid(row=row, column=2)
        row += 1

        self.ionFrame = ttk.Frame(self)
        self.ionFrame.grid(row=row, column=0, columnspan=3)
        self.__generateIonFrame__()
        self.numVar.trace('w', self.__generateIonFrame__)
        row += 1

        #Add the OK and cancel buttons
        box = ttk.Frame(self)
        w = ttk.Button(box, text="OK", width=10, command=self.__ok__)
        w.pack(side=tk.LEFT, padx=5, pady=5)
        w = ttk.Button(box, text="Cancel", width=10, command=self.__cancel__)
        w.pack(side=tk.LEFT, padx=5, pady=5)
        self.bind("<Return>", self.__ok__)
        self.bind("<Escape>", self.__cancel__)
        box.grid(row=row, columnspan=3)

    def __generateIonFrame__(self, *args):
        """Populate the entries for making ion frames"""
        try:
            n = int(self.numVar.get())
        except ValueError:
            return

        # variables for storing results
        self.ionAVars = []
        self.ionZVars = []
        self.ionFVars = []
        for child in self.ionFrame.winfo_children():
            child.destroy()
        #self.ionFrame.grid_forget()

        # Loop and make a UI element and var for each ion, and for each A/Z/F
        for i in range(n):
            AVar = tk.StringVar()
            AEntry = ttk.Entry(self.ionFrame, textvariable=AVar, width=6)
            AEntry.grid(row=i, column=0)
            self.ionAVars.append(AVar)

            ZVar = tk.StringVar()
            ZEntry = ttk.Entry(self.ionFrame, textvariable=ZVar, width=6)
            ZEntry.grid(row=i, column=1)
            self.ionZVars.append(ZVar)

            FVar = tk.StringVar()
            FEntry = ttk.Entry(self.ionFrame, textvariable=FVar, width=6)
            FEntry.grid(row=i, column=2)
            self.ionFVars.append(FVar)

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
            self.__apply__()
            # Sanity checks:
            # number of ions must be positive
            if self.result[0] <= 0:
                return False
            totalF = 0.
            # All A,Z,F must make sense:
            for i in range(self.result[0]):
                A = self.result[1][i]
                Z = self.result[2][i]
                F = self.result[3][i]
                if A <= 0 or Z <= 0 or A < Z or F <= 0 or F > 1:
                    return False
                totalF += F
            # Total F must add up to 1 (more or less):
            if math.fabs(totalF-1.) > 0.01:
                return False
        except:
            return False
        return True

    def __apply__(self):
        """Set the result"""
        n = int(self.numVar.get())
        A = []
        Z = []
        F = []
        for i in range(n):
            A.append( float(self.ionAVars[i].get()) )
            Z.append( float(self.ionZVars[i].get()) )
            F.append( float(self.ionFVars[i].get()) )
        self.result = [n, A, Z, F]