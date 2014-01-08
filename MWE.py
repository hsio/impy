# minimum working example.

from datetime import *

import tkinter as tk
import tkinter.ttk as ttk
from impy.gui.WindowManager import WindowManager

class Application(tk.Tk):
    """Analysis and database application for the NIF WRF data"""

    def __init__(self):
        super(Application, self).__init__(None)
        self.configure(background='#eeeeee')
        self.grid()
        self.createWidgets()
        self.minsize(150,200)
        self.title('IMPY')

        self.wm = WindowManager(self.winfo_screenwidth(), self.winfo_screenheight())
        self.wm.addWindow(self)

        # stretch the column to fill all space:
        tk.Grid.columnconfigure(self, 0, weight=1)
        tk.Grid.columnconfigure(self, 1, weight=1)
        tk.Grid.columnconfigure(self, 2, weight=1)

        # add a key binding to close:
        self.bind('<Escape>', self.quit)

    def createWidgets(self):
        self.DBInfoButton = ttk.Button(self, text="Open and run", command=self.run)
        self.DBInfoButton.grid(row=0, column=0)
        self.DBInfoButton = ttk.Button(self, text="Burn", command=self.showBurn)
        self.DBInfoButton.grid(row=1, column=0)

    def showBurn(self, *args):
        self.mod.display(type='GUI', wm=self.wm)

    def run(self, *args):
        t1 = datetime.now()
        print('Creating HYADES')
        from impy.implosions.Hyades import Hyades
        from tkinter.filedialog import askopenfilename
        FILEOPENOPTIONS = dict(defaultextension='.nc',
                       filetypes=[('Hyades netCDF','*.nc')],
                       multiple=False,
                       parent=None)
        filename = askopenfilename(parent=self)

        self.imp = Hyades(type='File', args=filename)
        self.imp.generate()

        t2 = datetime.now()
        print( '{:.2f}'.format((t2-t1).total_seconds()) + "s elapsed")


        t1 = datetime.now()
        print('Run fusion yield once:')
        from impy.modules.Burn import Burn
        self.mod = Burn()
        self.mod.run(self.imp)
        self.mod.display(type='CLI')
        t2 = datetime.now()
        print( '{:.2f}'.format((t2-t1).total_seconds()) + "s elapsed")
        t1 = datetime.now()


        t1 = datetime.now()
        print('Run fusion yield calculation 100x more:')
        for i in range(100):
            mod = Burn()
            mod.run(self.imp)
        t2 = datetime.now()
        print( '{:.2f}'.format((t2-t1).total_seconds()) + "s elapsed")

def main():
    app = Application()
    app.mainloop()

if __name__ == "__main__":
    main()