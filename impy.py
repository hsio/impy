
import platform
import os
import pkgutil

import tkinter as tk
import tkinter.ttk as ttk
from tkinter.filedialog import askdirectory
from impy.gui.WindowManager import WindowManager
from impy.gui.widgets.Progress_Dialog import Progress_Dialog
from impy.gui.widgets.Option_Prompt import Option_Prompt
import matplotlib, matplotlib.pyplot

from impy.util import ImplosionRunner
from impy.util import ModuleRunner

import impy.implosions
import impy.modules
from impy.implosions.Implosion import *
from impy.modules.Module import *

import multiprocessing

# auto import implosions and modules:
for importer, modname, ispkg in pkgutil.iter_modules(impy.implosions.__path__):
    module = __import__('impy.implosions.'+modname, fromlist="dummy")
    print("Imported implosion: ", module)
for importer, modname, ispkg in pkgutil.iter_modules(impy.modules.__path__):
    module = __import__('impy.modules.'+modname, fromlist="dummy")
    print("Imported module: ", module)


class impy(tk.Toplevel):
    """Post-processor for hydrodynamic simulation results."""
    __author__ = 'Alex Zylstra'
    __date__ = '2014-01-25'
    __version__ = '1.0.0'

    def __init__(self):
        super(impy, self).__init__(None)
        self.configure(background='#eeeeee')
        self.grid()
        self.modControlVars = dict()
        self.modControlChecks = dict()
        self.modules = dict()
        self.modRedisplay = dict()
        self.processes = []
        self.windows = []
        self.__createWidgets__()
        self.minsize(200,200)
        self.title('impy')

        self.wm = WindowManager(self.winfo_screenwidth(), self.winfo_screenheight())
        self.wm.addWindow(self)

        # stretch the column to fill all space:
        tk.Grid.columnconfigure(self, 0, weight=1)
        tk.Grid.columnconfigure(self, 1, weight=1)
        tk.Grid.columnconfigure(self, 2, weight=1)

        # add a key binding to close:
        self.bind('<Escape>', self.close)
        self.protocol("WM_DELETE_WINDOW", self.close)

        # Text for shortcuts:
        if platform.system() == 'Darwin':
            shortcutType = '⌘'
            shortcutModifier = 'Command-'
        else:
            shortcutType = 'Ctrl+'
            shortcutModifier = 'Control-'

        # Top-level menu bar:
        menubar = tk.Menu(self)
        # File menu:
        fileMenu = tk.Menu(menubar, tearoff=0)
        menubar.add_cascade(label='File', menu=fileMenu)
        helpMenu = tk.Menu(menubar, tearoff=0)
        menubar.add_cascade(label='Help', menu=helpMenu)
        self.config(menu=menubar)

        fileMenu.add_command(label='Open File\t\t' + shortcutType + 'O', command=self.openFile)
        self.bind('<' + shortcutModifier + 'o>', self.openFile)

        createMenu = tk.Menu(fileMenu, tearoff=0)
        fileMenu.add_cascade(label='Create...', menu=createMenu)

        # Auto generate creation options for each implosion type:
        for impType in allImplosions():
            createMenu.add_command(label='Create ' + impType.name(), command= lambda: self.open(impType.name()))

        fileMenu.add_command(label='Save\t\t\t' + shortcutType + 'S', command= lambda: self.save('CSV'))
        self.bind('<' + shortcutModifier + 's>', lambda *args: self.save('CSV'))
        fileMenu.add_command(label='Save pickle', command= lambda: self.save('pickle'))
        fileMenu.add_command(label='Save plots', command= lambda: self.save('plots'))
        fileMenu.add_command(label='Quit\t\t\t' + shortcutType + 'Q', command=self.close)
        self.bind('<' + shortcutModifier + 'q>', self.close)

        # Options in the help menu:
        helpMenu.add_command(label='About', command=self.__about__),

        self.__configureMatplotlib__()

    def __createWidgets__(self):
        """Create the UI elements for the main window"""
        label1 = ttk.Label(self, text='Type:')
        label1.grid(row=0, column=0)
        self.typeLabelVar = tk.StringVar()
        self.typeLabel = ttk.Label(self, textvar=self.typeLabelVar)
        self.typeLabel.grid(row=0, column=1, sticky='ns')
        impInfoButton = ttk.Button(self, text='Info', command=self.__impInfo__, width=4)
        impInfoButton.grid(row=0, column=2, sticky='ns')

        label2 = ttk.Label(self, text='File:')
        label2.grid(row=1, column=0)
        self.fileLabelVar = tk.StringVar()
        self.fileLabel = ttk.Label(self, textvar=self.fileLabelVar)
        self.fileLabel.grid(row=1, column=1, sticky='ns')

        bar1 = ttk.Separator(self)
        bar1.grid(row=2, column=0, columnspan=3, sticky='nsew')

        label3 = ttk.Label(self, text='Modules')
        label3.grid(row=3, column=0, columnspan=3, sticky='ns')
        # Add for each module:
        row=4
        for mod in allModules():
            label = ttk.Label(self, text=mod.name())
            label.grid(row=row, column=0)
            self.modControlVars[mod.name()] = tk.IntVar()
            check = ttk.Checkbutton(self, variable=self.modControlVars[mod.name()], command= lambda a=mod.name():self.__runModule__(a), state=tk.DISABLED)
            check.grid(row=row, column=1)
            self.modControlChecks[mod.name()] = check
            infoButton = ttk.Button(self, text='Info', command= lambda a=mod.name(): self.__modInfo__(a), width=4)
            infoButton.grid(row=row, column=2)
            row += 1

    def __configureMatplotlib__(self):
        # set matplotlib backend
        if matplotlib.get_backend() != 'tkagg':
            matplotlib.pyplot.switch_backend('TkAgg')
        matplotlib.pyplot.rc('font', **{'size':'8'})
        matplotlib.pyplot.rc('text', **{'usetex':False})
        matplotlib.rcParams['toolbar'] = 'None'

    def close(self, *args):
        """Handle closing the application."""
        matplotlib.pyplot.close("all")
        for p in self.processes:
            if p.is_alive():
                p.terminate()
        for w in self.windows:
            w.withdraw()
        self.withdraw()
        self.quit()

    def openFile(self, *args):
        """Open an implosion file. Attempts to handle conflicts (two implosions can open same type) gracefully."""
        # Get file types:
        filetypes = []
        correspondingImp = []
        for impType in allImplosions():
            for type in impType.getFileTypes():
                filetypes.append(type)
                correspondingImp.append(impType)

        # Prompt the user:
        from tkinter.filedialog import askopenfilename
        FILEOPENOPTIONS = dict(defaultextension='.nc',
                       filetypes=filetypes,
                       multiple=False,
                       parent=self)
        filename = askopenfilename(**FILEOPENOPTIONS)

        if filename == '':
            return

        # Create the implosion:
        self.filename = filename
        ext = '*.' + filename.split('.')[-1]
        # Find out which types this corresponds to:
        typeIndex = []
        for i in range(len(filetypes)):
            if filetypes[i][1] == ext:
                typeIndex.append(i)
        # If none:
        if len(typeIndex) == 0:
            return

        # Detect 'collisions', i.e. multiple implosions can open this file type
        if len(typeIndex) > 1:
            # Prompt the user to choose:
            from impy.gui.widgets.Option_Prompt import Option_Prompt
            opts = []
            for i in typeIndex:
                opts.append(correspondingImp[i].name())
            p = Option_Prompt(self, title='Choose implosion type', options=opts)

            # Figure out which was chosen:
            chosenName = p.result
            if p.result is None:
                return
            for i in typeIndex:
                if correspondingImp[i].name() == p.result:
                    typeIndex = i
                    break

        # If only one option, convert from list:
        else:
            typeIndex = typeIndex[0]

        # Now one is chosen:
        self.impType = correspondingImp[typeIndex].name()
        self.imp = correspondingImp[typeIndex](type='File', args=filename)

        # Set info to display:
        self.typeLabelVar.set(self.impType)
        self.fileLabelVar.set(os.path.split(self.filename)[-1])

        self.__runImplosion__()

    def open(self, type, *args):
        """Create a specified implosion using its own constructor.

        :param type: A string containing the name of the implosion class.
        """
        # Make the correct type:
        for impType in allImplosions():
            if impType.name() == type:
                self.impType = impType.name()
                self.imp = impType(type='GUI')
                break

        self.typeLabelVar.set(self.impType)
        try:
            self.fileLabelVar.set(os.path.split(self.imp.filename)[-1])
        except:
            self.fileLabelVar.set('')

        self.__runImplosion__()

    def save(self, saveType, *args):
        """Save results from all open modules to a directory.

        :param saveType: The type of thing to save: 'plots', 'csv', or 'pickle'.
        """
        # Get a directory to save to:
        dir = askdirectory(parent=self)
        if dir is None or dir == '':
            return
        # Create it if necessary:
        if not os.path.dirname(dir):
            os.makedirs(dir)

        # For saving plots, need to get the type:
        if saveType == 'plots':
            formats = ['']+list(matplotlib.pyplot.gcf().canvas.get_supported_filetypes().keys())
            prompt = Option_Prompt(self, title='Plot format', options=formats)

        # Loop over modules and save their stuff if the checkbox is selected:
        for mod in self.modules.keys():
            if self.modControlVars[mod].get() == 1:
                if saveType == 'plots':
                    if prompt.result is not None:
                        prefix = mod
                        self.modules[mod].savePlots(dir, prefix, prompt.result)
                else:
                    fname = os.path.join(dir, mod + '.' + saveType)
                    self.modules[mod].save(fname, type=saveType)

    def __runImplosion__(self):
        """Generate the implosion. Computation work done in a different process"""
        # Create a progress dialog:
        dialog = Progress_Dialog()
        dialog.set_text('Implosion generation...')
        dialog.update_idletasks()

        # Use helper function:
        parent_conn, child_conn = multiprocessing.Pipe()
        p = multiprocessing.Process(target=ImplosionRunner.run, args=(child_conn,))
        self.processes.append(p)

        # start the process and send implosion:
        p.start()
        try:
            parent_conn.send(self.imp)
        except:
            raise Exception('Implosion object passed to ImplosionRunner is not pickleable!')

        obj = None
        # Loop while the process is active:
        def callback():
            nonlocal dialog, p, parent_conn
            if dialog.cancelled:
                dialog.withdraw()
                p.terminate()
                return

            # Try to receive from the Pipe:
            if parent_conn.poll():
                # Update the progress, or end otherwise:
                obj = parent_conn.recv()
                if isinstance(obj, Exception):
                    from tkinter.messagebox import showerror
                    showerror('Error!', 'A problem occurred generating the implosion (class '+self.imp.name()+')\n'+obj.__str__())
                    dialog.withdraw()
                    p.terminate()
                    return
                elif isinstance(obj, float):
                    dialog.set(100*obj)
                elif isinstance(obj, Implosion):
                    # Pass info back to the main app:
                    self.imp = obj
                    self.after(10, self.__postImplosion__)
                    dialog.withdraw()
                    p.terminate()
                    return

            self.after(25, callback)

        self.after(10, callback)

    def __postImplosion__(self):
        """Tasks to execute after the implosion is generated."""
        # Enable modules:
        for key in self.modControlChecks.keys():
            self.modControlChecks[key].configure(state=tk.NORMAL)

        # Run any modules that were already checked in refresh mode
        for key in self.modControlChecks.keys():
            self.modRedisplay[key] = (self.modControlVars[key].get() == 1)
        for mod in allModules():
            if self.modRedisplay[mod.name()]:
                self.__runModule__(mod.name())

    def __runModule__(self, modName):
        """Run a specified module.

        :param mod: The name of the module class to run
        """
        mod = None
        for x in allModules():
            if x.name() == modName:
                mod = x

        if mod is None:
            return

        # Check whether we are loading or unloading:
        if self.modControlVars[mod.name()].get() == 1:
            # Prepare the module if not done already:
            if mod.name() not in self.modules.keys():
                self.modules[mod.name()] = mod(type='GUI')

            # Create a progress bar:
            dialog = Progress_Dialog()
            dialog.set_text('Running ' + mod.name() + '...')
            dialog.update_idletasks()

            # Run using helper function:
            parent_conn, child_conn = multiprocessing.Pipe()

            p = multiprocessing.Process(target=ModuleRunner.run, args=(child_conn,))
            self.processes.append(p)

            p.start()
            try:
                parent_conn.send(self.modules[mod.name()])
            except:
                print(mod.__dict__)
                raise Exception('Module object passed to ModuleRunner is not pickle-able!')

            try:
                parent_conn.send(self.imp)
            except:
                raise Exception('Implosion object passed to ModuleRunner is not pickle-able!')

            obj = None
            # Callpacks while the process is active:
            def callback():
                nonlocal p, dialog, parent_conn
                if dialog.cancelled:
                    p.terminate()
                    dialog.withdraw()
                    return

                # Try to receive from the pipe:
                if parent_conn.poll():
                    # Update the progress, or end otherwise:
                    obj = parent_conn.recv()
                    if isinstance(obj, Exception):
                        from tkinter.messagebox import showerror
                        showerror('Error!', 'A problem occurred running '+mod.name()+'\n'+obj.__str__())
                        # Uncheck from main window:
                        self.modControlVars[mod.name()].set(0)
                        dialog.withdraw()
                        p.terminate()
                        return
                    elif isinstance(obj, float) or isinstance(obj, int):
                        dialog.set(100*obj)
                    elif isinstance(obj, Module):
                        # Pass back to the main app:
                        self.modules[mod.name()].copy(obj)
                        self.after(10, lambda: self.__postModule__(obj))
                        dialog.withdraw()
                        p.terminate()
                        return  # Callback loop breaks here

                self.after(25, callback)  # loop

            # Start callback loop in a bit:
            self.after(10, callback)

        else:
            # If a module is unchecked, it is removed and deleted:
            self.__deleteModule__(mod)

    def __postModule__(self, mod):
        """Handle tasks to be done after completion of a module run."""
        self.modules[mod.name()].display(wm=self.wm, refresh=self.modRedisplay[mod.name()])
        try:
            self.modules[mod.name()].protocol(name='WM_DELETE_WINDOW', func=lambda: self.__deleteModule__(mod))
        except:
            pass

    def __deleteModule__(self, mod):
        """Delete a module"""
        if mod.name() in self.modules.keys():
            self.modules[mod.name()].withdraw()
            # remove all plots generated by that module:
            for fig in self.modules[mod.name()].getPlots():
                matplotlib.pyplot.close(fig)
            self.modules.__delitem__(mod.name())
            self.modControlVars[mod.name()].set(0)

    def __impInfo__(self):
        """Show some info about the current implosion type."""
        title = ''
        text = ''
        try:
            title = self.impType + ' Information'
            # Identify the implosion type:
            for impType in allImplosions():
                if impType.name() == self.impType:
                    text += 'Author: ' + impType.__author__ + '\n'
                    text += 'Date: ' + impType.__date__ + '\n'
                    text += 'Version: ' + impType.__version__ + '\n'
                    text += '\n'
                    text += self.imp.info()
        except:
            pass
        if text != '' and title != '':
            from tkinter.messagebox import showinfo
            showinfo(title=title, message=text)

    def __modInfo__(self, type):
        """Show some info about a selected module."""
        title = ''
        text = ''
        try:
            title = type + ' Information'
            mod = self.modules[type]
            text += 'Author: ' + mod.__author__ + '\n'
            text += 'Date: ' + mod.__date__ + '\n'
            text += 'Version: ' + mod.__version__ + '\n'
            text += '\n'
            text += mod.detailedInfo()
        except:
            pass
        if text != '' and title != '':
            from tkinter.messagebox import showinfo
            showinfo(title=title, message=text)

    def __about__(self):
        from tkinter.messagebox import showinfo
        title = 'impy'
        text = 'Author: ' + self.__author__ + '\n'
        text += 'Date: ' + self.__date__ + '\n'
        text += 'Version: ' + self.__version__ + '\n'
        showinfo(title=title, message=text)


if __name__ == "__main__":
    multiprocessing.freeze_support()
    root = tk.Tk()
    root.withdraw()
    app = impy()
    root.mainloop()