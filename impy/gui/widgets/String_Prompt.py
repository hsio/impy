__author__ = 'Alex Zylstra'
__date__ = '2014-01-06'
__version__ = '1.0.0'

import tkinter as tk
import tkinter.ttk as ttk

class String_Prompt(tk.Toplevel):
    """Implement a dialog window to prompt a user to input a string. The value can be retrieved by the `result` member::

        s = String_Prompt(...)
        s.result

    :param parent: The parent UI element
    :param title: (optional) A title to display on this window [default=None]
    :param text: (optional) Text to display next to the prompt [default=None]
    :param invalid: (optional) a list of invalid choices [default is empty list]

    :author: Alex Zylstra
    :date: 2014-01-06
    """

    def __init__(self, parent, title=None, text=None, default=None, invalid=None):
        """Initialize the dialog window"""
        super(String_Prompt, self).__init__(parent)
        self.transient(parent)
        self.parent = parent
        self.lift()

        self.grab_set()

        self.result = default
        self.invalid = invalid
        self.__create_widgets__(title, text, default)

        self.protocol("WM_DELETE_WINDOW", self.__cancel__)

        # a couple key bindings:
        self.bind('<Return>', self.__ok__)
        self.bind('<Escape>', self.__cancel__)

        self.configure(background='#eeeeee')

        self.wait_window(self)

    def __create_widgets__(self, title, text, default):
        """Create the UI"""
        if title is not None:
            self.title(title)

        if text is not None:
            label1 = ttk.Label(self, text=text)
            label1.pack()

        if default is not None:
            self.var = tk.StringVar(value=str(default))
            entry = ttk.Entry(self, textvariable=self.var)
            entry.pack()
            entry.focus_force()

        self.__make_buttons__()

    def __make_buttons__(self):
        """Add the OK and cancel buttons"""
        box = tk.Frame(self)

        w = ttk.Button(box, text="OK", width=10, command=self.__ok__)
        w.pack(side=tk.LEFT, padx=5, pady=5)
        w = ttk.Button(box, text="Cancel", width=10, command=self.__cancel__)
        w.pack(side=tk.LEFT, padx=5, pady=5)

        self.bind("<Return>", self.__ok__)
        self.bind("<Escape>", self.__cancel__)

        box.pack()

    def __ok__(self, event=None):
        """Handle activation of the OK button."""
        if not self.__validate__():
            print('not valid')
            return

        self.__apply__()
        self.withdraw()
        self.update_idletasks()

        self.__cancel__()

    def __cancel__(self, event=None):
        """Handle cancel button"""
        # put focus back to the parent window
        self.parent.focus_set()
        self.destroy()

    def __validate__(self):
        """Validate the selection, returns true if it is OK"""
        ret = (self.var.get() != '')
        # also check against list of invalid options if requested:
        if self.invalid is not None:
            for i in self.invalid:
                ret = ret and (self.var.get() != i)
        return ret

    def __apply__(self):
        """Set the result"""
        self.result = self.var.get()