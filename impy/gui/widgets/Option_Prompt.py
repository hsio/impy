
import tkinter as tk
import tkinter.ttk as ttk

class Option_Prompt(tk.Toplevel):
    """Implement a dialog window to prompt a user to select one of several options. The chosen result
    is accessible via the member variable `result`::

        o = Option_Prompt(...)
        o.result

    :param parent: The parent UI element
    :param title: (optional) A title to display on this window [default=None]
    :param text: (optional) Text to display next to the prompt [default=None]
    :param options: (optional) The list of options that the user can choose from. First element must be ''
    :param width: (optional) the width in characters for the drop-down menu [default=10]

    :author: Alex Zylstra
    :date: 2014-01-25
    """
    __author__ = 'Alex Zylstra'
    __date__ = '2014-01-25'
    __version__ = '1.0.0'

    def __init__(self, parent, title=None, text=None, options=[], width=10):
        """Initialize the dialog window"""
        super(Option_Prompt, self).__init__(parent)
        self.transient(parent)
        self.parent = parent
        self.lift()

        self.grab_set()

        self.result = None
        self.__create_widgets__(title, text, options, width)

        self.protocol("WM_DELETE_WINDOW", self.__cancel__)

        # a couple key bindings:
        self.bind('<Return>', self.__ok__)
        self.bind('<Escape>', self.__cancel__)

        self.configure(background='#eeeeee')

        self.wait_window(self)

    def __create_widgets__(self, title, text, options, width):
        """Create the UI"""
        if title is not None:
            self.title(title)

        if text is not None:
            label1 = ttk.Label(self, text=text)
            label1.pack()

        if options is not None:
            self.var = tk.StringVar('')
            menu = ttk.OptionMenu(self, self.var, *options)
            menu.configure(width=width)
            menu.pack()
            menu.focus_force()

        self.__make_buttons__()

    def __make_buttons__(self):
        """Add the OK and cancel buttons"""
        box = ttk.Frame(self)

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
        return self.var.get() != ''

    def __apply__(self):
        """Set the result"""
        self.result = self.var.get()