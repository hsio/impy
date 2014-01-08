__author__ = 'Alex Zylstra'
__date__ = '2014-01-06'
__version__ = '1.0.0'

# adapted from: http://stackoverflow.com/questions/13141259/expandable-collapsible-frame-in-python-tkinter

import tkinter as tk
import tkinter.ttk as ttk

class Collapsible_Frame(tk.Frame):
    """Implement a Tkinter frame which can be hidden/shown by the user.
    This is an extension of the built-in tkinter Frame.

    :param parent: The parent (i.e. containing) object.
    :param text: (optional) Text to display at the header.

    Any additional kwargs are passed directly to :py:class:`tkinter.Frame`'s constructor.

    Any elements that are intended to be added to the hideable area should be added to the `subFrame` member::

        c = Collapsible_Frame(self)
        ttk.Label(c, 'Hello World!')

    will create a label within the subframe. The subframe can be controlled as you wish.

    :author: Alex Zylstra
    :date: 2014-01-06
    """

    def __init__(self, parent, text='', **options):
        tk.Frame.__init__(self, parent, **options)

        self.configure(background='#eeeeee')

        self.show=tk.IntVar()
        self.show.set(0)
        self.titleFrame=ttk.Frame(self)
        self.titleFrame.pack(fill=tk.X, expand=1)
        ttk.Label(self.titleFrame, text=text).pack(side=tk.LEFT, fill=tk.X, expand=1)
        self.toggleButton=ttk.Checkbutton(self.titleFrame, width=2,text='+', command=self.toggle, variable=self.show, style='Toolbutton')
        self.toggleButton.pack(side=tk.LEFT)
        self.subFrame=tk.Frame(self, relief=tk.SUNKEN, borderwidth=1)

        # stretch the column to fill all space:
        tk.Grid.columnconfigure(self.subFrame, 0, weight=1)
        tk.Grid.columnconfigure(self.subFrame, 1, weight=1)
        tk.Grid.columnconfigure(self.subFrame, 2, weight=1)
        tk.Grid.columnconfigure(self.subFrame, 3, weight=1)

    def toggle(self):
        """Toggle the state of the frame, i.e. toggle between displayed and hidden."""
        if bool(self.show.get()):
            self.toggle_visible()
        else:
            self.toggle_hidden()

    def toggle_visible(self):
        """Set the frame to a visible state."""
        self.subFrame.pack(fill=tk.X, expand=1)
        self.toggleButton.configure(text='-')

    def toggle_hidden(self):
        """Set the frame to a hidden state."""
        self.subFrame.forget()
        self.toggleButton.configure(text='+')