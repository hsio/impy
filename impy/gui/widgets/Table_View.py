__author__ = 'Alex Zylstra'

import tkinter as tk
import tkinter.font as tkFont
import tkinter.ttk as ttk

def __sortby__(tree, col, descending):
    """Sort tree contents when a column is clicked on."""
    # grab values to sort
    data = [(tree.set(child, col), child) for child in tree.get_children('')]

    # reorder data
    data.sort(reverse=descending)
    for indx, item in enumerate(data):
        tree.move(item[1], '', indx)

    # switch the heading so that it will sort in the opposite direction
    tree.heading(col,
        command=lambda col=col: __sortby__(tree, col, int(not descending)))

class Table_Viewer_Frame(ttk.Frame):
    """Implement a tabular display inside a :py:class`tkinter.Frame`.

    :param columns: (optional) A :py:class:`tuple` containing :py:class:`str` names of each column
    :param data: (optional) The data to display. Should be a :py:class:`list` where each element is a :py:class:`tuple`
    corresponding to column values. The length of each tuple must match the length of the `columns` argument.
    :param parent: (optional) The parent of this frame [default=None]
    :param build: (optional) Whether to call the widget building functions in the constructor [default=True]

    :author: Alex Zylstra
    :date: 2014-01-11
    :version: 1.0.0
    """
    def __init__(self, columns=("Quantity", "Value"), data=[("","",)], parent=None, build=True):
        """Constructor"""
        super(Table_Viewer_Frame, self).__init__(parent)

        # initializations:
        self.tree = None
        self.tree_columns = columns  # the column headings
        self.tree_data = data # the tree data

        self.__setup_widgets__()

        if build:
            self.__build_tree__()

    def __setup_widgets__(self):
        self.pack(fill='both', expand=True)

        # a treeview with scrollbars.
        self.tree = ttk.Treeview(self, columns=self.tree_columns, show="headings")
        vsb = ttk.Scrollbar(self, orient="vertical", command=self.tree.yview)
        hsb = ttk.Scrollbar(self, orient="horizontal", command=self.tree.xview)
        self.tree.configure(yscrollcommand=vsb.set, xscrollcommand=hsb.set)
        self.tree.grid(column=0, row=0, sticky='nsew', in_=self)
        vsb.grid(column=1, row=0, sticky='ns', in_=self)
        hsb.grid(column=0, row=1, sticky='ew', in_=self)

        self.grid_columnconfigure(0, weight=1)
        self.grid_rowconfigure(0, weight=1)

    def __build_tree__(self):
        """Build the tree part of the widget"""
        # start by clearing everything already in the tree:
        for iid in self.tree.get_children():
            self.tree.delete(iid)

        # add in new data:
        for col in self.tree_columns:
            self.tree.heading(col, text=col,
                command=lambda c=col: __sortby__(self.tree, c, 0))
            # set up with width based on font
            self.tree.column(col, width=tkFont.Font().measure(col.title()))

        for item in self.tree_data:
            self.tree.insert('', 'end', values=item)

            # adjust columns lengths if necessary
            for indx, val in enumerate(item):
                ilen = tkFont.Font().measure(val)
                if self.tree.column(self.tree_columns[indx], width=None) < ilen:
                    self.tree.column(self.tree_columns[indx], width=ilen)

    def setData(self, columns, data):
        """Set the data to be displayed.

        :param columns: (optional) A :py:class:`tuple` containing :py:class:`str` names of each column
        :param data: (optional) The data to display. Should be a :py:class:`list` where each element is a :py:class:`tuple`
        corresponding to column values. The length of each tuple must match the length of the `columns` argument.
        """
        self.tree_columns = ("Quantity", "Value")  # the column headings
        self.tree_data = [("",      "",)]  # the tree data

        # Update the display:
        self.refresh()

    def refresh(self, *args):
        """Refresh the display"""
        self.__build_tree__()

class Table_Viewer(tk.Toplevel):
    """Implement a top-level window to display info in a tabular fashion.

    :param columns: (optional) A :py:class:`tuple` containing :py:class:`str` names of each column
    :param data: (optional) The data to display. Should be a :py:class:`list` where each element is a :py:class:`tuple`
    corresponding to column values. The length of each tuple must match the length of the `columns` argument.
    :param parent: (optional) The parent of this window [default=None]
    :param build: (optional) Whether to call the widget building functions in the constructor [default=True]
    :param title: (optional) The window title [default='']
    :param widgets: (optional) A list of tkinter :py:class`tkinter.Widget`s to be added to this window [default=None]

    :author: Alex Zylstra
    :date: 2014-01-10
    :version: 1.0.0
    """


    def __init__(self, columns=("Quantity", "Value"), data=[("","",)], parent=None, build=True, title='', widgets=[]):
        """Initialize the table."""
        super(Table_Viewer, self).__init__(parent)
        self.header_widgets = []  # control widgets to display at the top of the window

        self.protocol("WM_DELETE_WINDOW", self.close)

        self.configure(background='#eeeeee')

        # add header widgets to the GUI:
        for widget in widgets:
            widget.pack()

        self.frame = Table_Viewer_Frame(columns=columns, data=data, parent=self, build=build)

        self.title(title)

        # a couple key bindings:
        self.bind('<Escape>', self.close)
        self.bind('<Control-r>', self.frame.refresh)
        self.bind('<Command-r>', self.frame.refresh)

    def close(self, *args):
        """Close the window"""
        self.withdraw()