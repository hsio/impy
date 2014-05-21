
import tkinter as tk
import tkinter.ttk as ttk
import platform

class Progress_Dialog(tk.Toplevel):
    """Implement a generic progress bar dialog. A cancel event can be retrieved by the
    `cancelled` member. By default the range of the progress bar corresponds to 0-100.

    :param parent: (optional) the tkinter parent (usually should be None) [default=None]
    :param maximum: (optional) The maximum value, corresponding to the bar completely filled [default=100]
    """
    __author__ = 'Alex Zylstra'
    __date__ = '2014-01-23'
    __version__ = '0.1.0'

    def __init__(self, parent=None, maximum=100):
        """Initialize the progress dialog"""
        super(Progress_Dialog, self).__init__(parent)
        self.cancelled = False
        self.__createUI__(maximum)

        # Set window background
        if platform.system() == 'Darwin':
            self.configure(background='#E8E9E8')
        else:
            self.configure(background='#F1F1F1')

    def __createUI__(self, maximum):
        """Helper method to create the UI elements"""
        self.grid()
        self.label = ttk.Label(self, text="Importing WRF")
        self.label.grid(sticky='N', padx=2, pady=2)

        self.counter = tk.IntVar()
        self.progress_bar = ttk.Progressbar(self, variable=self.counter, maximum=maximum)
        self.progress_bar.grid(sticky='N', padx=2, pady=2)

        self.cancel_button = ttk.Button(self, text="Cancel", command=self.cancel)
        self.cancel_button.grid(sticky='S', padx=2, pady=2)

    def step(self, amount):
        """Increase the amount on the progress bar.

        :param amount: Incremental step
        """
        self.progress_bar.step(amount)

    def set(self, amount):
        """Set the value of the progress bar.

        :param amount: The value to set the progress bar to.
        """
        self.counter.set(amount)

    def cancel(self):
        """Attempt to cancel the operation."""
        self.cancelled = True

    def set_text(self, new_text):
        """Update the displayed info at the top of the window

        :param new_text: The new text string to display
        """
        self.label.configure(text=new_text)