__author__ = 'Alex Zylstra'
__date__ = '2014-03-26'
__version__ = '1.0.0'

import tkinter as tk
import numpy as np
import scipy, scipy.signal
import re

def parseGeometry(geometry):
    """Parse a tkinter geometry string.

    :param geometry: A geometry string from `Toplevel.geometry`
    :returns: A list containing [width, height, x, y]
    """
    m = re.match("(\d+)x(\d+)([-+]\d+)([-+]\d+)", geometry)
    if not m:
        raise ValueError("failed to parse geometry string")
    return list(map(int, m.groups()))

class WindowManager:
    """Window management functionality. This class acts as a manager for placement of all windows
    created by modules and implosions via one of a few methods.

    :param screenWidth: The width of the screen in pixels
    :param screenHeight: The height of the screen in pixels
    :author: Alex Zylstra
    :date: 2014-03-26
    """
    #: Step size in horizontal placement [pixels]
    dx = 10
    #: Step size in vertical placement [pixels]
    dy = 10

    def __init__(self, screenWidth, screenHeight):
        """Constructor for a new window manager."""
        self.windows = []
        self.screenWidth = screenWidth
        self.screenHeight = screenHeight

    def getLocation(self, width, height):
        """Get a location for a new window of specified size. Units in pixels.

        :param width: The width of the window
        :param height: The height of the window
        :returns: tuple containing coordinates (x,y)
        """
        usageMatrix = np.ones((self.screenHeight/20,self.screenWidth/20), dtype=np.float32)

        for w in self.windows:
            # Wrap in try/except block in case windows are totally deleted
            try:
                assert isinstance(w, tk.Toplevel) or isinstance(w, tk.Tk)
                if w.wm_state() == 'normal':
                    g = parseGeometry(w.geometry())
                    # width of extra window stuff including title bar:
                    borderx = w.winfo_rootx() - g[2]
                    bordery = w.winfo_rooty() - g[3]

                    x = np.floor(g[2]/20)
                    y = np.floor((g[3])/20)
                    w = np.ceil((g[0]+borderx)/20)
                    h = np.ceil((g[1]+bordery)/20)
                    usageMatrix[y:y+h,x:x+w] = 0
            except tk.TclError:
                self.windows.remove(w)
            except Exception as e:
                print('An exception occurred in WindowManager.getLocation: ', e)
                pass
            
        b = np.ones((height/20,width/20), dtype=np.float32)
        loc_ok = (scipy.signal.convolve2d(usageMatrix, b, 'valid') == int(width/20)*int(height/20))
        try:
            new_loc = np.argwhere(loc_ok)[0]
            return 20*new_loc[1], 20*new_loc[0]
        except Exception as e:
            print('An exception occurred in WindowManager.getLocation: ', e)
            return 0,0

    def addWindow(self, w, passive=False):
        """Add a new window.

        :param w: The new window (must be a :py:class:`tkinter.Toplevel` object)
        :param passive: (optional) if set to True, the `WindowManager` will not change the current placement of w.
        """
        assert isinstance(w, tk.Toplevel) or isinstance(w, tk.Tk)

        w.update_idletasks()

        # if requested, place the window:
        if not passive:
            # Get the width and height, then use the other function in this class to get a spot to put it:
            g = parseGeometry(w.geometry())
            width = g[0]
            height = g[1]
            x,y = self.getLocation(width, height)
            # place:
            w.geometry("%dx%d%+d%+d" % (width, height, x, y))
            w.update_idletasks()

        # add to the internal list of managed windows:
        self.windows.append(w)
