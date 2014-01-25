__author__ = 'Alex Zylstra'

from multiprocessing import Process, Pipe
from threading import Thread, Lock
import time
from impy.modules.Module import Module
from impy.implosions.Implosion import Implosion

def run(conn):
    """Function to handle running module in separate :py:class:`multithreading.Process`

    :param conn: A connection, i.e. one end of a `Pipe()`
    """
    # Need duck-checking instead of real type-checking...
    assert hasattr(conn, 'send') and hasattr(conn, 'recv')

    # Get the Module object from the pipe:
    mod = conn.recv()
    assert isinstance(mod, Module)
    # Get the Implosion object from the pipe:
    imp = conn.recv()
    assert isinstance(imp, Implosion)

    connLock = Lock()

    # Function for execution:
    def modRun():
        nonlocal imp,mod,conn, connLock
        try:
            mod.run(imp)
        except Exception as e:
            connLock.acquire()
            conn.send(e)
            connLock.release()

    t = Thread(target=modRun)
    t.start()

    while t.is_alive():
        connLock.acquire()
        conn.send(mod.progress())
        connLock.release()
        time.sleep(0.01)

    # When the thread is done, send the Module object back:
    conn.send(mod)

def ModuleRunner(app, imp, mod, dialog):
    """Take an implosion object and run the module run method in a separate process.

    :param app: A reference to the main application (:py:class:`impy`)
    :param imp: An :py:class:`Implosion` object
    :param mod: A :py:class:`Module` to run on the given implosion
    :param dialog: A dialog window (:py:class:`Progress_Dialog`)
    """
    assert isinstance(mod, Module)
    assert isinstance(imp, Implosion)

    parent_conn, child_conn = Pipe()

    p = Process(target=run, args=(child_conn,))
    app.processes.append(p)

    p.start()
    try:
        parent_conn.send(mod)
    except:
        raise Exception('Module object passed to ModuleRunner is not pickle-able!')

    try:
        parent_conn.send(imp)
    except:
        raise Exception('Implosion object passed to ModuleRunner is not pickle-able!')

    obj = None
    # Callpacks while the process is active:
    def callback():
        nonlocal p, dialog, parent_conn, app
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
                app.modControlVars[mod.name()].set(0)
                dialog.withdraw()
                p.terminate()
                return
            elif isinstance(obj, float) or isinstance(obj, int):
                dialog.set(100*obj)
            elif isinstance(obj, Module):
                # Pass back to the main app:
                app.modules[mod.name()].copy(obj)
                app.__postModule__(obj)
                dialog.withdraw()
                p.terminate()
                return  # Callback loop breaks here

        app.after(25, callback)  # loop

    # Start callback loop in a bit:
    app.after(10, callback)