__author__ = 'Alex Zylstra'

from multiprocessing import Process, Pipe
from threading import Thread, Lock
import time
from impy.implosions.Implosion import Implosion

def run(conn):
    """Function to handle running implosion generation in separate :py:class:`multithreading.Process`

    :param conn: A connection, i.e. one end of a `Pipe()`
    """
    # Need duck-checking instead of real type-checking...
    assert hasattr(conn, 'send') and hasattr(conn, 'recv')

    # Get the implosion object from the pipe:
    imp = conn.recv()
    assert isinstance(imp, Implosion)

    connLock = Lock()

    # Run in a separate thread in this process:
    def impRun():
        nonlocal imp, conn
        try:
            imp.generate()
        except Exception as e:
            connLock.acquire()
            conn.send(e)
            connLock.release()

    t = Thread(target=impRun)
    t.start()

    while t.is_alive():
        connLock.acquire()
        conn.send(imp.progress())
        connLock.release()
        time.sleep(0.01)

    # When the thread is done, send the Implosion object back:
    conn.send(imp)

def ImplosionRunner(app, imp, dialog):
    """Take an implosion object and run the generation method in a separate process.

    :param app: A reference to the main application (:py:class:`impy`)
    :param imp: An :py:class:`Implosion` object
    :param dialog: A dialog window (:py:class:`Progress_Dialog`)
    """
    assert isinstance(imp, Implosion)

    parent_conn, child_conn = Pipe()

    p = Process(target=run, args=(child_conn,))
    app.processes.append(p)

    p.start()
    try:
        parent_conn.send(imp)
    except:
        raise Exception('Implosion object passed to ImplosionRunner is not pickleable!')

    obj = None
    # Loop while the process is active:
    def callback():
        nonlocal dialog, p, parent_conn, app
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
                showerror('Error!', 'A problem occurred generating the implosion (class '+imp.name()+')\n'+obj.__str__())
                dialog.withdraw()
                p.terminate()
                return
            elif isinstance(obj, float):
                dialog.set(100*obj)
            elif isinstance(obj, Implosion):
                # Pass info back to the main app:
                app.imp = obj
                app.__postImplosion__()
                dialog.withdraw()
                p.terminate()
                return

        app.after(25, callback)

    app.after(10, callback)