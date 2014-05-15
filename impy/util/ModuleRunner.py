__author__ = 'Alex Zylstra'

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
