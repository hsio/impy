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