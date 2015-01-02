# impy - a post-processor for HYADES implosion simulations
# Copyright (c) Massachusetts Institute of Technology / Alex Zylstra
# Distributed under the MIT License

import sys
from cx_Freeze import setup, Executable
import scipy
import scipy.interpolate
import platform

# GUI applications require a different base on Windows (the default is for a
# console application).
base = None
if sys.platform == "win32":
    base = "Win32GUI"
    
# Dependencies are automatically detected, but it might need fine tuning.
if platform.system() == 'Darwin':
    options = {"packages": ["os","scipy","scipy.interpolate","impy","impy.gui","impy.implosions","impy.modules","impy.resources","impy.util","impy.gui.widgets"],
                                            "includes": ["numpy","scipy","scipy.interpolate","scipy.optimize","scipy.stats","matplotlib","matplotlib.pyplot","matplotlib.backends.backend_macosx","matplotlib.backends.backend_tkagg"],
                                            "excludes": [],
                                            "include_files": [('impy/implosions/LILAC_Materials.csv', 'LILAC_Materials.csv'),
                                                                ('/usr/local/lib/python3.3/site-packages/scipy/sparse/sparsetools/_csr.so','_csr.so'),
                                                                ('/usr/local/lib/python3.3/site-packages/scipy/sparse/sparsetools/_csc.so','_csc.so'),
                                                                ('/usr/local/lib/python3.3/site-packages/scipy/sparse/sparsetools/_coo.so','_coo.so'),
                                                                ('/usr/local/lib/python3.3/site-packages/scipy/sparse/sparsetools/_dia.so','_dia.so'),
                                                                ('/usr/local/lib/python3.3/site-packages/scipy/sparse/sparsetools/_bsr.so','_bsr.so'),
                                                                ('/usr/local/lib/python3.3/site-packages/scipy/sparse/sparsetools/_csgraph.so','_csgraph.so')],
                                            "optimize": 2
                                            }
    scripts = []
    
if platform.system() == 'Windows':
    options = {"packages": ["os","scipy","numpy","scipy.interpolate","impy","impy.gui","impy.implosions","impy.modules","impy.resources","impy.util","impy.gui.widgets"],
                                            "includes": ["numpy","scipy","scipy.interpolate","scipy.optimize","scipy.stats","matplotlib","matplotlib.pyplot","matplotlib.backends.backend_macosx","matplotlib.backends.backend_tkagg"],
                                            "excludes": [],
                                            "include_files": [('logo/logo.ico','favicon.ico'),
                                                                ('impy/implosions/LILAC_Materials.csv', 'LILAC_Materials.csv'),
                                                                ('C:\Python34\Lib\site-packages\scipy\sparse\_csparsetools.pyd','_csparsetools.pyd'),
                                                                ('C:\Python34\Lib\site-packages\scipy\sparse\_sparsetools.pyd','_sparsetools.pyd'),
                                                                ('C:\Python34\Lib\site-packages\scipy\sparse\csgraph\_min_spanning_tree.pyd','_min_spanning_tree.pyd'),
                                                                ('C:\Python34\Lib\site-packages\scipy\sparse\csgraph\_shortest_path.pyd','_shortest_path.pyd'),
                                                                ('C:\Python34\Lib\site-packages\scipy\sparse\csgraph\_tools.pyd','_tools.pyd'),
                                                                ('C:\Python34\Lib\site-packages\scipy\sparse\csgraph\_traversal.pyd','_traversal.pyd')
                                                                ],
                                            "optimize": 2,
                                             "include_msvcr": True,
                                             "icon": "logo/logo.ico"
                                            }
    scripts = ['win_postinst.py']

mac_options = {"iconfile": "logo/logo.icns"}
# http://msdn.microsoft.com/en-us/library/windows/desktop/aa371847(v=vs.85).aspx
shortcut_table = [
    ("DesktopShortcut",        # Shortcut
     "DesktopFolder",          # Directory_
     "impy",                   # Name
     "TARGETDIR",              # Component_
     "[TARGETDIR]impy.exe",    # Target
     None,                     # Arguments
     None,                     # Description
     None,                     # Hotkey
     None, # Icon
     None,                     # IconIndex
     None,                     # ShowCmd
     'TARGETDIR'               # WkDir
     )
    ]

# Now create the table dictionary
msi_data = {"Shortcut": shortcut_table}

# Change some default MSI options and specify the use of the above defined tables
msi_options = {'data': msi_data}


setup(  name = "impy",
        version = "0.2",
        description = "Implosion simulation post-processor",
    	packages=['impy', 'impy.gui', 'impy.gui.widgets', 'impy.implosions', 'impy.modules', 'impy.resources', 'impy.modules'],
        options = {"build_exe": options, "bdist_mac": mac_options, "bdist_msi": msi_options},
        scripts = scripts,
        executables = [Executable("impy.py", base=base, copyDependentFiles=True, compress=True)])
