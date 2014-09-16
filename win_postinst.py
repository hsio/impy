# impy - a post-processor for HYADES implosion simulations
# Copyright (c) Massachusetts Institute of Technology / Alex Zylstra
# Distributed under the MIT License

import os
import sys
import shutil
import impy

DESKTOP_FOLDER = get_special_folder_path("CSIDL_DESKTOPDIRECTORY")
APP_FOLDER = os.path.join(sys.prefix,'Program Files','impy')
NAME = 'impy.lnk'

if sys.argv[1] == '-install':
    create_shortcut(
        os.path.join(APP_FOLDER, 'impy.exe'), # program
        'Implosion post-processor', # description
        NAME, # filename
        impy.__file__, # parameters
        APP_FOLDER, # workdir
        os.path.join(APP_FOLDER, 'favicon.ico'), # iconpath
    )
    # move shortcut from current directory to DESKTOP_FOLDER
    shutil.move(os.path.join(os.getcwd(), NAME),
                os.path.join(DESKTOP_FOLDER, NAME))
    # tell windows installer that we created another
    # file which should be deleted on uninstallation
    file_created(os.path.join(DESKTOP_FOLDER, NAME))

if sys.argv[1] == '-remove':
    pass
    # This will be run on uninstallation. Nothing to do.
