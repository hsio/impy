# Global definitions for IO
# A. Zylstra 2012/05/24
import os

OutputDir = "Outputs"
#create the directory if it doesn't exist
d = os.path.dirname(OutputDir)
if not os.path.exists(OutputDir):
    os.makedirs(OutputDir)