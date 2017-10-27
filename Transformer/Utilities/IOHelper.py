# Transformer/Utilities/IOHelper.py


# -------
# Imports
# -------

import os;
import shutil;
import sys;
import warnings;


# ---------
# Functions
# ---------

def ClearDirectory(directoryPath, removeSubdirectories = False):
    for entry in os.listdir(directoryPath):
        # Convert to an absolute path.

        absPath = os.path.join(directoryPath, entry);

        if os.path.isfile(absPath):
            os.remove(absPath);

        elif removeSubdirectories and os.path.isdir(absPath):
            # Only remove subdirectories if removeSubdirectories is set.

            shutil.rmtree(absPath);

def OpenForCSVWriter(filePath):
    # If running on Windows, try to open the file with newline = '' set (Python >= 3) to stop the csv module inserting extra blank lines
    # If this is not possible, issue a RuntimeWarning.

    if sys.platform.startswith("win"):
        if sys.version_info.major >= 3:
            return open(filePath, 'w', newline = '');
        else:
            warnings.warn("CSV files output from Python < 3 on Windows platforms may have blank lines between rows.", RuntimeWarning);

    return open(filePath, 'w');
