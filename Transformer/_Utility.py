# Transformer/_Utility.py by J. M. Skelton


# -------
# Imports
# -------

import math;
import multiprocessing;

import numpy as np;


# ---------
# Functions
# ---------

def GetCPUCount():
    cpuCount = 1;

    # According to the documentation, cpu_count() may raise a NotImplementedError; if this happens, issue a warning.

    try:
        cpuCount = multiprocessing.cpu_count();
    except NotImplementedError:
        warnings.warn("multiprocessing.cpu_count() is not implemented on this platform -> the CPU count will default to 1.", RuntimeWarning);

    return cpuCount;
