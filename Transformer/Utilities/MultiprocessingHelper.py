# Transformer/Utilities/MultiprocessingHelper.py


# -------
# Imports
# -------

import multiprocessing;


# -------
# Classes
# -------

class Counter:
    def __init__(self, initialValue = 0, readLock = True):
        self._value = multiprocessing.Value('i', initialValue);

    def Current(self):
        with self._value.get_lock():
            return self._value.value;

    def Decrement(self, amount = 1):
        with self._value.get_lock():
            self._value.value -= amount;

    def Increment(self, amount = 1):
        with self._value.get_lock():
            self._value.value += amount;


# ---------
# Functions
# ---------

def CPUCount():
    cpuCount = 1;

    # According to the documentation, cpu_count() may raise a NotImplementedError; if this happens, issue a warning.

    try:
        cpuCount = multiprocessing.cpu_count();
    except NotImplementedError:
        warnings.warn("multiprocessing.cpu_count() is not implemented on this platform -> the CPU count will default to 1.", RuntimeWarning);

    return cpuCount;
