# Transformer/Screening/TotalEnergyCalculatorBase.py


# -------
# Imports
# -------

import multiprocessing;
import os;
import shutil;
import time;

try:
    # Python 2.x.

    from Queue import Empty;
except ImportError:
    # Python >= 3.

    from queue import Empty;

from Transformer.Utilities import IOHelper;
from Transformer.Utilities import MultiprocessingHelper;

# Try to import the tqdm module to display progress bars.

_TQDM = False;

try:
    import tqdm;

    _TQDM = True;
except ImportError:
    pass;


# -------------------------------
# TotalEnergyCalculatorBase Class
# -------------------------------

class TotalEnergyCalculatorBase(object):
    # -----------
    # Constructor
    # -----------

    def __init__(self, tempDirectory = None, useMP = False, mpNumProcesses = None):
        # If tempDir is not set, set it to the default value.
        # This includes the PID of the current process to (hopefully!) make it unique.

        if tempDirectory == None:
            tempDirectory = TotalEnergyCalculatorBase.DefaultTempDirectory.replace("<PID>", str(os.getpid()));

        # If useMP is set and mpNumProcesses is not supplied, set it to the CPU count.

        if useMP and mpNumProcesses == None:
            mpNumProcesses = mpNumProcesses = MultiprocessingHelper.CPUCount();

        # Store fields.

        self._tempDirectory = tempDirectory;

        self._useMP = useMP;
        self._mpNumProcesses = mpNumProcesses;

    # ---------------
    # Private Methods
    # ---------------

    def _CreateTempDirectories(self, numWorkerProcesses = None):
        tempDirectories = None;

        # Build and create a list of temporary directories, if required.

        if self.CalculatorRequiresTempDir():
            tempDirectoryBase = self._tempDirectory;

            if numWorkerProcesses == None:
                # If there are no worker processes, use the base directory.

                tempDirectories = [tempDirectoryBase];
            else:
                # If we are using multiple worker processes, use the base directory with an appended process number.

                tempDirectories = [
                    r"{0}.{1}".format(tempDirectoryBase, i + 1)
                        for i in range(0, numWorkerProcesses)
                    ];

            # Create directories.

            for tempDirectory in tempDirectories:
                if not os.path.isdir(tempDirectory):
                    os.makedirs(tempDirectory);
                else:
                    if len(os.listdir(tempDirectory)) != 0:
                        raise Exception("Error: Temporary directory \"{0}\" is not empty.".format(tempDirectory));

        # Return directories.

        if tempDirectories == None:
            return None;
        elif numWorkerProcesses == None:
            return tempDirectories[0];
        else:
            return tempDirectories;

    def _GetTotalEnergy(self, structure, degeneracy, raiseOnError, tempDirectory):
        # Pass to the GetTotalEnergy() method.

        totalEnergy = self.GetTotalEnergy(structure, degeneracy, raiseOnError = raiseOnError, tempDirectory = tempDirectory);

        # If using a temporary directory, clear it.

        if tempDirectory != None:
            IOHelper.ClearDirectory(tempDirectory, removeSubdirectories = True);

        # Return the calculated total energy.

        return totalEnergy;

    def _CalculateTotalEnergies(self, structures, degeneracies, raiseOnError, progressBar):
        # Create a temporary directory if required.

        tempDirectory = self._CreateTempDirectories();

        # Set up a primary iterator.

        iValues = range(0, len(structures));

        # If progressBar is set and the tqdm module is available, wrap the iterator in a tqdm progress bar.

        progressBar = progressBar and _TQDM;

        if progressBar:
            iValues = tqdm.tqdm(iValues);

        # Calculate a list of total energies.

        result = [
            self._GetTotalEnergy(structures[i], degeneracies[i], raiseOnError, tempDirectory)
                for i in iValues
            ];

        # If progressBar is set and the tqdm module was imported, print a blank line after the progress bar.

        if progressBar:
            print("");

        # Remove the temporary directory.

        if tempDirectory != None:
            shutil.rmtree(tempDirectory);

        # Return the result.

        return result;

    def _CalculateTotalEnergiesMP(self, structures, degeneracies, raiseOnError, progressBar):
        # Work out how many worker processes to run.

        numWorkerProcesses = min(
            self._mpNumProcesses // self.CalculatorNumThreads(), len(structures)
            );

        # Get a base temporary directory if required.

        tempDirectories = self._CreateTempDirectories(numWorkerProcesses = numWorkerProcesses);

        # Initialise Mapper objects for use with the QueueMap() routine from the MultiprocessingHelper module.

        mappers = [
            _CalculateTotalEnergiesMP_Mapper(self, raiseOnError, tempDirectory)
                for tempDirectory in tempDirectories
            ];

        # Calculate total energies.

        totalEnergies = MultiprocessingHelper.QueueMap(
            [item for item in zip(structures, degeneracies)], mappers, progressBar = progressBar
            );

        # Remove the temporary directories.

        if tempDirectories != None:
            for tempDirectory in tempDirectories:
                shutil.rmtree(tempDirectory);

        # Return the calculated total energies.

        return totalEnergies;

    # --------------
    # Pubilc Methods
    # --------------

    def CalculatorNumThreads(self):
        raise NotImplementedError("Error: CalculatorNumThreads() must be implemented in derived classes.");

    def CalulatorRequiresTempDir(self):
        raise NotImplementedError("Error: CalculatorRequiresTempDir() must be implemented in derived classes.");

    def GetTotalEnergy(self, structure, degeneracy, raiseOnError, tempDirectory):
        raise NotImplementedError("Error: GetTotalEnergy() must be implemented in derived classes.");

    def CalculateTotalEnergy(self, structure, degeneracy = 1, raiseOnError = True, tempDirectory = None):
        # Paraeter validation.

        if structure == None:
            raise Exception("Error: structure cannot be None.");

        # Create a temporary directory if required.

        tempDirectory = self._CreateTempDirectories();

        # Get the total energy.

        totalEnergy = self.GetTotalEnergy(structure, degeneracy, raiseOnError, tempDirectory);

        # Delete the temporary directory.

        if tempDirectory != None:
            shutil.rmtree(tempDirectory);

        # Return the total energy.

        return totalEnergy;

    def CalculateTotalEnergies(self, structures, degeneracies = None, raiseOnError = True, progressBar = True):
        # Parameter validation.

        if structures == None:
            raise Exception("Error: structures cannot be None.");

        if degeneracies != None and len(degeneracies) != len(structures):
            raise Exception("Error: If supplied, degeneracies must have the same length as structures.");

        # If degeneracies is not supplied, initialise it to a list of ones.

        if degeneracies == None:
            degeneracies = [1] * len(structures);

        # If the _useMP field is set and _mpNumProcesses is at least twice the value of CalculatorNumThreads(), calculate the total energies using the internal parallelisation; if not, fall back to the serial routine.

        totalEnergies = None;

        if self._useMP and (self._mpNumProcesses // self.CalculatorNumThreads() >= 2):
            totalEnergies = self._CalculateTotalEnergiesMP(structures, degeneracies, raiseOnError, progressBar);
        else:
            totalEnergies = self._CalculateTotalEnergies(structures, degeneracies, raiseOnError, progressBar);

        # Return the list of total energies.

        return totalEnergies;

    def CalculateTotalEnergiesGrouped(self, spacegroupGroups, raiseOnError = True, progressBar = True):
        # Parameter validation.

        if spacegroupGroups == None:
            raise Exception("Error: spacegroupGroups cannot be None.");

        # Merge the grouped structures and degeneracies into flat lists.

        spacegroups = [
            key for key in spacegroupGroups.keys()
            ];

        structuresFlat, degeneraciesFlat = [], [];

        for key in spacegroups:
            structures, degeneracies = spacegroupGroups[key];

            structuresFlat = structuresFlat + structures;
            degeneraciesFlat = degeneraciesFlat + degeneracies;

        # Pass into the CalculateTotalEnergies().

        totalEnergies = self.CalculateTotalEnergies(structuresFlat, degeneraciesFlat, raiseOnError = raiseOnError, progressBar = progressBar);

        # Group the total energies by spacegroup.

        totalEnergyGroups = { };

        pointer = 0;

        for key in spacegroups:
            structures, _ = spacegroupGroups[key];
            structureCount = len(structures);

            totalEnergyGroups[key] = totalEnergies[pointer:pointer + structureCount];
            pointer += structureCount;

        # Return the grouped total energies.

        return totalEnergyGroups;

    # -------------
    # Static Fields
    # -------------

    _CalculateTotalEnergiesMP_PollDelay = 0.1;
    _CalculateTotalEnergiesMP_QueueBurstFactor = 100;

    DefaultTempDirectory = "_CalculatorTemp.<PID>";


# --------------------------------------
# _CalculateTotalEnergiesMP_Mapper Class
# --------------------------------------

class _CalculateTotalEnergiesMP_Mapper(MultiprocessingHelper.MapperBase):
    def __init__(self, totalEnergyCalculator, raiseOnError, tempDirectory):
        self._totalEnergyCalculator = totalEnergyCalculator;

        self._raiseOnError = raiseOnError;
        self._tempDirectory = tempDirectory;

    def Map(self, arg):
        structure, degeneracy = arg;

        return self._totalEnergyCalculator._GetTotalEnergy(
            structure, degeneracy, raiseOnError = self._raiseOnError, tempDirectory = self._tempDirectory
            );
