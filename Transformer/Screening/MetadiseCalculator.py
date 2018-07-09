# Transformer/Screening/MetadiseCalculator.py


# -------
# Imports
# -------

import os;
import re;
import shutil;
import subprocess;
import warnings;

from Transformer.IO import StructureIO;
from Transformer.Screening.TotalEnergyCalculatorBase import TotalEnergyCalculatorBase;
from Transformer.Utilities import IOHelper;


# ------------------------
# MetadiseCalculator Class
# ------------------------

class MetadiseCalculator(TotalEnergyCalculatorBase):
    # -----------
    # Constructor
    # -----------

    def __init__(
        self,
        potenPath = None, speciesCharges = None,
        metadiseExe = None, tempDirectory = None,
        useMP = False, mpNumProcesses = None,
        atomicSymbolLookupTable = None, removeAtoms = None
        ):

        # Set default parameters.

        if metadiseExe == None:
            metadiseExe = MetadiseCalculator.DefaultMetadiseExe;

        # Only set potenPath if species charges are not supplied.

        if speciesCharges == None and potenPath == None:
            potenPath = MetadiseCalculator.DefaultPotenPath;

        # Parameter checking.

        if potenPath != None:
            if not os.path.isfile(potenPath):
                raise Exception("Error: The supplied potenPath \"{0}\" does not exist.".format(potenPath));

            # If both a poten file and a set of charges have been supplied, the poten file takes precedence -> print a warning.

            if speciesCharges != None:
                warnings.warn("If both a potential file and species charges are supplied, the former takes precedence and the potential will be read from file.", UserWarning);
        else:
            if speciesCharges == None:
                raise Exception("Error: One of potenPath or speciesCharges must be supplied.");

        # Call the base class constructor.

        super(MetadiseCalculator, self).__init__(
            tempDirectory = tempDirectory,
            useMP = useMP, mpNumProcesses = mpNumProcesses,
            atomicSymbolLookupTable = atomicSymbolLookupTable, removeAtoms = removeAtoms
            );

        # Store fields.

        self._potenPath = potenPath;
        self._speciesCharges = speciesCharges;

        self._metadiseExe = metadiseExe;

    # --------------
    # Pubilc Methods
    # --------------

    def CalculatorNumThreads(self):
        # Metadise is presently a serial code.

        return 1;

    def CalculatorRequiresTempDir(self):
        # Metadise does require a temporary directory.

        return True;

    def GetTotalEnergy(self, structure, degeneracy, raiseOnError, tempDirectory):
        atomicSymbolLookupTable = self._atomicSymbolLookupTable;

        # Copy/generate a poten.txt file in the temporary directory.

        if self._potenPath != None:
            shutil.copy(
                self._potenPath, os.path.join(tempDirectory, r"poten.txt")
                );
        else:
            speciesCharges = self._speciesCharges;

            # Get a list of atomic symbols for the atom types in the structure, and check there is an entry for each in the supplied species parameters.

            atomicSymbols, _ = structure.GetAtomicSymbolsCounts(atomicSymbolLookupTable = atomicSymbolLookupTable);

            for atomicSymbol in atomicSymbols:
                if atomicSymbol not in speciesCharges:
                    raise Exception("Error: Atom type '{0}' not defined in the supplied species parameters.".format(atomicSymbol));

            # Write out a poten file listing the supplied species charges and masses.

            with open(os.path.join(tempDirectory, "poten.txt"), 'w') as outputWriter:
                outputWriter.write("potential\n");

                outputWriter.write("    species\n");

                for symbol in atomicSymbols:
                    charge = speciesCharges[symbol];
                    outputWriter.write("        {0: >4}  core  {1: >8.3f}\n".format(symbol, charge));

                outputWriter.write("    ends\n");

                outputWriter.write("ends\n");

        # Write the structure as a VASP POSCAR file.

        StructureIO.WriteStructure(
            structure, os.path.join(tempDirectory, r"POSCAR.vasp"), atomicSymbolLookupTable = atomicSymbolLookupTable
            );

        # Write out a basic Metadise input file.

        with open(os.path.join(tempDirectory, r"input.txt"), 'w') as outputWriter:
            # Include block to read structure from a POSCAR-format file.

            outputWriter.write("vasp\n");
            outputWriter.write("    @include POSCAR.vasp\n");
            outputWriter.write("ends\n");
            outputWriter.write("\n");

            # Include block to read the potential from "poten.txt".

            outputWriter.write("@include poten.txt\n");
            outputWriter.write("\n");

            # Specify calculation (single-point energy/force).

            outputWriter.write("force\n");
            outputWriter.write("stop 99\n");
            outputWriter.write("\n");

            # Start/finish processing.

            outputWriter.write("start\n");
            outputWriter.write("stop\n");

        # Enter the temporary directory and call Metadise.

        workingDirectory = os.getcwd();

        os.chdir(tempDirectory);

        returnCode = subprocess.call(self._metadiseExe);

        os.chdir(workingDirectory);

        # Check the return code and try to extract a total energy from the summary file.

        totalEnergy = None;

        if returnCode == 0:
            summaryPath = os.path.join(tempDirectory, "summ_o0001.out");

            if os.path.isfile(summaryPath):
                with open(summaryPath, 'r') as inputReader:
                    for line in inputReader:
                        match = MetadiseCalculator._SummaryTotalEnregyRegex.search(line);

                        if match:
                            totalEnergy = float(match.group('energy'));

        # If raiseOnError is set, check whether we were able to extract a total energy.

        if totalEnergy == None and raiseOnError:
            raise Exception("Error: Unable to calculate a total energy - please check input/output in \"{0}\".".format(tempDirectory));

        # Clear the temporary directory.

        IOHelper.ClearDirectory(tempDirectory);

        # Return the total energy.

        return totalEnergy;

    # -------------
    # Static Fields
    # -------------

    _SummaryTotalEnregyRegex = re.compile(r"latt en=\s*(?P<energy>[+-]?\d+\.\d+)");

    DefaultPotenPath = r"./poten.txt";
    DefaultMetadiseExe = r"metadise";
