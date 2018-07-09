# Transformer/Framework/BatchIO.py


# -------
# Imports
# -------

import os;
import warnings;

from Transformer.IO import StructureSetIO;


# -----------------------------
# Batch Import/Export Functions
# -----------------------------

def ExportResultSet(resultSet, fileFormat = 'vasp', prefix = None, archiveDirectory = "./", atomicSymbolLookupTable = None, **kwargs):
    # Set up a generator.

    generator = ExportResultSetPassthrough(
        resultSet, fileFormat = fileFormat, prefix = prefix, archiveDirectory = archiveDirectory, atomicSymbolLookupTable = atomicSymbolLookupTable, **kwargs
        );

    # Run the generator to output the structure sets, but do nothing with the yield values.

    for _ in generator:
        pass;

def ExportResultSetPassthrough(resultSetGenerator, fileFormat = 'vasp', prefix = None, archiveDirectory = "./", atomicSymbolLookupTable = None, **kwargs):
    # prefix should not contain underscores; if it does, issue a warning and replace them with hyphens.

    if prefix != None:
        if '_' in prefix:
            warnings.warn("Underscores in prefix will be converted to hyphens.", UserWarning);

            prefix = prefix.replace('_', '-');

    # If archiveDirectory does not exist, create it.

    if not os.path.isdir(archiveDirectory):
        os.makedirs(archiveDirectory);

    # Loop over structure sets in the result set.

    for i, structureSet in enumerate(resultSetGenerator):
        # Determine a chemical formula.

        chemicalFormula = None;

        # We assume here that the supplied resultSet has come from one of the routines in this module, and thus that all structures in each set of spacegroup groups have the same composition.

        spacegroupGroups = structureSet.GetStructureSet();

        for structures, degeneracies in spacegroupGroups.values():
            chemicalFormula = structures[0].GetChemicalFormula(atomicSymbolLookupTable = atomicSymbolLookupTable);
            break;

        # Build a name for the archive.

        # If a prefix has been set, start with that.

        archiveName = "{0}_".format(prefix) if prefix != None else "";

        # Add the substitution number to the name.

        archiveName = archiveName + "{0:0>3}_".format(i + 1);

        # Finally, append the chemical formula and the file extension.

        archiveName = archiveName + chemicalFormula + ".tar.gz";

        # Export the structure set.

        StructureSetIO.ExportStructureSet(
            structureSet, os.path.join(archiveDirectory, archiveName), fileFormat = fileFormat, atomicSymbolLookupTable = atomicSymbolLookupTable, **kwargs
            );

        # Yield the structure set back to the caller.

        yield structureSet;


def ImportResultSet(prefix = None, archiveDirectory = "./", **kwargs):
    # If prefix is supplied, underscores are removed, if present, to mirror ExportAtomicSubstitutionResultSet().

    if prefix != None:
        if '_' in prefix:
            warnings.warn("Underscores in prefix will be converted to hyphens.", UserWarning);

            prefix = prefix.replace('_', '-');

    # Search archiveDirectory for .tar.gz files.

    inputFiles = [];

    for entry in os.listdir(archiveDirectory):
        absPath = os.path.join(archiveDirectory, entry);

        if os.path.isfile(absPath):
            if entry[-7:].lower() == ".tar.gz":
                inputFiles.append(entry);

    # Attempt to parse the file names and group them into result sets.

    resultSetGroups = { };

    # The format of the file names saved by ExportResultSet is "[<prefix>_]<number>_<chemical_formula>.tar.gz".

    for archiveFile in inputFiles:
        # Trim the .tar.gz extension and split at the underscore character.

        components = archiveFile[:-7].split('_');

        archivePrefix, archiveNumber, archiveChemicalFormula = None, None, None;

        if len(components) == 2:
            if components[0].isdigit():
                # Two elements -> archive number + chemical formula.

                archiveNumber = int(components[0]);
                archiveChemicalFormula = components[1];
            else:
                continue;

        elif len(components) == 3:
            if components[1].isdigit():
                # Three elements -> archive prefix, number and chemical formula.

                archivePrefix = components[0];
                archiveNumber = int(components[1]);
                archiveChemicalFormula = components[2];
            else:
                continue;

        else:
            continue;

        # Add prefix to resultSets if required.

        if archivePrefix not in resultSetGroups:
            resultSetGroups[archivePrefix] = { };

        # Set a key from the archive number and chemical formula.

        key = (archiveNumber, archiveChemicalFormula);

        # If the key is already present, it means archiveDirectory contains archives of multiple result sets that can't be separated.

        if key in resultSetGroups[archivePrefix]:
            # If a prefix was not supplied, or the archive prefix is equal to the target prefix, we cannot work out what to do without user input -> throw an error.

            if prefix == None or archivePrefix == prefix:
                raise Exception("Error: Multiple result sets in archive directory \"{0}\" cannot be separated - please specify the prefix manually or remove unwanted result sets.".format(archiveDirectory));

        resultSetGroups[archivePrefix][key] = archiveFile;

    # Check result sets were found.

    if len(resultSetGroups) == 0:
        raise Exception("Error: No result set archives found in archive directory \"{0}\".".format(archiveDirectory));

    if prefix != None:
        # If prefix is specified, check archives with that prefix were found.

        if prefix not in resultSetGroups:
            raise Exception("Error: Archive files with the prefix \"{0}\" were not found in archive directory \"{1}\".".format(prefix, archiveDirectory));
    else:
        # If not, check we only found archives with one prefix.

        if len(resultSetGroups) > 1:
            raise Exception("Error: Result set archives with multiple prefixes were found in archive directory \"{0}\" - please specify a prefix via the prefix keyword.".format(archiveDirectory));

        prefix = None;

        for key in resultSetGroups.keys():
            prefix = key;
            break;

    # Finally, check the group contain archives with the same number and different chemical formulae, and that the result sets are numbered sequentially from 1.

    resultSetGroup = resultSetGroups[prefix];

    archiveNumbers = [];

    for archiveNumber, _ in resultSetGroup.keys():
        if archiveNumber in archiveNumbers:
            raise Exception("Error: Archive directory \"{0}\" appears to contain multiple sets of results with the same prefix - please check.".format(archiveDirectory));

    archiveNumbers.sort();

    for i in range(0, len(archiveNumbers)):
        if archiveNumbers[i] != i + 1:
            raise Exception("Error: Archives appear to be missing from the specified result set in archive directory \"{0}\".".format(archiveDirectory));

    # Read archives.

    resultSet = [
        StructureSetIO.ImportStructureSet(os.path.join(archiveDirectory, resultSetGroup[key]), **kwargs)
            for key in sorted(resultSetGroup.keys(), key = lambda item : item[0])
        ];

    # Return result set.

    return resultSet;
