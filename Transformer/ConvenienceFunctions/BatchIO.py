# Transformer/ConvenienceFunctions/BatchIO.py by J. M. Skelton


# -------
# Imports
# -------

import os;
import re;
import tarfile;
import warnings;

from Transformer import IO;
from Transformer import Utilities;


# ----------------
# Module Constants
# ----------------

ExportFileFormats = {
    'vasp' : '.vasp',
    'aims' : '.geometry.in'
    };


# --------------------------
# Result Set Export Function
# --------------------------

_ExportResultSet_TemporaryFileName = r"_ExportTemp.tmp";

def ExportResultSet(resultSet, prefix = None, atomicSymbolLookupTable = None, workingDirectory = "./", fileFormat = 'vasp', spacegroupSubfolders = False):
    fileFormat = fileFormat.lower();

    # Check fileFormat.

    if fileFormat not in ExportFileFormats:
        raise Exception("Error: fileFormat '{0}' is not supported.".format(fileFormat));

    # prefix should not contain underscores; if it does, issue a warning and replace them with hyphens.

    if prefix != None:
        if '_' in prefix:
            warnings.warn("Underscores in prefix will be converted to hyphens.", UserWarning);

            prefix = prefix.replace('_', '-');

    # If workingDirectory does not exist, create it.

    if not os.path.isdir(workingDirectory):
        os.makedirs(workingDirectory);

    # Loop over items in the result set.

    for i, spacegroupGroups in enumerate(resultSet):
        # Sort the spacegropup keys and reorder them to descending symmetry order.

        keys = sorted(spacegroupGroups.keys())[::-1];

        # Determine a chemical formula.
        # We assume here that the supplied resultSet has come from one of the routines in this module, and thus that all structures in each set of spacegroup groups have the same composition.

        structures, _ = spacegroupGroups[keys[0]];

        chemicalFormula = structures[0].GetChemicalFormula(atomicSymbolLookupTable = atomicSymbolLookupTable);

        # Build a name for the archive.

        # If a prefix has been set, start with that.

        archiveName = "{0}_".format(prefix if prefix != None else "");

        # Add the substitution number to the name.

        archiveName = archiveName + "{0:0>3}_".format(i + 1);

        # Finally, append the chemical formula and the file extension.

        archiveName = archiveName + chemicalFormula;

        # Calculate the common divisor to normalise the degeneracies.

        mergedDegeneracies = [];

        for key in keys:
            _, degeneracies = spacegroupGroups[key];
            mergedDegeneracies = mergedDegeneracies + degeneracies;

        commonDivisor = Utilities.GetCommonDivisor(mergedDegeneracies);

        # Output and archive the structures.

        with tarfile.open(os.path.join(workingDirectory, "{0}.tar.gz".format(archiveName)), 'w:gz') as archiveFile:
            for key in keys:
                spacegroupNumber, spacegroupSymbol = key;

                # If spacegroupSubfolders is set, the structures will be divided into spacegroup subfolders.

                subfolderName = None;

                if spacegroupSubfolders:
                    subfolderName = "{0}-{1}".format(
                        spacegroupNumber, spacegroupSymbol.replace('/', '_')
                        );

                structures, degeneracies = spacegroupGroups[key];

                # Write out each structure and add to the archive.

                for i, (structure, degeneracy) in enumerate(zip(structures, degeneracies)):
                    # Give each structure a title line that includes the chemical formula, spacegroup and normalised degeneracy.

                    structure.SetName(
                        "{0} : SG = {1} ({2}), rel. weight = {3}".format(chemicalFormula, spacegroupNumber, spacegroupSymbol, degeneracy // commonDivisor)
                        );

                    # Generate a file name from the chemical formula, spacegroup number and structure number.

                    fileName = "{0}_SG-{1}_{2:0>4}{3}".format(chemicalFormula, spacegroupNumber, i + 1, ExportFileFormats[fileFormat]);

                    if fileFormat == 'vasp':
                        IO.WritePOSCARFile(structure, _ExportResultSet_TemporaryFileName);
                    elif fileFormat == 'aims':
                        IO.WriteAIMSGeometryFile(structure, _ExportResultSet_TemporaryFileName);

                    archiveFile.add(
                        _ExportResultSet_TemporaryFileName, arcname = "{0}/{1}/{2}".format(archiveName, subfolderName, fileName) if subfolderName != None else "{0}/{1}".format(archiveName, fileName)
                        );

                    # Delete the temporary file once added to the archive.

                    os.remove(_ExportResultSet_TemporaryFileName);


# ---------------------------
# Result Set Import Functions
# ---------------------------

_ImportResultSetArchive_StructureNameRegex = re.compile(r"(?P<chemical_formula>[a-zA-Z0-9]+) \: SG \= (?P<space_group_number>\d+) \((?P<space_group_symbol>[a-zA-Z0-9/_-]+)\)\, rel\. weight \= (?P<degeneracy>\d+)");

_ImportResultSetArchive_TemporaryFileName = r"_ImportTemp.tmp";

def ImportResultSetArchive(filePath):
    structures, degeneracies = [], [];

    with tarfile.open(filePath, 'r:gz') as archiveFile:
        # Loop over file paths in the archive.

        for archivePath in archiveFile.getnames():
            # Extract the file name from the path.

            fileName = os.path.split(archivePath)[-1];

            # Infer the file type from the extension and check support.

            fileType = None;

            for key, fileExtension in ExportFileFormats.items():
                if len(fileName) > len(fileExtension) and fileName[-len(fileExtension):].lower() == fileExtension:
                    fileType = key;
                    break;

            if fileType not in ExportFileFormats:
                raise Exception("Error: Archive file \"{0}\" contains files with an unknown file type.'".format(filePath));

            # Temporarily extract the file to read in.

            with archiveFile.extractfile(archivePath) as inputReader:
                with open(_ImportResultSetArchive_TemporaryFileName, 'wb') as outputWriter:
                    outputWriter.write(inputReader.read());

            # Read the structure from the file.

            structure = None;

            if fileType == 'vasp':
                structure = IO.ReadPOSCARFile(_ImportResultSetArchive_TemporaryFileName);
            elif fileType == 'aims':
                structure = IO.ReadAIMSGeometryFile(_ImportResultSetArchive_TemporaryFileName);

            structures.append(structure);

            # If we have not stopped trying to retrieve degeneracies, attempt to retrieve one from the structure or input file, depending on the file type.

            if degeneracies != None:
                degeneracy = None;

                if fileType == 'vasp':
                    # For file formats that can store the structure name, the degeneracy can be obtained from the structure object.

                    match = _ImportResultSetArchive_StructureNameRegex.search(structure.GetName());

                    if match:
                        degeneracy = int(match.group('degeneracy'));

                elif fileType == 'aims':
                    # For other file formats, the name is written as a comment.

                    with open(_ImportExportTemporaryFileName, 'r') as inputReader:
                        for line in inputReader:
                            match = _ImportResultSetArchive_StructureNameRegex.search(line);

                            if match:
                                degeneracy = int(match.group('degeneracy'));
                                break;

                if degeneracy == None:
                    # If we fail to retrieve a degeneracy, stop attempting to retrieve them anad issue a warning.

                    warnings.warn("Degeneracies could not be extracted from one or more input files -> all degeneracies will be set to 1.", UserWarning);

                    degeneracies = None;
                else:
                    # If not, add the degeneracy to the list.

                    degeneracies.append(degeneracy);

            # Remove the temporary file.

            os.remove(_ImportResultSetArchive_TemporaryFileName);

        # If degeneracies could not be retrieved, set it to a list of ones.

        if degeneracies == None:
            degeneracies = [1] * len(structures);

    # Return the lists of structures and degeneracies.

    return structures, degeneracies;

def ImportResultSet(prefix = None, directory = "./"):
    # If prefix is supplied, underscores are removed, if present, to mirror ExportAtomicSubstitutionResultSet().

    if prefix != None:
        if '_' in prefix:
            warnings.warn("Underscores in prefix will be converted to hyphens.", UserWarning);

            prefix = prefix.replace('_', '-');

    # Search directory for .tar.gz files.

    inputFiles = [];

    for entry in os.listdir(directory):
        absPath = os.path.join(directory, entry);

        if os.path.isfile(absPath):
            if entry[-7:].lower() == ".tar.gz":
                inputFiles.append(entry);

    # The format of the file names saved by ExportAtomicSubstitutionResultSet is "[<prefix>_]<number>_<chemical_formula>.tar.gz".

    resultSetExports = { };

    for archiveFile in inputFiles:
        # Trim the .tar.gz extension and split at the underscore character.

        components = archiveFile[:-7].split('_');

        archivePrefix, archiveNumber, archiveChemicalFormula = None, None, None;

        if len(components) == 2:
            # Two elements -> archive number + chemical formula.

            if components[0].isdigit():
                archiveNumber = int(components[0]);
                archiveChemicalFormula = components[1];
            else:
                continue;

        elif len(components) == 3:
            # Three elements -> archive prefix, number and chemical formula.

            if components[1].isdigit():
                archivePrefix = components[0];
                archiveNumber = int(components[1]);
                archiveChemicalFormula = components[2];
            else:
                continue;

        else:
            continue;

        # Add prefix to resultSets if required.

        if archivePrefix not in resultSetExports:
            resultSetExports[archivePrefix] = { };

        # Set a key from the archive number and chemical formula.

        key = (archiveNumber, archiveChemicalFormula);

        # If the key is already present, it means directory contains archives of multiple result sets that can't be separated.

        if key in resultSetExports[archivePrefix]:
            # If a prefix was not supplied, or the archive prefix is equal to the target prefix, we cannot work out what to do without user input -> throw an error.

            if prefix == None or archivePrefix == prefix:
                raise Exception("Error: Multiple result sets in input directory \"{0}\" cannot be separated - please specify the prefix manually or remove unwanted result sets from the input directory.".format(directory));

        resultSetExports[archivePrefix][key] = archiveFile;

    # Check result sets were found.

    if len(resultSetExports) == 0:
        raise Exception("Error: No result set archive files found in input directory \"{0}\".".format(directory));

    if prefix != None:
        # If prefix is specified, check archives with that prefix were found.

        if prefix not in resultSetExports:
            raise Exception("Error: Archive files with the prefix \"{0}\" were not found in input directory \"{1}\".".format(prefix, directory));
    else:
        # If not, check we only found archives with one prefix.

        if len(resultSetExports) > 1:
            raise Exception("Error: Result set archives with multiple prefixes were found in input directory \"{0}\" - please specify a prefix via the prefix keyword.".format(directory));

        prefix = [key for key in resultSetExports.keys()][0];

    # Finally, check the result sets are numbered sequentially and do not contain archives with the same number and different chemical formulae.

    resultSetExport = resultSetExports[prefix];

    archiveNumbers = [];

    for archiveNumber, _ in resultSetExport.keys():
        if archiveNumber in archiveNumbers:
            raise Exception("Error: Input directory \"{0}\" appears to contain multiple sets of results with the same prefix - please check.".format(directory));

    # Extract data from archives.

    resultSet = [];

    for key in sorted(resultSetExport.keys(), key = lambda item : item[0]):
        # Import the archive file.

        archiveFileName = resultSetExport[key];

        structures, degeneracies = ImportResultSetArchive(
            os.path.join(directory, archiveFileName)
            );

        # Add the structures and degeneracies to the result set.

        resultSet.append(
            (structures, degeneracies)
            );

    # Return the result set.

    return resultSet;
