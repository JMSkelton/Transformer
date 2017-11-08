# Transformer/IO/StructureSetIO.py


# -------
# Imports
# -------

import fractions;
import os;
import re;
import tarfile;
import warnings;

from Transformer.IO import StructureIO;

from Transformer.StructureSet import StructureSet;


# ---------------------------------
# StructureSet Batch-Export Routine
# ---------------------------------

_ExportStructureSet_TemporaryFileName = r"_ExportTemp.tmp";

def ExportStructureSet(structureSet, archivePath, fileFormat = 'vasp', atomicSymbolLookupTable = None, spacegroupSubfolders = False, normaliseDegeneracies = False):
    # Check the temporary file is not present.

    if os.path.isfile(_ExportStructureSet_TemporaryFileName):
        raise Exception("Error: Temporary file \"{0}\" already exists.".format(_ExportStructureSet_TemporaryFileName));

    # Check file format is supported and get the default extension.

    fileFormat = fileFormat.lower();

    fileExtension = None;

    for supportedFileFormat, defaultExtension, _, writeSupport in StructureIO.SupportedFileFormats:
        if fileFormat == supportedFileFormat and writeSupport:
            fileExtension = defaultExtension;
            break;

    if fileExtension == None:
        raise Exception("Error: File format '{0}' is not supported for writing.".format(fileFormat));

    # Get the archive name from the path.

    archiveName = None;

    _, tail = os.path.split(archivePath);

    if tail.lower().endswith(".tar.gz"):
        archiveName = tail[:-7];
    else:
        archiveName, _ = os.path.splitext(tail);

    # If normaliseDegeneracies is set, calculate a common divisor to normalise the degeneracies.

    commonDivisor = None;

    if normaliseDegeneracies:
        _, degeneraciesFlat = structureSet.GetStructureSetFlat();
        commonDivisor = _GetCommonDivisor(degeneraciesFlat);

    # Output and archive the structures.

    with tarfile.open(archivePath, 'w:gz') as archiveFile:
        spacegroupGroups = structureSet.GetStructureSet();

        for (spacegroupNumber, spacegroupSymbol), (structures, degeneracies) in spacegroupGroups.items():
            # If spacegroupSubfolders is set, the structures will be divided into spacegroup subfolders.

            subfolderName = None;

            if spacegroupSubfolders:
                subfolderName = "{0}-{1}".format(
                    spacegroupNumber, spacegroupSymbol.replace('/', '_')
                    );

            # Generate a template for the structure title lines.
            # Each structure will be gien a name that includes the chemical formula, spacegroup and normalised degeneracy.

            titleLineTemplate = None;

            if commonDivisor != None:
                titleLineTemplate = "{{0}} : SG = {0} ({1}), rel. weight = {{1}}".format(spacegroupNumber, spacegroupSymbol);
            else:
                titleLineTemplate = "{{0}} : SG = {0} ({1}), degeneracy = {{1}}".format(spacegroupNumber, spacegroupSymbol);

            # Generate a template for the file names.
            # Each file will be named with the spacegroup number and structure number.

            fileNameTemplate = "SG-{0:0>3}_{{0:0>4}}{1}".format(spacegroupNumber, fileExtension);

            # Write out each structure and add to the archive.

            for i, (structure, degeneracy) in enumerate(zip(structures, degeneracies)):
                # Set the structure name (title line).

                if commonDivisor != None:
                    degeneracy = degeneracy // commonDivisor;

                structure.SetName(
                    titleLineTemplate.format(structure.GetChemicalFormula(atomicSymbolLookupTable = atomicSymbolLookupTable), degeneracy)
                    );

                # Write out the structure to a temporary file.

                fileName = fileNameTemplate.format(i + 1);

                StructureIO.WriteStructure(
                    structure, _ExportStructureSet_TemporaryFileName, fileFormat = fileFormat, atomicSymbolLookupTable = atomicSymbolLookupTable
                    );

                archiveFile.add(
                    _ExportStructureSet_TemporaryFileName, arcname = "{0}/{1}/{2}".format(archiveName, subfolderName, fileName) if subfolderName != None else "{0}/{1}".format(archiveName, fileName)
                    );

                # Delete the temporary file once added to the archive.

                os.remove(_ExportStructureSet_TemporaryFileName);


# ---------------------------------
# StructureSet Batch-Import Routine
# ---------------------------------

_ImportStructureSet_TemporaryFileName = r"_ImportTemp.tmp";

_ImportStructureSet_StructureNameRegex = re.compile(r"(?P<chemical_formula>[a-zA-Z0-9]+) \: SG \= (?P<space_group_number>\d+) \((?P<space_group_symbol>[a-zA-Z0-9/_-]+)\)\, (rel\. weight|degeneracy) \= (?P<degeneracy>\d+)");

def ImportStructureSet(filePath, fileFormat = None, atomTypeNumberLookupTable = None):
    # Check the temporary file is not present.

    if os.path.isfile(_ImportStructureSet_TemporaryFileName):
        raise Exception("Error: Temporary file \"{0}\" already exists.".format(_ImportStructureSet_TemporaryFileName));

    structures, degeneracies = [], [];

    with tarfile.open(filePath, 'r:gz') as archiveFile:
        # Loop over file paths in the archive.

        for archivePath in archiveFile.getnames():
            # Extract the file name from the path.

            _, fileName = os.path.split(archivePath);

            # If a file format is not supplied, try to determine one from the file name.
            # If the archive was written using ExportStructureSet(), the default file extensions should have been used, so it shouldn't be necessary to supply a file format.

            fileFormatCurrent = None;

            if fileFormat == None:
                fileFormatCurrent = StructureIO._TryGetFileFormat(fileName);

            # Temporarily extract the file to read in.

            with archiveFile.extractfile(archivePath) as inputReader:
                with open(_ImportStructureSet_TemporaryFileName, 'wb') as outputWriter:
                    outputWriter.write(inputReader.read());

            # Read the structure from the file.

            try:
                structure = StructureIO.ReadStructure(
                    _ImportStructureSet_TemporaryFileName, fileFormat = fileFormatCurrent, atomTypeNumberLookupTable = atomTypeNumberLookupTable
                    );

                structures.append(structure);

                # If we have not stopped trying to retrieve degeneracies, attempt to retrieve one from the structure or input file, depending on the file type.

                if degeneracies != None:
                    degeneracy = None;

                    if fileFormatCurrent == 'vasp':
                        # For file formats that can store the structure name, the degeneracy can be obtained from the structure object.

                        match = _ImportStructureSet_StructureNameRegex.search(structure.GetName());

                        if match:
                            degeneracy = int(match.group('degeneracy'));
                    else:
                        # For other file formats, the name is written as a comment.

                        with open(_ImportStructureSet_TemporaryFileName, 'r') as inputReader:
                            for line in inputReader:
                                match = _ImportStructureSet_StructureNameRegex.search(line);

                                if match:
                                    degeneracy = int(match.group('degeneracy'));
                                    break;

                    # If we were able to retrieve a degeneracy, add it to the list.
                    # If not, stop attempting to retrieve them anad issue a warning.

                    if degeneracy != None:
                        degeneracies.append(degeneracy);
                    else:
                        warnings.warn("Degeneracies could not be extracted from one or more input files -> all degeneracies will be set to 1.", UserWarning);

                        degeneracies = None;

            finally:
                # Make sure the temporary file is removed if an error occurs while reading.

                if os.path.isfile(_ImportStructureSet_TemporaryFileName):
                    os.remove(_ImportStructureSet_TemporaryFileName);

    # Return the strucutures and degeneracies as a StructureSet object.
    # Since this function expects to read an exported structure set, we skip the initial merging when constructing the set.

    return StructureSet(structures = structures, degeneracies = degeneracies, noInitialMerge = True);


# -----------------
# Utility Functions
# -----------------

# This code is based on the discussion in https://www.rookieslab.com/posts/cpp-python-code-to-find-gcd-of-a-list-of-numbers.

def _GetCommonDivisor(integers):
    divisor = integers[0];

    for i in integers[1:]:
        divisor = fractions.gcd(divisor, i);

    return divisor;
