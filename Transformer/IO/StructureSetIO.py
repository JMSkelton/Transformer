# Transformer/IO/StructureSetIO.py


# -------
# Imports
# -------

import fractions;
import io;
import os;
import re;
import sys;
import tarfile;
import warnings;

from Transformer.IO import StructureIO;

from Transformer.StructureSet import StructureSet;

from Transformer.Utilities import MultiprocessingHelper;

# Try to import the tqdm module for displaying progress bars in the ImportStructureSet routine.

_TQDM = False;

try:
    import tqdm;

    _TQDM = True;
except ImportError:
    pass;


# ---------------------------
# StructureSet Export Routine
# ---------------------------

def ExportStructureSet(structureSet, archivePath, fileFormat = 'vasp', atomicSymbolLookupTable = None, spacegroupSubfolders = False, normaliseDegeneracies = False):
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

                # We want to write out the structure and add it to the archive without going via a temporary file, as this is slow and could make the function unsafe for multithreaded use.
                # To do this, we encode the file text into an in-memory BytesIO stream, which behaves like a file opened in binary ('rb') mode, as required by tarfile.addfile().

                fileObj = None;

                if sys.version_info.major < 3:
                    # In Python 2.x, strings are represented by default in binary (as a bytes object), and can be be written directly into a BytesIO stream.

                    fileObj = io.BytesIO();

                    StructureIO.WriteStructure(
                        structure, fileObj, fileFormat = fileFormat, atomicSymbolLookupTable = atomicSymbolLookupTable
                        );

                    # Seek back to the beginning after writing.

                    fileObj.seek(0);
                else:
                    # Newer versions of Python store strings as Unicode objects, which need to be encoded before being written into a binary stream.
                    # We therefore first compile the file text in a temporary StringIO buffer, then encode the value to unicode and build a BytesIO around the resulting binary data.
                    # Although this involves holding two copies of the file contents in memory, it's still significantly faster than writing and reading a temporary file.

                    with io.StringIO() as stringBuffer:
                        StructureIO.WriteStructure(
                            structure, stringBuffer, fileFormat = fileFormat, atomicSymbolLookupTable = atomicSymbolLookupTable
                            );

                        fileObj = io.BytesIO(
                            stringBuffer.getvalue().encode('utf8')
                            );

                try:
                    fileName = fileNameTemplate.format(i + 1);

                    fileInfo = tarfile.TarInfo(
                        name = "{0}/{1}/{2}".format(archiveName, subfolderName, fileName) if subfolderName != None else "{0}/{1}".format(archiveName, fileName)
                        );

                    fileInfo.size = len(
                        fileObj.getvalue()
                        );

                    archiveFile.addfile(
                        tarinfo = fileInfo, fileobj = fileObj
                        );
                finally:
                    if fileObj != None:
                        # Release the memory held by the BytesIO stream.

                        fileObj.close();


# ---------------------------
# StructureSet Import Routine
# ---------------------------

_ImportStructureSet_StructureNameRegex = re.compile(r"(?P<chemical_formula>[a-zA-Z0-9]+) \: SG \= (?P<space_group_number>\d+) \((?P<space_group_symbol>[a-zA-Z0-9/_-]+)\)\, (rel\. weight|degeneracy) \= (?P<degeneracy>\d+)");

def _ImportStructureSet_PerformSymmetryAnalysis(structure):
    structure.GetSpacegroup();

    return structure;

def ImportStructureSet(filePath, fileFormat = None, atomTypeNumberLookupTable = None, progressBar = False, useMP = False, mpNumProcesses = None):
    # If progressBar is set and the tqdm module is not available, issue a RuntimeWarning and reset it.

    if progressBar and not _TQDM:
        warings.warn("The tqdm module could not be imported -> progressBar will be reset to False.", RuntimeWarning);

        progressBar = False;

    structures, degeneracies = [], [];

    with tarfile.open(filePath, 'r:gz') as archiveFile:
        # Loop over file paths in the archive.

        members = archiveFile.getmembers();

        iValues = range(0, len(members));

        if progressBar:
            iValues = tqdm.tqdm(iValues);

        for i in iValues:
            member = members[i];

            if not member.isfile():
                # Skip any members that are not files.

                continue;

            # Extract the file name from the path.

            fileName = os.path.split(member.name)[-1];

            # If a file format is not supplied, try to determine one from the file name.
            # If the archive was written using ExportStructureSet(), the default file extensions should have been used, so it shouldn't be necessary to supply a file format.

            fileFormatCurrent = None;

            if fileFormat == None:
                fileFormatCurrent = StructureIO._GetCheckFileFormat(fileName, fileFormat, 'r');

            # Extract the file text into an in-memory buffer.

            stringBuffer = None;

            # tarfile.extractfile() can only be used with a with statement in newer versions of Python.

            inputReader = archiveFile.extractfile(member.name);

            try:
                stringBuffer = io.StringIO(
                    inputReader.read().decode('utf8')
                    );
            finally:
                inputReader.close();

            try:
                # Read the structure from the memory stream.

                structure = StructureIO.ReadStructure(
                    stringBuffer, fileFormat = fileFormatCurrent, atomTypeNumberLookupTable = atomTypeNumberLookupTable
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

                        stringBuffer.seek(0);

                        for line in stringBuffer:
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
                # Make sure the buffer is disposed of.

                if stringBuffer != None:
                    stringBuffer.close();

    # Constructing a StructureSet object from the list of structures requires a symmetry analysis to be performed on each structure to compute the spacegroups.
    # We use the QueueMapFunction() in the MultiprocessingHelper module to run through each structure through a dummy _ImportStructureSet_PerformSymmetryAnalysis() function to trigger the symmetry analysis and cache the result.
    # If useMP is set, the symmetry analyses can be performed in parallel; if not, setting maxNumProcesses = 1 causes QueueMapFunction() to fall back to a serial routine.

    structures = MultiprocessingHelper.QueueMapFunction(
        _ImportStructureSet_PerformSymmetryAnalysis, structures, maxNumProcesses = mpNumProcesses if useMP else 1, progressBar = progressBar
        );

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
