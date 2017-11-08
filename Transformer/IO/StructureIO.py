# Transformer/IO/StructureIO.py


# -------
# Imports
# -------

from Transformer.IO import _AIMS, _VASP;


# ---------
# Constants
# ---------

# Tuples of (formatCode, defaultExtension, readSupport, writeSupport).

SupportedFileFormats = [
    ('aims', '.geometry.in', True, True),
    ('vasp',        '.vasp', True, True)
    ];


# ----------------
# Public Functions
# ----------------

def ReadStructure(filePath, fileFormat = None, atomTypeNumberLookupTable = None):
    if fileFormat == None:
        # If a file format is not supplied, try and determine one from the file path.

        fileFormat = _TryGetFileFormat(filePath);
    else:
        # If one is, convert it to lower case to match the convention in SupportedFileFormats.

        fileFormat = fileFormat.lower();

    # Check the selected file format is supported for reading.

    for supportedFileFormat, _, readSupport, _ in SupportedFileFormats:
        if fileFormat == supportedFileFormat and not readSupport:
            raise Exception("Error: File format '{0}' is not supported for reading.".format(fileFormat));

    # Dispatch to different reader functions depending on the selected file format.

    if fileFormat == 'aims':
        return _AIMS.ReadGeometryInFile(filePath, atomTypeNumberLookupTable = atomTypeNumberLookupTable);

    elif fileFormat == 'vasp':
        return _VASP.ReadPOSCARFile(filePath, atomTypeNumberLookupTable = atomTypeNumberLookupTable);

    else:
        # Catch-all, just in case.

        raise NotImplementedError("Error: An import routine for the file format '{0}' has not yet been implemented.".format(fileFormat));

def WriteStructure(structure, filePath, fileFormat = None, atomicSymbolLookupTable = None):
    # Determine/lower case fileFormat.

    if fileFormat == None:
        fileFormat = _TryGetFileFormat(filePath);
    else:
        fileFormat = fileFormat.lower();

    # Check the format is supported for writing.

    for supportedFileFormat, _, _, writeSupport in SupportedFileFormats:
        if fileFormat == supportedFileFormat and not writeSupport:
            raise Exception("Error: File format '{0}' is not supported for writing.".format(fileFormat));

    # Dispatch to the appropriate writer function.

    if fileFormat == 'aims':
        _AIMS.WriteGeometryInFile(
            structure, filePath, atomicSymbolLookupTable = atomicSymbolLookupTable
            );

    elif fileFormat == 'vasp':
        _VASP.WritePOSCARFile(
            structure, filePath, atomicSymbolLookupTable = atomicSymbolLookupTable
            );

    else:
        raise NotImplementedError("Error: An export routine for the file format '{0}' has not yet been implemented.".format(fileFormat));


# -----------------
# Utility Functions
# -----------------

def _TryGetFileFormat(filePath):
    # Check filePath against the default extensions for the supported file types listed in SupportedFileFormats.

    for supportedFileFormat, defaultExtension, _, _ in SupportedFileFormats:
        if filePath.lower().endswith(defaultExtension):
            return supportedFileFormat;

    # If the file type cannot be determined, raise an error.

    raise Exception("Error: A file format could not be determined from the file path.");
