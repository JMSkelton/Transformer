# Transformer/IO/StructureIO.py


# -------
# Imports
# -------

import os;

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

def GetFileTypeFromExtension(filePath):
    fileFormat = _TryGetFileFormat(filePath);

    for supportedFileFormat, defaultExtension, _, _ in SupportedFileFormats:
        if fileFormat == supportedFileFormat:
            return (fileFormat, defaultExtension);

def ReadStructure(filePathOrObj, fileFormat = None, atomicSymbolLookupTable = None):
    # _GetCheckFileFormat() attempts to automatically determine a file format if one is not supplied, and also checks the format is supported for reading.

    fileFormat = _GetCheckFileFormat(filePathOrObj, fileFormat, mode = 'r');

    # File-like object to read from.

    fileObj = None;

    # Variable to keep track of whether fileObj was opened within this routine.

    fileObjOpened = False;

    # If filePathOrObj is a string, assume it specifies a file path and set fileObj to an open file for reading.
    # If not, assume filePathOrObj implements a file-like interface and use as is.

    if isinstance(filePathOrObj, str):
        fileObj = open(filePathOrObj, 'r');
        fileObjOpened = True;
    else:
        fileObj = filePathOrObj;

    try:
        # Dispatch to different reader functions depending on the selected file format.

        if fileFormat == 'aims':
            return _AIMS.ReadGeometryInFile(
                fileObj, atomicSymbolLookupTable = atomicSymbolLookupTable
                );

        elif fileFormat == 'vasp':
            return _VASP.ReadPOSCARFile(
                fileObj, atomicSymbolLookupTable = atomicSymbolLookupTable
                );

        else:
            # Catch-all, just in case.

            raise NotImplementedError("Error: An import routine for the file format '{0}' has not yet been implemented.".format(fileFormat));

    finally:
        # If we opened a file, make sure it gets closed.

        if fileObjOpened:
            fileObj.close();

def WriteStructure(structure, filePathOrObj, fileFormat = None, atomicSymbolLookupTable = None):
    # Determine and/or check the file format is supported for writing.

    fileFormat = _GetCheckFileFormat(filePathOrObj, fileFormat, mode = 'w');

    # If required, set up a file-like object to write to.

    fileObj = None;
    fileObjOpened = False;

    if isinstance(filePathOrObj, str):
        fileObj = open(filePathOrObj, 'w');
        fileObjOpened = True;
    else:
        fileObj = filePathOrObj;

    # Dispatch to the appropriate writer function.

    try:
        if fileFormat == 'aims':
            _AIMS.WriteGeometryInFile(
                structure, fileObj, atomicSymbolLookupTable = atomicSymbolLookupTable
                );

        elif fileFormat == 'vasp':
            _VASP.WritePOSCARFile(
                structure, fileObj, atomicSymbolLookupTable = atomicSymbolLookupTable
                );

        else:
            raise NotImplementedError("Error: An export routine for the file format '{0}' has not yet been implemented.".format(fileFormat));

    finally:
        if fileObjOpened:
            fileObj.close();


# -----------------
# Utility Functions
# -----------------

def _GetCheckFileFormat(filePathOrObj, fileFormat, mode):
    if fileFormat == None:
        # If filePathOrObj is a string, assume it specifies a file path and match the ending against the default extensions of the formats listed in SupportedFileFormats.

        if isinstance(filePathOrObj, str):
            _, tail = os.path.split(filePathOrObj);

            tail = tail.lower();

            for supportedFileFormat, defaultExtension, _, _ in SupportedFileFormats:
                if tail.endswith(defaultExtension):
                    fileFormat = supportedFileFormat;

        # If we were unable to determine a file format automatically, throw an error.

        if fileFormat == None:
            raise Exception("Error: A file format could not be automatically determined.");
    else:
        fileFormat = fileFormat.lower();

    # Check the file format is supported for the desired mode.

    if mode == 'r':
        # Check for reading support.

        for supportedFileFormat, _, readSupport, _ in SupportedFileFormats:
            if fileFormat == supportedFileFormat and not readSupport:
                raise Exception("Error: File format '{0}' is not supported for reading.".format(fileFormat));
    elif mode == 'w':
        # Check for writing support.

        for supportedFileFormat, _, _, writeSupport in SupportedFileFormats:
            if fileFormat == supportedFileFormat and not writeSupport:
                raise Exception("Error: File format '{0}' is not supported for writing.".format(fileFormat));
    else:
        # Catch all - should only be for debugging purposes.

        raise Exception("Error: Unknown mode '{0}'.".format(mode))

    # Return the file format.

    return fileFormat;
