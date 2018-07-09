# Transformer/Structure.py


# ----------------
# Module Docstring
# ----------------

""" Contains the primary Structure class used to represent crystal structures. """


# -------
# Imports
# -------

import math;
import re;

import numpy as np;

try:
    # spglib version <= 1.8.x

    import spglib as spg;
except ImportError:
    # Newer spglib versions.

    from pyspglib import spglib as spg;

from Transformer import Constants;


# ---------------
# Structure Class
# ---------------

class Structure:
    """
    Lightweight class for storing crystal structures.

    Internally, the structure data is stored in NumPy arrays to keep a low memory footprint.
    The atom-type numbers and positions are stored in a structured array and kept in sort order, allowing for efficient comparison between structures.

    Class methods expose a core set of symmetry routines and an extensive set of utility functions.

    Properties that are non-trivial to compute or are likely to be accessed repeatedly are automatically cached for performance.

    Usage notes:
        Use the Get*/Set* methods to access and update the structure data; modifying the internal arrays directly may cause other Transformer routines to do unexpected things!
        Most derived properties (e.g. symmetry properties) are returned "as is", and need to be copied before modification.
    """

    # -----------
    # Constructor
    # -----------

    def __init__(self, latticeVectors, atomPositions, atomTypes, name = None, atomicSymbolLookupTable = None):
        """
        Class constructor.

        Arguments:
            latticeVectors -- must be convertible to a 3x3 NumPy matrix.
            atomPositions --must be convertible to an Nx3 NumPy matrix.
            atomTypes -- may be integer atom-type numbers or atomic symbols.

        Keyword arguments:
            name -- a name for the structure.
            atomicSymbolLookupTable -- optional { type_number : symbol } dictionary used to map non-standard atom types to integer type numbers.
        """

        if len(atomTypes) != len(atomPositions):
            raise Exception("Error: The lengths of atomPositions and atomTypes are not consistent.");

        # Convert atom types to atom-type numbers.

        atomTypeNumbers = [
            AtomTypeToAtomTypeNumber(atomType, atomicSymbolLookupTable) for atomType in atomTypes
            ];

        for typeNumber in atomTypeNumbers:
            if typeNumber == None:
                raise Exception("Error: One or more atom types could not be converted to atom-type numbers.");

        # The atom data is stored in a NumPy structured array.
        # We also create typed views into the data for ease of manipulation.

        numAtoms = len(atomTypeNumbers);

        atomData = np.zeros(numAtoms, dtype = Structure._AtomDataType);

        atomTypeNumbersView = atomData.view(dtype = np.int64).reshape((numAtoms, 4))[:, 0];
        atomTypeNumbersView[:] = atomTypeNumbers;

        atomPositionsView = atomData.view(dtype = np.float64).reshape((numAtoms, 4))[:, 1:];
        atomPositionsView[:] = atomPositions;

        # Clamp (fractional) atom positions to the range [0, 1].

        atomPositionsView[:] %= 1.0;

        atomData.sort();

        # Set fields.

        self.SetLatticeVectors(latticeVectors);

        self._atomData = atomData;

        self._atomTypeNumbersView = atomTypeNumbersView;
        self._atomPositionsView = atomPositionsView;

        self._name = name;

        self._ResetLazyPropertyFields();

    # ---------------
    # Private Methods
    # ---------------

    def _ResetLazyPropertyFields(self):
        """ Resets internal fields used for caching/lazy initialisation of non-trivial properties. """

        self._symmetryAnalysisTolerance = None;

        self._pSpacegroup = None;
        self._pSymmetryOperations = None;
        self._pUniqueAtomIndices = None;

        self._pAtomPositionsCartesian = None;

        self._pChemicalFormula = None;
        self._pNeighbourTable = None;

    def _PerformSymmetryAnalysis(self, tolerance, spacegroupOnly = False):
        """ Performs a symmetry analysis using spglib and the supplied tolerance. """

        # If no symmetry tolerance is provided, use the default value.

        if tolerance == None:
            tolerance = Structure.DefaultSymmetryEquivalenceTolerance;

        # Only perform the symmetry analysis if it hasn't been done yet or if the tolerance has changed.

        performAnalysis = tolerance != self._symmetryAnalysisTolerance;

        if not performAnalysis:
            if spacegroupOnly:
                performAnalysis = self._pSpacegroup == None;
            else:
                performAnalysis = self._pSymmetryOperations == None;

        if performAnalysis:
            # if spacegroupOnly is set, we do not perform and cache a full symmetry dataset - for large sets of high-symmetry structures this can consume a lot of memory.

            if spacegroupOnly:
                # Call the get_spacegroup() routine from spglib.

                result = spg.get_spacegroup(
                    (self._latticeVectors, self._atomPositionsView, self._atomTypeNumbersView),
                    symprec = tolerance
                    );

                if result == None:
                    raise Exception("Error: spglib get_spacegroup() routine returned None.");

                match = Structure._SpacegroupRegex.match(result);

                if match:
                    self._pSpacegroup = (
                        (int(match.group('spacegroup_number')), match.group('spacegroup_symbol').strip())
                        );
                else:
                    raise Exception("Error: Failed to parse string returned by spglib get_spacegroup().")

            else:
                # Call the get_symmetry_dataset() routine from spglib.

                result = spg.get_symmetry_dataset(
                    (self._latticeVectors, self._atomPositionsView, self._atomTypeNumbersView),
                    symprec = tolerance
                    );

                if result == None:
                    raise Exception("Error: spglib get_symmetry_dataset() routine returned None.");

                # spglib returns symmetry operations in rotation + translation format.
                # The translation vectors sometimes show numerical noise, which can sometimes cause structures related by symmetry to be detected as inequivalent.
                # We therefore round them based on the set tolerance and clamp to the range [0, 1].

                roundDecimals = -1 * int(
                    math.floor(math.log10(tolerance))
                    );

                translations = [
                    np.round(translation, roundDecimals) % 1.0
                        for translation in result['translations']
                    ];

                self._pSpacegroup = (result['number'], result['international']);
                self._pSymmetryOperations = [item for item in zip(result['rotations'], translations)];

                # Store the unique atom indices along with the number of equivalent sites.

                uniqueAtomIndices = np.unique(result['equivalent_atoms']);

                self._pUniqueAtomIndices = (
                    [index for index in uniqueAtomIndices],
                    [len(np.where(result['equivalent_atoms'] == index)[0]) for index in uniqueAtomIndices]
                    );

            # Store the symmetry tolerance used to call the spglib functions.

            self._symmetryAnalysisTolerance = tolerance;

    # ----------------
    # Property Getters
    # ----------------

    def GetAtomCount(self):
        """ Return the number of atoms in the structure. """

        return len(self._atomData);

    def GetLatticeVectors(self):
        """ Return a copy of the lattice vectors as a list of NumPy vectors. """

        return [np.copy(row) for row in self._latticeVectors[:]];

    def GetAtomTypeNumbers(self):
        """ Return a copy of the atom-type numbers. """

        return [atomTypeNumber for atomTypeNumber in self._atomTypeNumbersView[:]];

    def GetAtomPositions(self):
        """ Return a copy of the atom positions as a list of NumPy vectors. """

        return [np.copy(row) for row in self._atomPositionsView[:]];

    def GetLatticeVectorsNumPy(self, copy = True):
        """
        Return the internal 3x3 NumPy matrix representation of the lattice vectors.

        Keyword arguments:
            copy -- if True (default), return a copy of the data.
        """

        if copy:
            return np.copy(self._latticeVectors);
        else:
            return self._latticeVectors;

    def GetAtomTypeNumbersNumPy(self, copy = True):
        """
        Return the internal NumPy view of the atom-type numbers.

        Keyword arguments:
            copy -- if True (default), return a copy of the data.
        """

        if copy:
            return np.copy(self._atomTypeNumbersView);
        else:
            return self._atomTypeNumbersView;

    def GetAtomPositionsNumPy(self, copy = True):
        """
        Return the internal NumPy view of the atom positions.

        Keyword arguments:
            copy -- if True (default), return a copy of the data.
        """

        if copy:
            return np.copy(self._atomPositionsView);
        else:
            return self._atomPositionsView;

    def GetAtomDataNumPy(self, copy = True):
        """
        Return the internal NumPy representation of the atom data.
        This is a 1D structured array of type Structure._AtomDataType.

        Keyword arguments:
            copy -- if True (default), return a copy of the data.
        """

        if copy:
            return np.copy(self._atomData);
        else:
            return self._atomData;

    def GetName(self):
        """ Return the structure name, if set; otherwise, returns a chemical formula. """

        name = self._name;

        if name != None and name != "":
            return name;
        else:
            # If name is empty, compute and return a chemical formula.

            return self.GetChemicalFormula();

    # ----------------
    # Property Setters
    # ----------------

    def SetName(self, name):
        """ Set a name for the structure. """

        self._name = name;

    def SetLatticeVectors(self, latticeVectors):
        """
        Set the lattice vectors.

        Arguments:
            latticeVectors -- must be convertible to a 3x3 NumPy matrix.
        """

        latticeVectors = np.array(latticeVectors, dtype = np.float64);

        dim1, dim2 = latticeVectors.shape;

        if dim1 != dim2 != 3:
            raise Exception("Error: latticeVectors must be convertible to a 3x3 NumPy array.");

        self._latticeVectors = latticeVectors;

    def SetAtoms(self, indexOrIndices, atomTypes, atomPositions, atomicSymbolLookupTable = None):
        """
        Set type(s) and/or position(s) of selected atoms.

        Arguments:
            indexOrIndices -- an index or list of indices of atoms to update.
            atomTypes -- an atom type or list of atom types to update for the selected atoms.
            atomPositions -- a NumPy vector or list of NumPy vectors to update for the selected atoms.

        Keyword arguments:
            atomicSymbolLookupTable -- { type_number : symbol } dictionary used to map non-standard atom types to integer type numbers.

        Notes:
            If indexOrIndices is a list, atomTypes and atomPositions, if supplied, must also be lists.
            Atom types may be specified by symbols or integer atom-type numbers.
            If an atom type and/or position is set to None, its value is left unchanged.
            Since updating the same atom twice is likely a mistake, doing so raises an error.
        """

        indices = None;

        if isinstance(indexOrIndices, list):
            indices = indexOrIndices;

            # If a list of indices is supplied, it should not contian duplicate indices.

            indicesCheck = set();

            for index in indices:
                if index in indicesCheck:
                    raise Exception("Error: If a list of indices is supplied, it cannot contain duplicate indices.");
                else:
                    indicesCheck.add(index);

            # If a list of indices is supplied, atomTypes and atomPositions, if supplied, must also be lists.

            if atomTypes != None:
                if not (isinstance(atomTypes, list) and len(atomTypes) == len(indices)):
                    raise Exception("Errror: If a list of indices is supplied, atomTypes, if supplied, must be a list and must contain the same number of entries.");
            else:
                atomTypes = [None] * len(indices);

            if atomPositions != None:
                if not (isinstance(atomPositions, list) and len(atomPositions) == len(indices)):
                    raise Exception("Error: If a list of indices is supplied, atomPositions, if supplied, must be a list and must contain the same number of entries.");
            else:
                atomPositions = [None] * len(indices);

        else:
            # If a single index is supplied, assume atomTypes and atomPositions are scalars and convert everything to lists.

            indices = [indexOrIndices];

            atomTypes = [atomTypes];
            atomPositions = [atomPositions];

        numAtoms = len(self._atomData);

        for index in indices:
            if index >= numAtoms:
                raise Exception("Error: Index {0} is out of range for number of atoms {1}.".format(index, numAtoms));

        # Convert atom types to atom-type numbers.

        for i, atomType in enumerate(atomTypes):
            typeNumber = AtomTypeToAtomTypeNumber(atomType, atomicSymbolLookupTable);

            if typeNumber == None:
                raise Exception("Error: Unable to convert atom-type '{0}' to an atom-type number.".format(atomType));

            atomTypes[i] = typeNumber;

        for i, index in enumerate(indices):
            # Set atoms.

            if atomTypes != None and atomTypes[i] != None:
                # Substitute atom.

                self._atomTypeNumbersView[index] = atomTypes[i];

            if atomPositions != None and atomPositions[i] != None:
                # Update atom position.

                self._atomPositionsView[index] = [x % 1.0 for x in atomPositions[i]];

        # Resort the atoms and reset cached property fields.

        self._atomData.sort();

        self._ResetLazyPropertyFields();

    def SwapAtoms(self, atomTypes1, atomTypes2, atomicSymbolLookupTable = None):
        """
        Swaps all atoms of type(s) specified in atomTypes1 for the corresponding type(s) in atomTypes2.

        Arguments:
            atomTypes1 -- an atom type or list of atom types to swap.
            atomTypes2 -- an atom type or list of atom types to swap to.

        Keyword arguments:
            atomicSymbolLookupTable -- { type_number : symbol } dictionary used to map non-standard atom types to integer type numbers.

        Notes:
            If atomTypes1 is a list, atomTypes2 must also be a list.
            Updating the atoms is performed in one go, so it is possible to swap two types of atoms with this method, i.e. atomTypes1 = [a, b] and atomTypes2 = [b, a] will work as expected.
        """

        if atomTypes1 == None:
            raise Exception("Error: atomTypes1 cannot be None.");

        if atomTypes2 == None:
            raise Exception("Error: atomTypes2 cannot be None.")

        # Convert atom types to atom-type numbers.

        atomTypeNumbers1, atomTypeNumbers2 = [], [];

        if isinstance(atomTypes1, list):
            atomTypeNumbers1 = [
                AtomTypeToAtomTypeNumber(atomType, atomicSymbolLookupTable) for atomType in atomTypes1
                ];
        else:
            atomTypeNumbers1 = [AtomTypeToAtomTypeNumber(atomTypes1, atomicSymbolLookupTable)];

        if isinstance(atomTypes2, list):
            atomTypeNumbers2 = [
                AtomTypeToAtomTypeNumber(atomType, atomicSymbolLookupTable) for atomType in atomTypes2
                ];
        else:
            atomTypeNumbers2 = [AtomTypeToAtomTypeNumber(atomTypes2, atomicSymbolLookupTable)];

        if len(atomTypeNumbers1) != len(atomTypeNumbers2):
            raise Exception("Error: atomTypes1 and atomTypes2 must contain the same number of atom types.");

        # Get indices.

        atomTypeNumbers = self._atomTypeNumbersView;

        indices, atomTypes = [], [];

        for typeNumber1, typeNumber2 in zip(atomTypeNumbers1, atomTypeNumbers2):
            # Both type numbers should not be None, and the first type number should be present in the structure.

            if typeNumber1 == None or typeNumber2 == None:
                raise Exception("Error: One or more atom types were not recognised or were set to None.");

            swapIndices, = np.where(atomTypeNumbers == typeNumber1);

            # Check there are atoms of typeNumber1 in the structure.

            if len(swapIndices) > 0:
                indices.extend(swapIndices);

                atomTypes = atomTypes + [typeNumber2] * len(swapIndices);

        # Pass to SetAtoms() method.

        self.SetAtoms(indices, atomTypes, None);

    def DeleteAtoms(self, indexOrIndices):
        """
        Delete selected atom(s).

        Arguments:
            indexOrIndices -- an index or list of indices of atoms to delete.
        """

        if indexOrIndices == None:
            raise Exception("Error: indexOrIndices cannot be None.");

        indices = None;

        if isinstance(indexOrIndices, list):
            indices = indexOrIndices;
        else:
            indices = [indexOrIndices];

        deleteIndices = set(indices);

        if len(deleteIndices) > 0:
            atomData = self._atomData;

            numAtoms = len(atomData);
            numAtomsNew = numAtoms - len(deleteIndices);

            # Build a new array of atom data with the marked atoms removed, and update the internal fields

            atomDataNew = np.zeros(
                numAtomsNew, dtype = Structure._AtomDataType
                );

            pointer = 0;

            for i in range(0, numAtoms):
                if i not in deleteIndices:
                    atomDataNew[pointer] = atomData[i];
                    pointer += 1;

            self._atomData = atomDataNew;

            self._atomTypeNumbersView = atomDataNew.view(dtype = np.int64).reshape((numAtomsNew, 4))[:, 0];
            self._atomPositionsView = atomDataNew.view(dtype = np.float64).reshape((numAtomsNew, 4))[:, 1:];

    # --------
    # Symmetry
    # --------

    def GetSymmetryAnalysisTolerance(self):
        """ Return the symmetry-analysis tolerance. """

        if self._symmetryAnalysisTolerance != None:
            return self._symmetryAnalysisTolerance;
        else:
            # If a symmetry tolerance has not been set, return the default value that would be used when getting symmetry properties.

            return Structure._DefaultSymmetryAnalysisTolerance;

    def GetSpacegroup(self, tolerance = None):
        """
        Get the spacegroup as a (spacegroup_number, spacegroup_symbol) tuple.

        Keyword arguments:
            tolerance -- tolerance for performing the symmetry analysis.
        """

        # If the spacegroup has not already been computed, find and cache it.

        self._PerformSymmetryAnalysis(tolerance = tolerance, spacegroupOnly = True);

        return self._pSpacegroup;

    def GetSymmetryOperations(self, tolerance = None):
        """
        Get the symmetry operations as a list of (rotation, translation) tuples.

        Keyword arguments:
            tolerance -- tolerance for performing the symmetry analysis.

        Return value:
            A list of (rotation, translation) tuples.
            Rotations are stored as 3x3 NumPy integer matrices, and translations as NumPy vectors.
        """

        # Perform the symmetry analysis if required.

        self._PerformSymmetryAnalysis(tolerance = tolerance);

        return self._pSymmetryOperations;

    def GetUniqueAtomIndices(self, tolerance = None):
        """
        Get a list of representitive unique atom indices and the numbers of equivalent sites.

        Keyword arguments:
            tolerance -- tolerance for performing the symmetry analysis.

        Return value:
            A tuple containing (representitive_indices, num_equivalent_sites) lists.
        """

        # Perform the symmetry analysis if required.

        self._PerformSymmetryAnalysis(tolerance = tolerance);

        return self._pUniqueAtomIndices;

    # ----------
    # Comparison
    # ----------

    def CompareStructure(self, structure, tolerance = None):
        """
        Return True if structure has the same number of atoms, lattice vectors, atom-type numbers and positions.
        Lattice vectors and positions are compared to within a set tolerance.

        Keyword arguments:
            tolerance -- tolerance for comparing lattice vectors and positions.
        """

        # If no equivalence tolerance is supplied, use the default.

        if tolerance == None:
            tolerance = Structure.DefaultSymmetryEquivalenceTolerance;

        # Fastest check: compare numbers of atoms.

        if len(self._atomData) != structure.GetAtomCount():
            return False;

        # Next fastest check: compare lattice vectors.

        if not (np.abs(self._latticeVectors - structure.GetLatticeVectorsNumPy(copy = False)) <= tolerance).all():
            return False;

        # Next fastest check: compare atom-type numbers.

        atomTypeNumbers1 = self._atomTypeNumbersView;
        atomTypeNumbers2 = structure.GetAtomTypeNumbersNumPy(copy = False);

        if not (atomTypeNumbers1 == atomTypeNumbers2).all():
            return False;

        # Slowest check: compare atom positions.

        atomPositions1 = self._atomPositionsView;
        atomPositions2 = structure.GetAtomPositionsNumPy(copy = False);

        return (np.abs(atomPositions1 - atomPositions2) <= tolerance).all();

    # ---------------
    # Utility Methods
    # ---------------

    def GetAtomPositionsCartesian(self):
        """ Return the atom positions in Cartesian coordinates as a list of NumPy vectors. """

        if self._pAtomPositionsCartesian == None:
            v1, v2, v3 = self._latticeVectors[:];

            atomPositionsCartesian = [];

            for x, y, z in self._atomPositionsView[:]:
                atomPositionsCartesian.append(
                    x * v1 + y * v2 + z * v3
                    );

            self._pAtomPositionsCartesian = atomPositionsCartesian;

        return self._pAtomPositionsCartesian;

    def GetAtomTypeNumberPlaceholder(self):
        """ Return a negative unused atom-type number. """

        return min(0, np.min(self._atomTypeNumbersView)) - 1;

    def GetAtomTypeNumbersCounts(self):
        """ Return a tuple of (type_numbers, atom_counts) lists. """

        # Since the atom-type numbers are maintained in sort order, this is a NumPy one-liner.

        typeNumbers, atomCounts = np.unique(self._atomTypeNumbersView, return_counts = True);

        # For consistency with other methods, convert to lists when returning.

        return (list(typeNumbers), list(atomCounts));

    def GetAtomicSymbolsCounts(self, atomicSymbolLookupTable = None):
        """
        Return a tuple of (atomic_symbols, atom_counts) lists.

        Keyword arguments:
            atomicSymbolLookupTable -- a dictionary mapping atom-type numbers to custom atomic symbols.
        """

        typeNumbers, atomCounts = self.GetAtomTypeNumbersCounts();

        # Convert the atom-type numbers to symbols.

        symbols = [];

        for typeNumber in typeNumbers:
            # Try to look up the symbol from either the lookup table or the Constants module.

            symbol  = AtomTypeNumberToAtomicSymbol(typeNumber, atomicSymbolLookupTable);

            if symbol != None:
                symbols.append(symbol);
            else:
                # If this fails, convert the type number to a string.

                symbols.append(
                    str(typeNumber)
                    );

        # Return the list of symbols and counts.

        return (symbols, atomCounts);

    def GetAtomIndexRanges(self):
        """ Return a dictionary of tuples mapping atom-type numbers to index ranges. """

        atomTypeNumbers = self._atomTypeNumbersView;

        atomIndexRanges = { };

        currentStartIndex, currentAtomTypeNumber = None, None;

        for i, atomTypeNumber in enumerate(atomTypeNumbers):
            if atomTypeNumber != currentAtomTypeNumber:
                if currentAtomTypeNumber != None:
                    atomIndexRanges[currentAtomTypeNumber] = (currentStartIndex, i);

                currentStartIndex = i;
                currentAtomTypeNumber = atomTypeNumber;

        atomIndexRanges[currentAtomTypeNumber] = (currentStartIndex, len(atomTypeNumbers));

        return atomIndexRanges;

    def GetChemicalFormula(self, atomicSymbolLookupTable = None):
        """
        Build a chemical formula from the atomic composition of the structure.

        Keyword arguments:
            atomicSymbolLookupTable -- a lookup table mapping atom-type numbers to atomic symbols.
        """

        # If a lookup table is not provided and _pChemicalFormula has been initialised, return a copy of that.

        if atomicSymbolLookupTable == None and self._pChemicalFormula != None:
            return self._pChemicalFormula;

        # Build the chemical formula.

        atomicSymbols, atomCounts = self.GetAtomicSymbolsCounts(atomicSymbolLookupTable = atomicSymbolLookupTable);

        chemicalFormula = "";

        for symbol, count in zip(atomicSymbols, atomCounts):
            # If the symbol starts with a digit or a +/- character, wrap it in parentheses -- otherwise, the chemical formula will be difficult to read.

            if symbol[0] in "0123456789+-" or symbol[-1] in "0123456789+-":
                symbol = "({0})".format(symbol);

            chemicalFormula += "{0}{1}".format(symbol, count);

        # Only store a copy of the chemical formula if not using a lookup table.

        if atomicSymbolLookupTable == None:
            self._pChemicalFormula = chemicalFormula;

        return chemicalFormula;

    def GetLatticeParameters(self):
        """ Compute and return a tuple of (a, b, c, \alpha, \beta, \gamma, V) lattice parameters. """

        latticeVectors = self._latticeVectors;

        # Lattice constants are the lengths of the three lattice vectors.

        cellLengths = [
            np.linalg.norm(v) for v in latticeVectors
            ];

        # Cell angles are calculated using a . b = |a||b|cos(\theta) -> \theta = acos([a . b] / |a||b|).
        # \alpha, \beta and \gamma are the angles beteen (b, c), (a, c) and (a, b), respectively.

        cellAngles = [
            (180.0 / math.pi) * math.acos(np.dot(latticeVectors[i1], latticeVectors[i2]) / (cellLengths[i1] * cellLengths[i2]))
                for i1, i2 in [(1, 2), (0, 2), (0, 1)]
            ];

        # V = a . (b x c).

        volume = np.dot(latticeVectors[0], np.cross(latticeVectors[1], latticeVectors[2]));

        # Convert to a tuple and return.

        return tuple(cellLengths) + tuple(cellAngles) + (volume, );

    def GetNeighbourTable(self):
        """ Compute and return an NxN NumPy matrix of interatomic distances (a neighbour table). """

        if self._pNeighbourTable == None:
            atomPositions = self._atomPositionsView;

            numAtoms = len(atomPositions);

            # Build an NxN array with the vectors (in fractional coordinates) between all pairs of atoms.

            vectors = np.zeros(
                (numAtoms, numAtoms, 3), dtype = np.float64
                );

            for i in range(0, numAtoms):
                vectors[i, :] = atomPositions[:] - atomPositions[i];

            # Apply periodic boundary conditions.

            vectors[vectors < -0.5] += 1.0;
            vectors[vectors >= 0.5] -= 1.0;

            # Convert fractional to Cartesian coordinates (a NumPy one-liner...!).

            vectors = np.einsum('ijk,kl', vectors, self._latticeVectors);

            # Convert vectors to distances.

            self._pNeighbourTable = np.sqrt(np.sum(vectors ** 2, axis = 2));

        return self._pNeighbourTable;

    # --------------------
    # Manipulation Methods
    # --------------------

    def Clone(self):
        """ Return a deep copy as a new Structure object. """

        # Given the way the constructor works, all the fields except the name (which ought to be an immutable sttring anyway) should be deep copied.

        return Structure(self._latticeVectors, self._atomPositionsView, self._atomTypeNumbersView, name = self._name);

    def GetSupercell(self, dimA = 1, dimB = 1, dimC = 1):
        """ Return a dimA x dimB x dimC supercell expansion as a new Structure object. """

        if dimA < 1 or dimB < 1 or dimC < 1:
            raise Exception("Error: Supercell dimensions must be >= 1.");

        if dimA % 1 != 0 or dimB % 1 != 0 or dimC % 1 != 0:
            raise Exception("Error: Supercell dimensions must be integers.");

        # If dimA = dimB = dimC = 1, we can simply return a clone.

        if dimA == 1 and dimB == 1 and dimC == 1:
            return self.Clone();

        atomTypeNumbers = self._atomTypeNumbersView;

        numAtoms = len(atomTypeNumbers);

        # Copy and rescale atom positions.

        atomPositions = np.copy(self._atomPositionsView);

        atomPositions[:] /= np.array([dimA, dimB, dimC], dtype = np.float64);

        # Create arrays for new atom-type numbers and positions.

        newAtomTypeNumbers = np.zeros(
            dimA * dimB * dimC * numAtoms, dtype = np.int64
            );

        newAtomPositions = np.zeros(
            (len(newAtomTypeNumbers), 3), dtype = np.float64
            );

        # Generate new atom positions and atomic numbers.

        pointer = 0;

        for i in range(0, dimA):
            for j in range(0, dimB):
                for k in range(0, dimC):
                    translationVector = np.array(
                        [float(i) / dimA, float(j) / dimB, float(k) / dimC], dtype = np.float64
                        );

                    newAtomTypeNumbers[pointer:pointer + numAtoms] = atomTypeNumbers;
                    newAtomPositions[pointer:pointer + numAtoms] = atomPositions + translationVector;

                    pointer += numAtoms;

        # Generate new lattice vectors.

        newLatticeVectors = np.copy(self._latticeVectors);

        newLatticeVectors[0, :] *= dimA;
        newLatticeVectors[1, :] *= dimB;
        newLatticeVectors[2, :] *= dimC;

        # If the current structure has a name, generate a new one by appending "(<dimA>x<dimB>x<dimC> SC)" to it.

        newName = None;

        if self._name != None:
            newName = "{0} ({1}x{2}x{3} SC)".format(self._name, dimA, dimB, dimC);

        return Structure(newLatticeVectors, newAtomPositions, newAtomTypeNumbers, name = newName);

    def GetUniqueAtomicSubstitutions(self, atomType1, atomType2, tolerance = None):
        """
        Return lists of unique structures obtained by substituting atomType1 with atomType2 and associated degeneracies.
        If atomType2 is None, vacancies will be created by deleting atoms of type atomType1.

        Arguments:
            atomType1, atomType2 -- may be atom-type numbers or atomic symbols.

        Keyword arguments:
            tolerance -- tolerance to be used for identifying equivalent atoms.

        Return value:
            A tuple of (structures, degeneracies) lists containing the unique structures and associated degeneracies.
        """

        # Convert atomTypes to a atom-type numbers.

        atomTypeNumber1 = AtomTypeToAtomTypeNumber(atomType1);
        atomTypeNumber2 = AtomTypeToAtomTypeNumber(atomType2);

        atomTypeNumbers = self._atomTypeNumbersView;

        # Get the indices and site degeneracies of the unique atoms in the structure.

        uniqueAtomIndices, siteDegeneracies = self.GetUniqueAtomIndices(tolerance = tolerance);

        # Generate substituted structures.

        structures, degeneracies = [], [];

        for atomIndex, siteDegeneracy in zip(uniqueAtomIndices, siteDegeneracies):
            if atomTypeNumbers[atomIndex] == atomTypeNumber1:
                newStructure = self.Clone();

                # If atomTypeNumber2 is None, delete the atom.

                if atomTypeNumber2 != None:
                    newStructure.SetAtoms(atomIndex, atomTypeNumber2, None);
                else:
                    newStructure.DeleteAtoms(atomIndex);

                structures.append(newStructure);
                degeneracies.append(siteDegeneracy);

        return (structures, degeneracies);

    # -------------
    # Static Fields
    # -------------

    """ Structured data type used to store atom data internally. """

    _AtomDataType = [('TypeNumber', 'i8'), ('PosX', 'f8'), ('PosY', 'f8'), ('PosZ', 'f8')];

    """ Regex to parse string returned by the spglib get_spacegroup() function. """

    _SpacegroupRegex = re.compile(
        r"(?P<spacegroup_symbol>.+) \((?P<spacegroup_number>\d+)\)$"
        );

    """ Default tolerance for symmetry analysis and structure comparisons. """

    DefaultSymmetryEquivalenceTolerance = 1.0e-5;


# ---------
# Functions
# ---------

def AtomTypeToAtomTypeNumber(atomType, atomicSymbolLookupTable = None):
    """ Convert atomType, which may be an integer or an atomic symbol, to an atom-type number. """

    try:
        # Assume atomType is (convertible to) an integer.

        return int(atomType);
    except (TypeError, ValueError):
        # If that fails, assume atomType is an atomic symbol and look up the corresponding atomic number.

        if atomicSymbolLookupTable != None:
            # If a lookup table is supplied, search it for the symbol.

            for typeNumber, symbol in atomicSymbolLookupTable.items():
                if atomType == symbol:
                    return typeNumber;

        # If that fails, look up the symbol in the periodic table.

        return Constants.SymbolToAtomicNumber(atomType);

def AtomTypeNumberToAtomicSymbol(atomTypeNumber, atomicSymbolLookupTable = None):
    """ Convert an integer atomTypeNumber to an atomic symbol. """

    # If a lookup table is supplied, search for the atom-type number.

    if atomicSymbolLookupTable != None:
        if atomTypeNumber in atomicSymbolLookupTable:
            return atomicSymbolLookupTable[atomTypeNumber];

    # If not, look up the type number in the periodic table.

    return Constants.AtomicNumberToSymbol(atomTypeNumber);
