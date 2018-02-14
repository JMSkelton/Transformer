# Transformer/Framework/FilterBases.py


# -------
# Imports
# -------

import warnings;

from Transformer import Structure;

# Try to import the tqdm module to display progress bars in some of the slower routines.

_TQDM = False;

try:
    import tqdm;

    _TQDM = True;
except ImportError:
    pass;


# ----------------
# FilterBase Class
# ----------------

class FilterBase(object):
    # -----------
    # Constructor
    # -----------

    def __init__(self):
        # Initialise fields.

        self._substitutions = None;

        self._tolerance = None;

        self._useMP = False;
        self._mpNumProcesses = None;

        self._substitutionIndex = None;

    # -------
    # Methods
    # -------

    def IsMPSafe(self):
        # The safe default is to return False; however, since this really _should_ be set in a derived class, we issue a warning if this method is called.

        warnings.warn("IsMPSafe() should be overridden in derived classes.", RuntimeWarning);

        return False

    def RequiresFilteredStructures(self):
        # The safe default is to return True; again, this should be set in a derived class, so we issue a warning if this method gets called.

        warnings.warn("RequiresFilteredStructures() should be overridden in derived classes.", RuntimeWarning);

        return True;

    def Initialise(self, substitutions, tolerance, printProgressUpdate, useMP, mpNumProcesses):
        # Update fields.

        self._substitutions = substitutions;

        self._tolerance = tolerance;

        self._printProgressUpdate = printProgressUpdate;

        self._useMP = useMP;
        self._mpNumProcesses = mpNumProcesses;

    def OnStartSubstitution(self, index):
        # Update the _substitutionIndex field.

        self._substitutionIndex = index;

    def TestSubstitutedStructure(self, structure, degeneracy):
        # Accept all structures (i.e. do nothing).

        return True;

    def SetFilteredStructureSet(filteredStructureSet, degeneracies):
        # Do nothing.

        pass;

    def FinaliseMergedStructureSet(self, structureSet):
        # Do nothing.

        pass;


# ------------------------
# CoverageFilterBase Class
# ------------------------

class CoverageFilterBase(FilterBase):
    # -----------
    # Constructor
    # -----------

    def __init__(self, coverage = 'full'):
        # Call the base class constructor.

        super(CoverageFilterBase, self).__init__();

        # Parameter validation.

        if coverage not in CoverageFilterBase.CoverageLevels:
            raise Exception(
                "Error: coverage must be one of {0}.".format(", ".join(CoverageFilterBase.CoverageLevels))
                );

        # Store the coverage level.

        self._coverage = coverage;

        # Initialise empty fields for storing the current and last sets of filtered structures passed to SetFilteredStructures().

        self._lastFilteredStructureSet = None;
        self._currentFilteredStructureSet = None;

    # -------
    # Methods
    # -------

    def RequiresFilteredStructures(self):
        # If the coverage level is set to 'low', or the final substitution is being processed, the filtered structures are not required.

        return self._coverage != 'low' and self._substitutionIndex < len(self._substitutions) - 1;

    def OnStartSubstitution(self, index):
        # If the coverage level is not set to 'low', update the last filtered structure set and call the base class method.

        if self._coverage != 'low':
            self._lastFilteredStructureSet = self._currentFilteredStructureSet;
            self._currentFilteredStructureSet = None;

            self._substitutionIndex = index;

    def SetFilteredStructureSet(self, filteredStructureSet):
        # If the coverage level is not set to 'low', store the filtered structures and degeneracies.

        if self._coverage != 'low':
            self._currentFilteredStructureSet = filteredStructureSet;

    def FinaliseMergedStructureSet(self, structureSet):
        # For all bar the 'low' coverage level, we need to apply the current substitution to the last set of filtered structures (if available).

        if self._coverage != 'low':
            lastFilteredStructureSet = self._lastFilteredStructureSet;

            if lastFilteredStructureSet != None and lastFilteredStructureSet.GetStructureCount() > 0:
                structures, degeneracies = self._lastFilteredStructureSet.GetStructureSetFlat();

                printProgressUpdate = self._printProgressUpdate;

                if printProgressUpdate:
                    # Since the user will lilely be using a derived class, we get the class name for the message (rather than printing "CoverageFilterBase").

                    print("{0}: Finalising merged structure set for {1} coverage.".format(self.__class__.__name__, self._coverage));

                # Decide whether to display a tqdm progress bar while generating new child structures.

                progressBar = printProgressUpdate and _TQDM;

                # Load the current substitution and convert atom types to atom-type numbers.

                atomType1, atomType2 = self._substitutions[self._substitutionIndex];

                atomTypeNumber1 = Structure.AtomTypeToAtomTypeNumber(atomType1);
                atomTypeNumber2 = Structure.AtomTypeToAtomTypeNumber(atomType2);

                # Load the symmetry tolerance.

                tolerance = self._tolerance;

                # Load the current filtered structure set if required.

                currentFilteredStructureSet = None;

                if self._coverage == 'full':
                    currentFilteredStructureSet = self._currentFilteredStructureSet;

                # Count the number of structures added to the merged and filtered structure sets.

                addCount, addCountFiltered = 0, 0;

                if progressBar:
                    print("");

                # Set up a primary iterator; if progressBar is set, wrap the iterator in a tqdm progress bar.

                iValues = range(0, len(structures));

                if _TQDM:
                    iValues = tqdm.tqdm(iValues);

                for i in iValues:
                    structure, degeneracy = structures[i], degeneracies[i];

                    # Generate substituted child structures.

                    newStructures, newDegeneracies = structure.GetUniqueAtomicSubstitutions(atomTypeNumber1, atomTypeNumber2, tolerance = tolerance);

                    # The degeneracies obtained from the Structure.GetUniqueAtomicSubstitutions() function need to be multiplied by that of the parent structure.

                    newDegeneracies = [
                        newDegeneracy * degeneracy for newDegeneracy in newDegeneracies
                        ];

                    # Test each structure with the TestSubstitutedStructure() method.

                    for newStructure, newDegeneracy in zip(newStructures, newDegeneracies):
                        if self.TestSubstitutedStructure(newStructure, newDegeneracy):
                            if structureSet.Add(newStructure, newDegeneracy):
                                addCount += 1;

                        elif currentFilteredStructureSet != None:
                            if currentFilteredStructureSet.Add(newStructure, newDegeneracy):
                                addCountFiltered += 1;

                if progressBar:
                    print("");

                if printProgressUpdate:
                    if addCount > 0:
                        print("{0}: Added {1} structure(s) to merged structure set.".format(self.__class__.__name__, addCount));

                    if addCountFiltered > 0:
                        print("{0}: Added {1} structure(s) to the filtered list.".format(self.__class__.__name__, addCount));

                    if addCount > 0 or addCountFiltered > 0:
                        print("");

    # -------------
    # Static Fields
    # -------------

    CoverageLevels = ['low', 'medium', 'full'];


# -----------------------
# RankingFilterBase Class
# -----------------------

class RankingFilterBase(CoverageFilterBase):
    # -----------
    # Constructor
    # -----------

    def __init__(self, coverage = 'full', cutoffMode = 'percentage', cutoff = 100):
        # Call the base class constructor.

        super(RankingFilterBase, self).__init__(coverage = coverage);

        # Validate and store cutoff settings.

        cutoffMode = cutoffMode.lower();

        if cutoffMode not in ['percentage', 'number']:
            raise Exception("Error: Unknown cutoff mode '{0}'.".format(cutoffMode));

        if cutoffMode == 'percentage' and (cutoff <= 0.0 or cutoff > 100.0):
                raise Exception("Error: For the 'percentage' ranking mode, cutoff should be > 0 and < 100.");
        elif cutoffMode == 'number':
            cutoff = int(cutoff);

            if cutoff <= 0:
                raise Exception("Error: For the 'number' ranking mode, cutoff should be > 0.");

        self._cutoffMode = cutoffMode;
        self._cutoff = cutoff;

    # -------
    # Methods
    # -------

    def IsMPSafe(self):
        return True;

    def CalculateScores(self, structureSet):
        raise Exception("Error: CalculateScores() must be overridden in derived classes.");

    def FinaliseMergedStructureSet(self, structureSet):
        # Call the base class method to update the structure set according to the coverage level (if required).

        super(RankingFilterBase, self).FinaliseMergedStructureSet(structureSet);

        # Pass the structure set to the virtual CalculateScores() method to score the structures.

        scores = self.CalculateScores(structureSet);

        # Rank the structures.

        scoresFlat = [];

        for spacegroup in scores.keys():
            for i, score in enumerate(scores[spacegroup]):
                scoresFlat.append(
                    (spacegroup, i, score)
                    );

        scoresFlat.sort(key = lambda item : item[2]);

        # Obtain the index at which to split the ranked structures into keep and remove sets.

        splitIndex = None;

        cutoffMode, cutoff = self._cutoffMode, self._cutoff;

        if cutoffMode == 'percentage':
            splitIndex = round(
                (cutoff / 100.0) * len(scoresFlat)
                );

            splitIndex = max(int(splitIndex), 1);

        elif cutoffMode == 'number':
            splitIndex = min(cutoff, len(scoresFlat));

        # Collect a dictionary of structures to remove.

        removeIndices = { };

        for spacegroup, structureNumber, _ in scoresFlat[splitIndex:]:
            if spacegroup not in removeIndices:
                removeIndices[spacegroup] = [];

            removeIndices[spacegroup].append(structureNumber);

        # If the coverage setting requires filtered structures, copy them into a new structure set and pass it to the base class SetFilteredStructureSet() method.

        if self.RequiresFilteredStructures():
            structureSetDict = structureSet.GetStructureSet();

            filteredStructures, filteredDegeneracies = [], [];

            for spacegroup, indices in removeIndices.items():
                structures, degeneracies = structureSetDict[spacegroup];

                filteredStructures = filteredStructures + [structures[index] for index in indices];
                filteredDegeneracies = filteredDegeneracies + [degeneracies[index] for index in indices];

            removedStructures = structureSet.CloneNew(
                structures = filteredStructures, degeneracies = filteredDegeneracies, noInitialMerge = True
                );

            self.SetFilteredStructureSet(removedStructures);

        # Remove the marked structures from the current set.

        structureSet.Remove(removeIndices);

    # -------------
    # Static Fields
    # -------------

    CutoffModes = ['percentage', 'number'];
