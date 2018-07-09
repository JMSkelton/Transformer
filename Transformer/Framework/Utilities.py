# Transfomer/Framework/Utilities.py


# ----------------
# Module Docstring
# ----------------

""" Utilities for (re)processing sets of structures in parallel. """


# -------
# Imports
# -------

import multiprocessing;

from Transformer import StructureSet;

from Transformer.Utilities import MultiprocessingHelper;


# --------------------------
# MergeStructureSets Routine
# --------------------------

def MergeStructureSets(structureSets, useMP = False, mpNumProcesses = None, printProgressUpdate = True, inPlace = False):
    """
    Merge the list of structure sets in structureSets using the StructureSet.UpdateUnion() method.

    Arguments:
        structureSets -- list of StructureSet objects to merge.

    Keyword arguments:
        useMP -- passed to StructureSet.Union().
        mpNumProcesses -- passed to StructureSet.Union().
        printProgressUpdate -- if True (default), print status messages during merging.
        inPlace -- if True, this routine will work directly on the StructureSet objects in structureSets, which will cause the first set in structureSets to be modified (default: False).

    Return value:
        StructureSet object obtained after merging the items in structureSets.

    Notes:
        inPlace = True is a performance optimisation and should be used only when the supplied structure sets are no longer needed.
        If inPlace is not set (the default), the structure sets are cloned before merging using the StructureSet.CloneNew() routine.
    """

    # If inPlace is not set, clone the structure sets.

    if not inPlace:
        structureSetsNew = [];

        for structureSet in structureSets:
            structures, degeneracies = structureSet.GetStructureSetFlat();

            structureSetsNew.append(
                structureSet.CloneNew(structures = structures, degeneracies = degeneracies, noInitialMerge = True)
                );

        structureSets = structureSetsNew;

    # If there is only one structure set, no need to do anything.

    if len(structureSets) == 1:
        return structureSets[0];

    # Setting up a tqdm-based progress bar for the reduction would be (a) fiddly, and (b) not particularly informative.
    # If a progress bar is requested, we print a set of status messages instead.

    if printProgressUpdate:
        numStructures = sum(
            structureSet.GetStructureCount() for structureSet in structureSets
            );

        print("MergeStructureSets: Merging {0} structure sets w/ {1} structure(s)".format(len(structureSets), numStructures));

    numStructures = sum(
        structureSet.GetStructureCount()
            for structureSet in structureSets
        );

    # Format string for printing status messages.

    formatString = None;

    if printProgressUpdate:
        formatString = "{{0: >{0},}}".format(len("{0:,}".format(numStructures)));

    # Perform a sequence of unions to merge the structure sets.

    structureSetRef = structureSets[0];

    numStructuresMerged = structureSetRef.GetStructureCount();

    for structureSetAdd in structureSets[1:]:
        structureSetRef.UpdateUnion(
            structureSetAdd, useMP = useMP, mpNumProcesses = mpNumProcesses
            );

        numStructuresMerged += structureSetAdd.GetStructureCount();

        if printProgressUpdate:
            numRemaining = numStructures - numStructuresMerged;

            statusMessage = None;

            if numRemaining > 0:
                statusMessage = "MergeStructureSets: Merged {0} -> {1} structure, {2} remaining".format(
                    formatString.format(numStructuresMerged), formatString.format(structureSetRef.GetStructureCount()), formatString.format(numRemaining)
                    );
            else:
                statusMessage = "MergeStructureSets: Merged {0} -> {1} structure".format(
                    formatString.format(numStructuresMerged), formatString.format(structureSetRef.GetStructureCount())
                    );

            print(statusMessage);

    if printProgressUpdate:
        print("");

    return structureSetRef;

def _MergeStructureSets_MapFunction(args):
    """ Method run by worker processes when performing parallel merging. """

    # Unpack arguments.

    structureSet1, structureSet2 = args;

    # Merge the second structure set into the first.

    structureSet1.UpdateUnion(structureSet2);

    # Clear symmetry-expansions cache before returning.

    structureSet1.ClearSymmetryExpansionsCache();

    # Return the merged structure set.

    return structureSet1;


# ---------------------------
# ReduceStructureList Routine
# ---------------------------

def ReduceStructureList(structures, degeneracies = None, useMP = False, mpNumProcesses = None, printProgressUpdate = True, **kwargs):
    """
    Reduce the list of structures and, optionally, degeneracies, by merging them into StructureSet objects set up using the supplied keyword arguments.

    Arguments:
        structures -- list of structures to reduce.

    Keyword arguments:
        degeneracies -- optional list of degeneracies for each structure in structures (default: None).
        useMP -- if True, perform parts of the reduction in parallel using process-based multithreading (default: False).
        mpNumProcesses -- if useMP is set, specifies the number of processes to use for reduction (default: automatically determined from the number of CPU cores and the list of structures).
        printProgressUpdate -- if True, print status messages.
        **kwargs -- passed through to the StructureSet constructor.

    Return value:
        A StructureSet object containing the reduced structure set and associated degeneracies.

    Notes:
        This method provides a convenience function for reducing a list of structures in parallel using similar process-based threading to the AtomicSubstitutions() routine and its derivatives.
        Its main purpose is to allow a structure set built using one symmetryExpansion setting to be "re-reduced" with a stricter one (e.g. symmetryExpansion = 'fast' -> 'full').
        The serial code path selected by setting useMP = False simply calls the StructureSet class constructor and is provided for API consistency.
    """

    if structures == None or len(structures) == 0:
        raise Exception("Error: structures cannot be None and must contain at least one structure.")

    if degeneracies != None and len(degeneracies) != len(structures):
        raise Exception("Error: If supplied, degeneracies must be the same length as structures.");

    # If there is only one structure in the input list, there's no point in taking the parallel code path.

    if len(structures) == 1:
        useMP = False;

    if useMP:
        # If not supplied, set mpNumProcesses.

        if mpNumProcesses == None:
            mpNumProcesses = min(
                len(structures), MultiprocessingHelper.CPUCount()
                );

        # Set up one _StructureSetAccumulator object per worker process.

        accumulators = [
            _StructureSetAccumulator(**kwargs)
                for i in range(0, mpNumProcesses)
            ];

        # Use the MultiprocessingHelper.QueueAccumulate() routine to merge the structures into structure sets.

        structureSets = MultiprocessingHelper.QueueAccumulate(
            [item for item in zip(structures, degeneracies)], accumulators, progressBar = printProgressUpdate
            );

        # Merge the structure sets using the MergeStructureSets() utility routine.

        return MergeStructureSets(
            structureSets, useMP = True, mpNumProcesses = mpNumProcesses, printProgressUpdate = printProgressUpdate, inPlace = True
            );

    else:
        # Serial code path.

        return StructureSet.StructureSet(
            structures = structures, degeneracies = degeneracies, **kwargs
            );

# Helper class for grouping a list of structures into sets.

class _StructureSetAccumulator(MultiprocessingHelper.AccumulatorBase):
    """
    Implementation of the MultiprocessingHelper.AccumulatorBase class for building structure sets in parallel.
    This class is used as part of the parallel code path in the ReduceStructureList() routine.
    """

    def __init__(self, **kwargs):
        """
        Class constructor.

        Keyword arguments:
            **kwargs are passed to the StructureSet constructor when setting up the internal structure set.
        """

        # Initialise an internal StructureSet object with the supplied kwargs.

        self._structureSet = StructureSet.StructureSet(
            **kwargs
            );

    def Accumulate(self, item):
        """ Implementation of the AccumulatorBase.Accumulate() method. """

        # Unpack input item and add the structure and degeneracy to the internal structure set.

        structure, degeneracy = item;

        self._structureSet.Add(
            structure, degeneracy = degeneracy
            );

    def Finalise(self):
        """ Implementtion of the AccumulatorBase.Finalise() method. """

        structureSet = self._structureSet;

        # If using, clear the symmetry-expansions cache to save memory.

        structureSet.ClearSymmetryExpansionsCache();

        # Return the internal structure set.

        return self._structureSet;
