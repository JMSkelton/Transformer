# Transformer/Framework/Filters.py


# -------
# Imports
# -------

from Transformer.Framework.FilterBases import CoverageFilterBase, RankingFilterBase;


# ----------------------
# SpacegroupFilter Class
# ----------------------

class SpacegroupFilter(CoverageFilterBase):
    # -----------
    # Constructor
    # -----------

    def __init__(self, removeSpacegroups, coverage = 'full'):
        # Parameter validation.

        if removeSpacegroups == None:
            raise Exception("removeSpacegroups must not be None.")

        # Call the base class constructor.

        super(SpacegroupFilter, self).__init__(coverage = coverage);

        # Store the list of spacegroups to remove.

        self._removeSpacegroups = set(
            int(spacegroup) for spacegroup in removeSpacegroups
            );

    # -------
    # Methods
    # -------

    def IsMPSafe(self):
        return True;

    def TestSubstitutedStructure(self, structure, degeneracy):
        # Get the spacegroup number for the supplied structure.

        spacegroupNumber, _ = structure.GetSpacegroup(tolerance = self._tolerance);

        # Return False for (i.e. reject) structures with spacegroups in the list of groups to remove.

        return spacegroupNumber not in self._removeSpacegroups;


# ------------------
# EnergyFilter Class
# ------------------

class EnergyFilter(RankingFilterBase):
    def __init__(self, energyCalculator, coverage = 'full', cutoffMode = 'percentage', cutoff = 100):
        # Parameter validation.

        if energyCalculator == None:
            raise Exception("Error: energyCalculator must not be None.");

        # Call the base class constructor.

        super(EnergyFilter, self).__init__(coverage = coverage, cutoffMode = cutoffMode, cutoff = cutoff);

        # Store the energy calculator.

        self._energyCalculator = energyCalculator;

    def CalculateScores(self, structureSet):
        # Pass the structure set to the CalculateTotalEnergiesGrouped() method of the energy calculator to calculate the energies.

        return self._energyCalculator.CalculateTotalEnergiesGrouped(
            structureSet.GetStructureSet(), raiseOnError = True, progressBar = self._printProgressUpdate
            );
