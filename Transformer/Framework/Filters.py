# Transformer/Framework/Filters.py


# -------
# Imports
# -------

from Transformer.Framework.FilterBases import CoverageFilterBase;


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

        # Store the list of spacegroups to remove and the coverage level.

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
