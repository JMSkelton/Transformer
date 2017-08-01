# Transformer/DevelopmentTools.py by J. M. Skelton


# -------
# Imports
# -------

from Transformer import StructureTools;
from Transformer import ConvenienceFunctions;


# ---------
# Functions
# ---------

def MapResultSetStructures(parentStructure, structures1, structures2, tolerance = None):
    # Take the symmetry operations from the parent structure.

    symmetryOperations = parentStructure.GetSymmetryOperations();

    # Map structures in the first result set to those in the second.

    structureMappings = [];

    for i, structureRef in enumerate(structures1):
        # Expand the reference structure to a set of comparison structures by applying the parent symmetry operations.

        compareStructures = [
            structureRef.GetSymmetryTransform(rotation, translation)
                for rotation, translation in symmetryOperations
            ];

        # Remove duplicates to avoid unnecessary comparisons.

        compareStructures, _ = StructureTools.MergeStructureSet(
            compareStructures, tolerance = tolerance
            );

        # Compare each transformed structure to each structure in the second result set.

        indices = [];

        for j, structure in enumerate(structures2):
            for compareStructure in compareStructures:
                if compareStructure.CompareCell(structure):
                    indices.append(j);
                    break;

        structureMappings.append(indices);

    return structureMappings;
