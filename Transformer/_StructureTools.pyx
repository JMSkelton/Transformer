# Transformer/_StructureTools.pyx by J. M. Skelton


# -------
# Imports
# -------

from Transformer import Structure;

import numpy as np;
cimport numpy as np;


# ------------------
# Types Declarations
# ------------------

cdef struct AtomData:
    np.int64_t type_number;
    np.float64_t x, y, z;


# ---------------------
# Function Declarations
# ---------------------

# C implementation of a comparison between a set of positions and a set of symmetry-transformed references using a set tolerance.
# Using a loop allows the function to return as soon as a match is found, which can save a significant amount of work compared to a NumPy array implementation.

cpdef bint _MergeStructureSet_ComparePositions(np.float64_t [:, :] compare_positions, np.float64_t [:, :, :] transformed_ref_positions, np.float64_t tolerance):
    cdef size_t num_sym_ops, num_atoms;

    num_atoms = compare_positions.shape[0];
    num_sym_ops = transformed_ref_positions.shape[0];

    cdef size_t i, j, k;
    cdef bint match;

    for i in range(num_sym_ops):
        match = True;

        for j in range(num_atoms):
            for k in range(0, 3):
                if abs(compare_positions[j, k] - transformed_ref_positions[i, j, k]) >= tolerance:
                    match = False;
                    break;

            if not match:
                break;

        if match:
            return True;

    return False;

# C implementation of the _MergeStructureSet_GenerateSymmetryTransformedPositions() routine from the StructureTools module.

cpdef np.float64_t [:, :, :] _MergeStructureSet_GenerateSymmetryTransformedPositions(structure, list symmetry_operations):
    # Get the atom data from the supplied Structure object.

    cdef np.ndarray[dtype = AtomData, ndim = 1] atom_data = structure.GetAtomDataNumPy(copy = False);

    cdef size_t num_atoms = atom_data.shape[0];
    cdef size_t num_sym_ops = len(symmetry_operations);

    # Typed array to hold the transformed structures.

    cdef np.ndarray[dtype = AtomData, ndim = 2] transformed_structures = np.zeros(
        (num_sym_ops, num_atoms), dtype = Structure.Structure._AtomDataType
        );

    # Loop counters.

    cdef size_t i, j;

    # Typed variables to receive the rotation matrices and translations in the list of symmetry operations.

    cdef np.ndarray[dtype = np.int32_t, ndim = 2] rotation_matrix;
    cdef np.ndarray[dtype = np.float64_t, ndim = 1] translation_vector;

    # AtomData structures to hold a local copy of the initial atom data and a working copy for applying the transformation.

    cdef AtomData atom_init, atom_trans;

    for i, (rotation_matrix, translation_vector) in enumerate(symmetry_operations):
        for j in range(num_atoms):
            atom_init = atom_data[j];

            atom_trans = atom_init;

            # Apply rotation (3x3 matrix multiplying a three-component vector).

            atom_trans.x = rotation_matrix[0, 0] * atom_init.x + rotation_matrix[0, 1] * atom_init.y + rotation_matrix[0, 2] * atom_init.z;
            atom_trans.y = rotation_matrix[1, 0] * atom_init.x + rotation_matrix[1, 1] * atom_init.y + rotation_matrix[1, 2] * atom_init.z;
            atom_trans.z = rotation_matrix[2, 0] * atom_init.x + rotation_matrix[2, 1] * atom_init.y + rotation_matrix[2, 2] * atom_init.z;

            # Apply the translation and clamp to the range [0, 1].

            atom_trans.x = (atom_trans.x + translation_vector[0]) % 1.0;
            atom_trans.y = (atom_trans.y + translation_vector[1]) % 1.0;
            atom_trans.z = (atom_trans.z + translation_vector[2]) % 1.0;

            transformed_structures[i, j] = atom_trans;

        # Sort the transformed positions.

        transformed_structures[i].sort();

    # Return a view of the transformed positions.

    return transformed_structures.view(dtype = np.float64).reshape((num_sym_ops, num_atoms, 4))[:, :, 1:];
