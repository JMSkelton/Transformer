# Transformer/_StructureSet.pyx


# -------
# Imports
# -------

import numpy as np;
cimport numpy as np;

from Transformer.Structure import Structure;


# ------------------
# Types Declarations
# ------------------

cdef struct _AtomData:
    np.int64_t type_number;
    np.float64_t x, y, z;


# ---------------------
# Function Declarations
# ---------------------

# Cython implementation of a comparison between a set of positions and a set of symmetry-transformed references using a set tolerance.
# Using a loop allows the function to return as soon as a match is found, which can save a significant amount of work compared to a NumPy array implementation.

cpdef bint _CompareAtomPositions(np.float64_t [:, :] compare_positions, np.float64_t [:, :, :] transformed_positions, np.float64_t tolerance, list compare_atom_index_ranges):
    cdef size_t num_sym_ops = transformed_positions.shape[0];

    cdef size_t i, j1, j2, j, k;
    cdef bint match;

    for i in range(num_sym_ops):
        match = True;

        for j1, j2 in compare_atom_index_ranges:
            for j in range(j1, j2):
                for k in range(0, 3):
                    if abs(compare_positions[j, k] - transformed_positions[i, j, k]) >= tolerance:
                        match = False;

                        # Breaks out of k loop.

                        break;

                if not match:
                    # Breaks out of j loop.

                    break;

            if not match:
                # Breaks out of compare_position_index_ranges enumeration.

                break;

        if match:
            return True;

    return False;

# Cython implementation of the _GenerateSymmetryTransformedPositions() routine from the Merging module.

cpdef np.float64_t [:, :, :] _GenerateSymmetryTransformedPositions(structure, list symmetry_operations, np.float64_t tolerance):
    # Get the atom data from the supplied Structure object.

    cdef np.ndarray[dtype = _AtomData, ndim = 1] atom_data = structure.GetAtomDataNumPy(copy = False);

    cdef size_t num_atoms = atom_data.shape[0];
    cdef size_t num_sym_ops = len(symmetry_operations);

    # Typed array to hold the transformed structures.

    cdef np.ndarray[dtype = _AtomData, ndim = 2] transformed_structures = np.zeros(
        (num_sym_ops, num_atoms), dtype = Structure._AtomDataType
        );

    # Loop counters.

    cdef size_t i, j;

    # Typed variables to receive the rotation matrices and translations in the list of symmetry operations.

    cdef np.ndarray[dtype = np.int32_t, ndim = 2] rotation_matrix;
    cdef np.ndarray[dtype = np.float64_t, ndim = 1] translation_vector;

    # AtomData structures to hold a local copy of the initial atom data and a working copy for applying the transformation.

    cdef _AtomData atom_init, atom_trans;

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

            # Workaround for fractional coordinates not wrapped to 0.0 by the % operator.

            if abs(atom_trans.x - 1.0) < tolerance:
                atom_trans.x = 0.0;

            if abs(atom_trans.y - 1.0) < tolerance:
                atom_trans.y = 0.0;

            if abs(atom_trans.z - 1.0) < tolerance:
                atom_trans.z = 0.0;

            transformed_structures[i, j] = atom_trans;

        # Sort the transformed positions.

        transformed_structures[i].sort();

    # Return a view of the transformed positions.

    return transformed_structures.view(dtype = np.float64).reshape((num_sym_ops, num_atoms, 4))[:, :, 1:];
