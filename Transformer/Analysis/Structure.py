# Transformer/Analysis/Structure.def


# ----------------
# Module Docstring
# ----------------

""" Contains routines for analysing structures. """


# -------
# Imports
# -------

import math;
import warnings;

import numpy as np;

from Transformer import Structure;


# ---------
# Functions
# ---------

def Gaussian(x, a, mu, sigma):
    """
    Return G(x) = (a / (\sigma * sqrt(2\pi))) * exp(-(x - \mu) ^ 2 / (2 \sigma ^ 2)).
    """

    return (a / (sigma * math.sqrt(2.0 * math.pi))) * np.exp(-1.0 * (x - mu) ** 2 / (2.0 * sigma ** 2));

def RadialDistributionFunction(structure, supercell = None, atomType1 = None, atomType2 = None, rMax = None, binWidth = 0.01, broadening = None):
    """
    Calculate the radial distribution function g(r) for the supplied structure.

    Keyword arguments:
        supercell -- calculat RDF using a (dim1, dim2, dim3) supercell expansion of structure.
        atomType1, atomType2 -- specify atom type(s) to calculate partial RDFs (default: None = all atoms/total RDF).
        rMax -- radius to calculate RDF to (default: half the length of the shortest lattice vector in the structure/supercell).
        binWidth -- resolution of interatomic-distance histogram (default: 0.01).
        broadening -- if supplied, specifies the width of a Gaussian function to broaden the histogram after binning.

    Return value:
        A tuple of (r, g_r, intg_r) giving r, g(r) and its numerical integral intg(r).

    Notes:
        If rMax is over half the length of the shortest lattice vector in the structure/supercell, some pair distances may be missing from the histogram; this can be rectified by specifying a (larger) supercell expansion.
        During normalisation, the r = 0 bin is set to zero; binWidth should therefore be small enough to ensure that pair distances are not placed into this bin.
    """

    # If supercell is set, generate a supercell expansion of the structure.

    if supercell is not None:
        structure = structure.GetSupercell(*supercell);

    # If supplied, check atom types and convert to index ranges for specifying pair distances.

    numAtoms = structure.GetAtomCount();
    atomIndexRanges = structure.GetAtomIndexRanges();

    # To avoid copy/pasting code (!).

    atomIndices = [None, None];

    for i, atomType in enumerate([atomType1, atomType2]):
        if atomType is not None:
            atomTypeNumber = Structure.AtomTypeToAtomTypeNumber(atomType);

            if atomTypeNumber is None:
                # atomType is not an integer and does not map to an atomic number.

                raise Exception("Error: Unable to convert atomic symbol {0} to an atom-type number.".format(atomType));

            if atomTypeNumber not in atomIndexRanges:
                raise Exception("Error: Atom type {0} was not found in the supplied structure.".format(atomType));

            i1, i2 = atomIndexRanges[atomTypeNumber];

            atomIndices[i] = np.array(
                [i for i in range(i1, i2)], dtype = np.int
                );
        else:
            # If an atom-type number is not specified, default to all atoms.

            atomIndices[i] = np.array(
                [i for i in range(0, numAtoms)], dtype = np.int
                );

    atomIndices1, atomIndices2 = atomIndices;

    numAtoms1, numAtoms2 = len(atomIndices1), len(atomIndices2);

    # The maximum distance over which we can include all unique interatomic distances is half the shortest lattice vector.

    maxDist = min(
        np.linalg.norm(v) / 2.0 for v in structure.GetLatticeVectors()
        );

    if rMax is not None:
        if rMax <= 0.0:
            raise Exception("Error: If supplied, rMax must be greater than 0.");

        if rMax > maxDist:
            warnings.warn("For r > min(|a|,|b|,|c|) / 2 the RDF may not include all the unique interatomic distances.", UserWarning);
    else:
        rMax = maxDist;

    # Bin interatomic distances into a histogram.

    if binWidth >= rMax:
        raise Exception("Error: binWidth must be less than rMax.");

    neighbourTable = structure.GetNeighbourTable();

    binEdges = np.arange(
        0.0, rMax + 1.0e-5, binWidth
        );

    # Mask distances in the neighbour table to be binned.

    maskIndices1 = np.zeros((numAtoms, numAtoms), np.bool);
    maskIndices1[atomIndices1, :] = True;

    maskIndices2 = np.zeros((numAtoms, numAtoms), np.bool);
    maskIndices2[:, atomIndices2] = True;

    mask = np.logical_and(maskIndices1, maskIndices2);

    # Unset the diagonal elements ("self distances").

    mask[np.eye(numAtoms, dtype = np.bool)] = False;

    # If the bin width is too large, some interatomic distances may be counted in the r = 0 bin, which is assumed to be zero during normalisation.
    # We therefore check bin width against the masked neighbour table, and issue a warning if this is the case.

    if (neighbourTable[mask] <= binWidth).any():
        warnings.warn("One or more interatomic distances are smaller than the histogram bin width -- check the input structure and recalculate with a finer resolution if required.", UserWarning);

    bins, binEdges = np.histogram(
        neighbourTable[mask], bins = binEdges
        );

    # g(r) is defined as the number of atoms in a shell of r + dr centered on a reference atom; the r axis therefore corresponds to the left edges of the bins.

    rdfR = binEdges[:-1];

    # Calculate \rho = N_B / V.

    _, _, _, _, _, _, v = structure.GetLatticeParameters();

    rho = numAtoms2 / v;

    # To obtain g(r) we normalise the histogram by \rho * N_A * 4 \pi r ^ 2 dr.
    # g(r = 0) _should_ be zero -> to avoid a divide-by-zero error we simply leave the first bin to zero.

    rdfGR = np.zeros_like(rdfR);

    rdfGR[1:] = bins[1:] / (rho * numAtoms1 * 4.0 * math.pi * rdfR[1:] ** 2 * binWidth);

    # Apply a broadening if required.

    if broadening is not None:
        if broadening <= 0.0 or 10.0 * broadening >= rMax:
            raise Exception("Error: If supplied, broadening should be > 0 and < 1/10 rMax.");

        convX = np.arange(-5.0 * broadening, 5.0 * broadening + 1.0e-5, binWidth);

        convY = Gaussian(convX, 1.0, 0.0, broadening);
        convY /= convY.sum();

        rdfGR = np.convolve(rdfGR, convY, mode = 'same');

        # Sanity check.

        assert len(rdfGR) == len(rdfR);

    # Calculate the cumulative numerical integral.

    rdfIntGR = np.cumsum(
        rdfGR * rho * 4.0 * math.pi * rdfR ** 2 * binWidth
        );

    # Return the RDF and its integral.

    return (rdfR, rdfGR, rdfIntGR);
