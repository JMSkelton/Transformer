# Example_CZTS-Vacancies.py by J. M. Skelton


# -------
# Imports
# -------

# Import routines from Transformer.

from Transformer.IO import ReadPOSCARFile;
from Transformer.ConvenienceFunctions import AtomicSubstitutions, ExportAtomicSubstitutionResultSet;


# ----
# Main
# ----

# Read the input structure.

structure = ReadPOSCARFile("Cu2ZnSnS4.vasp");

# Generate a 2x2x1 supercell (64 atoms).

supercell = structure.GetSupercell((2, 2, 1));

# Generate structures with 1-2 Cu vacancies.

# As far as Transformer is concerned, this is an atomic-substitution problem in which the atom being substituted is being replaced by None (nothing).

substitutions = [('Cu', None)] * 2;

# Use the AtomicSubstitutions convenience function to generate the defective structures.

_, resultSet = AtomicSubstitutions(
    supercell, substitutions
    );

# Export the results.

ExportAtomicSubstitutionResultSet(
    resultSet, prefix = "CZTS-Cu-Vacancies", workingDirectory = r"Example_CZTS-Vacancies"
    );

# Generate structures with a single Zn or Sn defect.

for vacancyAtom in 'Zn', 'Sn':
    substitutions = [(vacancyAtom, None)];

    _, resultSet = AtomicSubstitutions(
        supercell, substitutions
        );

    ExportAtomicSubstitutionResultSet(
        resultSet, prefix = r"CZTS-{0}-Vacancy".format(vacancyAtom), workingDirectory = r"Example_CZTS-Vacancies"
        );

# Generate structures with 1-4 S vacancies.
# The final round of substitutions (4th S vacancy) takes a few minutes to run; the example can be sped up by lowering the number of substitutions.

substitutions = [('S', None)] * 4;

_, resultSet = AtomicSubstitutions(
    supercell, substitutions
    );

ExportAtomicSubstitutionResultSet(
    resultSet, prefix = r"CZTS-S-Vacancies", workingDirectory = r"Example_CZTS-Vacancies"
    );

# Generate structures with Schottky defects (i.e. cation vacancies with balanced S vacancies to preserve charge neutrality).

for cation, numCationVacancies, numAnionVacancies in ('Cu', 2, 1), ('Zn', 1, 1), ('Sn', 1, 2):
    substitutions = [(cation, None)] * numCationVacancies;
    substitutions = substitutions + [('S', None)] * numAnionVacancies;

    # When performing the substitutions, use the storeIntermediate argument to retain only the initial and final sets of structures.

    _, resultSet = AtomicSubstitutions(
        supercell, substitutions, storeIntermediate = [0, len(substitutions)]
        );

    ExportAtomicSubstitutionResultSet(
        resultSet, prefix = r"CZTS-{0}-SchottkyDefect".format(cation), workingDirectory = r"Example_CZTS-Vacancies"
        );
