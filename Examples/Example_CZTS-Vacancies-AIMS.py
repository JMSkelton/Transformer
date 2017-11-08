# Example_CZTS-Vacancies-AIMS.py by J. M. Skelton


# -------
# Imports
# -------

# Import routines from Transformer.

from Transformer.IO.StructureIO import ReadStructure;

from Transformer.Framework.BatchIO import ExportResultSet;
from Transformer.Framework.Core import AtomicSubstitutions;


# ----
# Main
# ----

# Read the input structure.

structure = ReadStructure("Cu2ZnSnS4.geometry.in");

# Generate a 2x2x1 supercell (64 atoms).

supercell = structure.GetSupercell(2, 2, 1);

# Generate structures with 1-2 Cu vacancies.

# In Transformer, vacancy creation is simply a substitution problem where atoms are substituted with None (i.e. nothing).

substitutions = [('Cu', None)] * 2;

# Use the AtomicSubstitutions convenience function to generate the defective structures.

resultSet = AtomicSubstitutions(
    supercell, substitutions
    );

# Export the results.

ExportResultSet(
    resultSet, prefix = "CZTS-Cu-Vacancies", archiveDirectory = r"Example_CZTS-Vacancies", fileFormat = 'aims'
    );

# Generate structures with a single Zn or Sn defect.

for vacancyAtom in 'Zn', 'Sn':
    substitutions = [(vacancyAtom, None)];

    resultSet = AtomicSubstitutions(
        supercell, substitutions
        );

    ExportResultSet(
        resultSet, prefix = r"CZTS-{0}-Vacancy".format(vacancyAtom), archiveDirectory = r"Example_CZTS-Vacancies", fileFormat = 'aims'
        );

# Generate structures with up to four S vacancies.

substitutions = [('S', None)] * 4;

resultSet = AtomicSubstitutions(
    supercell, substitutions
    );

ExportResultSet(
    resultSet, prefix = r"CZTS-S-Vacancies", archiveDirectory = r"Example_CZTS-Vacancies", fileFormat = 'aims'
    );

# Generate structures with Schottky defects (i.e. cation vacancies with balanced S vacancies to preserve charge neutrality).

for cation, numCationVacancies, numAnionVacancies in ('Cu', 2, 1), ('Zn', 1, 1), ('Sn', 1, 2):
    substitutions = [(cation, None)] * numCationVacancies;
    substitutions = substitutions + [('S', None)] * numAnionVacancies;

    # Wrapping the substitution sequence in a list will cause AtomicSubstitutions to treat it as a single operation with respect to generting output structures.

    resultSet = AtomicSubstitutions(
        supercell, [substitutions]
        );

    ExportResultSet(
        resultSet, prefix = r"CZTS-{0}-SchottkyDefect".format(cation), archiveDirectory = r"Example_CZTS-Vacancies", fileFormat = 'aims'
        );
