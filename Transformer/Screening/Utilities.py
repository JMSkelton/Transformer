# Transformer/Screening/Utilities.py


# -------
# Imports
# -------

import csv;
import math;

from Transformer.Utilities import IOHelper;


# ---------
# Functions
# ---------

def RankEnergies(totalEnergyGroups):
    if totalEnergyGroups == None:
        raise Exception("Error: totalEnergyGroups cannot be None.");

    energiesFlat = [];

    for key, totalEnergies in totalEnergyGroups.items():
        energiesFlat = energiesFlat + [
            (key, i, totalEnergy) for i, totalEnergy in enumerate(totalEnergies)
            ];

    # Sort by energy, then spacegroup, then structure number.

    energiesFlat.sort(key = lambda item : (item[2], item[1], item[0]));

    return energiesFlat;

def PrintRankedEnergies(rankedEnergies, energyUnits = "eV", maxPrint = None):
    # If maxPrint is set, truncate rankedEnergies if required.

    if maxPrint != None:
        if maxPrint < len(rankedEnergies):
            rankedEnergies = rankedEnergies[:maxPrint]

    # Put energy units into the third column title.

    column3Header = "E_0 [{0}]".format(energyUnits);

    # Work out fiels lengths for table formatting.

    maxStructureNumber = max(
        structureNumber for _, structureNumber, _ in rankedEnergies
        );

    field2Length = max(
        len(str(maxStructureNumber + 1)), 16
        );

    maxAbsEnergy = max(
        math.fabs(energy) for _, _, energy in rankedEnergies
        );

    field3Length = max(
        len("{0:.3f}".format(maxAbsEnergy)) + 1, len(column3Header), 16
        );

    # Print the formatted table.

    headerRowFormatCode = "{{0: ^16}} | {{1: ^{0}}} | {{2: ^{1}}}".format(field2Length, field3Length);

    headerRow = headerRowFormatCode.format("Spacegroup", "Structure #", column3Header);

    print(headerRow);
    print('-' * len(headerRow));

    dataRowFormatCode = r"{{0: <3}} {{1: <12}} | {{2: >{0}}} | {{3: >{1}.3f}}".format(field2Length, field3Length);

    for (spacegroupNumber, spacegroupSymbol), structureNumber, totalEnergy in rankedEnergies:
        print(dataRowFormatCode.format(spacegroupNumber, spacegroupSymbol, structureNumber + 1, totalEnergy));

    print("");

def ExportRankedEnergiesToCSV(rankedEnergies, filePath, energyUnits = "eV"):
    # The OpenForCSVWriter() helper routine takes care of a bug in the csv module that causes extra newline characters to be written on Windows.

    with IOHelper.OpenForCSVWriter(filePath) as outputWriter:
        outputWriterCSV = csv.writer(outputWriter, delimiter = ',', quotechar = '\"', quoting = csv.QUOTE_ALL);

        # Write header row.

        outputWriterCSV.writerow(
            ["Spacegroup #", "Spacegroup", "Structure #", "E [{0}]".format(energyUnits)]
            );

        # Output data.

        for (spacegroupNumber, spacegroupSymbol), structureNumber, totalEnergy in rankedEnergies:
            outputWriterCSV.writerow(
                [spacegroupNumber, spacegroupSymbol, structureNumber + 1, totalEnergy]
                );

def ImportRankedEnergiesFromCSV(filePath):
    rankedEnergies = [];

    with open(filePath, 'r') as inputReader:
        inputReaderCSV = csv.reader(inputReader);

        # Skip header row.

        next(inputReaderCSV);

        for row in inputReaderCSV:
            # Files exported using ExportRankedEnergiesToCSV() have rows in the format (spacegroup_number, spacegroup_symbol, structure_number, total_energy).
            # Lists returned by RankEnergies() have the format ((spacegroup_number, spacegroup_symbol), structure_number, total_energy).
            # The structure numbers in the CSV file are one-based, whereas the indices returned by RankEnergies() are zero-based.

            rankedEnergies.append(
                ((int(row[0]), row[1]), int(row[2]) - 1, float(row[3]))
                );

    return rankedEnergies;
