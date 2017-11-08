# Transformer/IO/_Common.py


# -------
# Imports
# -------

from Transformer import Constants;


# ---------
# Functions
# ---------

def AtomicSymbolsToAtomTypeNumbers(atomicSymbols, atomTypeNumberLookupTable = None):
    # If a lookup table has been supplied, ensure the keys are in title case for consistent comparison.

    if atomTypeNumberLookupTable != None:
        atomTypeNumberLookupTable = {
            key.title() : value
                for key, value in atomTypeNumberLookupTable.items()
            };

    # Convert atomic symbols to atom-type numbers.

    atomTypeNumbers = [];

    for i, symbol in enumerate(atomicSymbols):
        # Convert the symbol to title case.

        symbol = symbol.title();

        if atomTypeNumberLookupTable != None and symbol in atomTypeNumberLookupTable:
            # If a lookup table has been supplied, check this first.

            atomTypeNumbers.append(
                atomTypeNumberLookupTable[symbol]
                );
        else:
            # Next, try the periodic table in the Constants module.

            typeNumber = Constants.SymbolToAtomicNumber(symbol);

            if typeNumber != None:
                atomTypeNumbers.append(typeNumber);
            else:
                # If we still need a type number, try to convert the symbol to an integer.
                # If the conversion fails, it would be best for the user to define the symbol themselves via a lookup table -> raise an error.

                try:
                    atomTypeNumbers.append(
                        int(symbol)
                        );
                except ValueError:
                    raise Exception("Error: Unknown atomic symbol '{0}'.".format(symbol));

    # Return the converted atom-type numbers.

    return atomTypeNumbers;
