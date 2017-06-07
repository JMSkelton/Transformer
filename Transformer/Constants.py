# Transformer/Constants.py by J. M. Skelton


# --------------
# Periodic Table
# --------------

# The symbol 'X' is assigned to zero -- spglib functions will happily accept this as a wildcard.

PeriodicTable = [
     'X',
     'H', 'He',
    'Li', 'Be',                                                                                                                                                  'B',  'C',  'N',  'O',  'F', 'Ne',
    'Na', 'Mg',                                                                                                                                                 'Al', 'Si',  'P',  'S', 'Cl', 'Ar',
     'K', 'Ca', 'Sc',                                                                                     'Ti',  'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr',
    'Rb', 'Sr',  'Y',                                                                                     'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te',  'I', 'Xe',
    'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta',  'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn',
    'Fr', 'Ra', 'Ac', 'Th', 'Pa',  'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds', 'Rg', 'Cn', 'Nh', 'Fl', 'Mc', 'Lv', 'Ts', 'Og'
    ];

def SymbolToAtomicNumber(symbol):
    for i, testSymbol in enumerate(PeriodicTable):
        # Make sure symbol has the correct casing.

        if symbol.title() == testSymbol:
            return i;

    raise Exception("Error: Atomic symbol '{0}' not found in the internal periodic table.".format(symbol));

def AtomicNumberToSymbol(atomicNumber):
    if atomicNumber < len(PeriodicTable):
        return PeriodicTable[atomicNumber];
    else:
        raise Exception("Error: Atomic number {0} not assigned a symbol in the internal periodic table.".format(atomicNumber));
