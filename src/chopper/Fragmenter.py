from rdkit.Chem import BRICS
from rdkit import Chem

#
# (1) Determine where BRICS is going to cleave (FindBRICSBonds)
# (2) Determine if any of the bonds are C.2=C.2 double bonds. If so, remove bond from list of bonds to cleave
# (3) Do we have any fragments that will be linkers that are linker-linker adjacent? Remove those cleave bonds
#


#
# Does the bond consist of a carbon (C) and dummy created by BRICS (*)?
#
# @input: bond indicated by indices (startAtomIndex)
# @input: Original RDKit molecule
#
def isCCBond(startAtom, endAtom):
    return startAtom.GetSymbol() == constants.CARBON_ATOM_STR and \
           endAtom.GetSymbol() == constants.CARBON_ATOM_STR

#
# Does this rdkit molecule (fragment) have a double bond
#
def isDoubleBond(bond):
   return bond.GetbondType() == Chem.rdchem.BondType.DOUBLE

#
# Does this rdkit molecule (fragment) have a C.2 = C.2 double bond
#
# A bond of the form ((start, end), ('str', 'str'))
#
def isC2C2DoubleBond(level, bond, rdkit_mol):

   #
   # Get the bond object from this molecule
   #
   startAtomIndex = bond[0][0]
   endAtomIndex = bond[0][1]
   
   if not isinstance(startAtomIndex, int) or not isinstance(startAtomIndex, int):
      utilities.emitError(level, f'Internal error: Expected (int, int) bond format not found: {bond}')

   actual_mol_bond = rdkit_mol.GetBondBetweenAtoms(startAtomIndex, endAtomIndex)
   
   if actualmolBond is None:
       utilities.emitError(level, f'Internal error: Bond {(startAtomIndex, endAtomIndex)} not found as expected.')

   #
   # Do we have a double bond? And, is C.2=C.2
   return isDoubleBond(actual_mol_bond) and \
          isCCBond(rdkit_mol.GetAtomWithIdx(startAtomIndex), rdkit_mol.GetAtomWithIdx(endAtomIndex))

#
# @input: list of bonds that would be cleaved
# @input: rdkit molecule object
# @output: (valid cleave bonds, removed cleave bonds)
#
def AnalyzeCleaveC2C2DoubleBonds(level, bonds, rdkit_mol):

    cleaveBonds = []
    removedBonds = []

    for bond in bonds:
        (removedBonds if isC2C2DoubleBond(level, bond, rdkit_mol) else cleaveBonds).append(bond)

    return cleaveBonds, removedBonds

#
# Is the proposed fragment a linker (by our definition)?
# @input: RDKIT molecule
# @output True if fragment contains fewer than (4) atoms excluding dummy atoms
#
def isLinker(level, rdkit_mol):
    # Get the atom types 
    atomSyms = [atom.GetSymbol() for atom in rdkit_mol.GetAtoms()]

    dummyCount = atomsSyms.count(DUMMY_ATOM_STR)
    rCount = atomsSyms.count(RADICAL_ATOM_STR)

    if rCount != 0:
        utilities.emitError(level, "R defines dummy atoms (unexpectedly) {rCount} vs. {dummyCount}")

    # A linker is defined as having fewer than X atoms (that are non-dummy and non-radical)
    return len(rdkit_mol.GetAtoms()) - (dummyCount + rCount) <= constants.LINKER_MAXIMUM_NUM_ATOMS

#
# Do we have any cleave bonds that would connect two linkers together?
#
def AnalyzeCleaveLinkerLinkerBonds(level, bonds, mol):

    # If we cleaved right now with N bonds, we should acquire an ORDERED list of N+1 fragments
    #
    # We should have fragments (F) and bonds (B) with their indices
    #
    #    |-----------------|---------|-----|-------------|
    #            F0        B0   F1   B1 F2 B2    F3
    #
    try:
        broken_mol = BRICS.BreakBRICSBonds(m, bonds)
        fragments = Chem.BRICS.GetMolFrags(broken_mol, asMols = True, sanitizeFrags = False)
    except:
        utilities.emitError(f'Linker-Linker analysis failed to fragment with BRICS; returning all cleave bonds')
        return cleaveBonds, []

    if len(bonds) + 1 != len(fragments):
        emitError(level, f'Unepxected number of fragments when cleaved; expected {len(bonds) + 1}, found {len(fragments)}')

    #
    # Identify subsequences of linkers; remove cleave bonds if a subsequence is found
    #
    cleaveBonds = []
    removedBonds = []

    for index_f in range(len(fragments) - 1):
        (removedBonds if isLinker(level, fragments[i]) and isLinker(level, fragments[i+1]) else cleaveBonds).append(bonds[i])

    return cleaveBonds, removeBonds

#
# @input:  rdkit molecule that has been stripped of hydrogens
# @output: List of rdkit objects (fragments)
#
def identifyCleaveBonds(level, mol):

    #
    # Acquire the bonds that will be cleaved.
    #
    cleaveBonds = list[Chem.BRICS.FindBRICSBonds(mol)] # FindBRICSBonds is a generator
    
    #
    # Check cleave bonds for C.2=C.2 double bond
    #
    validCleaveBonds, removed = AnalyzeCleaveC2C2DoubleBonds(level, cleaveBonds, mol)

    res_len = (len(validCleaveBonds) + len(removed))
    if len(cleaveBonds) != len(validCleaveBonds) + len(removed):
        utilities.emitError(level, f'Internal C.2=C.2 cleave analysis error; expected {len(cleaveBonds)}, but found {res_len}.')

    if len(cleaveBonds) != len(validCleaveBonds):
        utilities.emit(level, f'Removed C.2=C.2 double bond(s) {" ".join(removed)} from being cleaved')

    #
    # Check cleave bonds for linker-linker connections
    #
    finalCleaveBonds, removed = AnalyzeCleaveLinkerLinkerBonds(level, validCleaveBonds, mol)

    res_len = (len(finalCleaveBonds) + len(removed))
    if len(validCleaveBonds) != len(finalCleaveBonds) + len(removed):
        utilities.emitError(level, f'Internal Linker-Linker cleave analysis error; expected {len(validCleaveBonds)}, but found {res_len}.')

    if len(finalCleaveBonds) != len(validCleaveBonds):
        utilities.emit(level, f'Removed Linker-Linker bond(s) {" ".join(removed)} from being cleaved')

    return finalCleaveBonds


#
# @input:  rdkit molecule that has been stripped of hydrogens
# @output: List of rdkit objects (fragments)
#
def fragment(level, mol):

    cleaveBonds = identifyCleaveBonds(level, mol.Get)

    try:
        utilities.emit(f'Cleaving with {cleaveBonds}')

        broken_mol = BRICS.BreakBRICSBonds(m, cleaveBonds)
    
        return Chem.BRICS.GetMolFrags(broken_mol, asMols = True, sanitizeFrags = False)

    except:
    
        utilities.emitError(f'Molecule failed to fragment with BRICS')

    return fragments