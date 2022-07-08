from rdkit import Chem
import fragmenter
import constants
        
def reconnectC2C2DoubleBonds(rdkit_mol_list):
    raise NotImplementedError

#
# Return the number of atoms that are
#    (1) dummy atoms or
#    (2) radicals or
#    (3) hydrogens
#
# @input: rdkit fragment object
# @output: True if it meets Brick criteria (False indicates a Linker)
#
def CountAtomsToBeRemoved(rdkit_fragment):

    count = 0
    for atom in m.GetAtoms():

        sym = atom.GetSymbol()

        if sym == constants.DUMMY_ATOM_STR or \
           sym == constants.RADICAL_ATOM_STR or \
           sym == constants.HYDROGEN_ATOM_STR:
            count = count + 1

    return count

#
# Brick criteria:
#    Number of 'meaningful' atoms is 4 or more
#
# @input: rdkit fragment object
# @output: True if it meets Brick criteria (False indicates a Linker)
#
def isBrick(rdkit_frag):
    return rdkit_frag.GetNumAtoms() - CountAtomsToBeRemoved(rdkit_frag) >= constants.BRICK_MINIMUM_NUM_ATOMS


#
# @input: list of rdkit molecule fragment objects
# @output: a tuple of two lists (bricks, linkers) both lists are rdkit objects
#
def identifyBricksAndLinkers(rdkit_fragments):

    bricks = []
    linkers = []
    for rdkit_frag in rdkit_fragments:

        if isBrick(rdkit_frag):
            bricks.append(rdkit_frag)

        else:
            linkers.append(rdkit_frag)
            
    return bricks, linkers

#
# Main chop routine:
#     (1) Fragment with BRICS
#     (2) Reconnect C.2 = C.2 double bonds broken by BRICS
#     (3) 
#
# @input: A Molecule object (which contains an rdkit molecule object)
# @output: (Brick object list, Linker object list)
#
def chop(molecule):

    #
    # (1) Fragment
    #
    try: 
        rdkit_fragments = fragmenter.fragment(molecule.GetRDKitObject())

    except RDKitError:
        return []

    #
    # (2) Reconnect C.2 = C.2 Double Bonds
    #
    rdkit_fragments = reconnectC2C2DoubleBonds(rdkit_fragments)

    if len(mfl) != len(mfl2):
        print("Reconnected a double bond between", len(mfl) - len(mfl2) + 1, "fragments")

    #
    # (3) Distribute the fragments into bricks / linkers
    #
    bricks, linkers = identifyBricksAndLinkers(fragments)

    #
    # (4) Create 'local' molecule representations
    #
    return localizeMolecules(bricks, linkers, molecule)
       
#
# Chop many molecules
#
# @input: list of Molecule objects
#

def chopall(molList):

    for mol in molList:

        print(f"Processing molecule{molecule.GetFileName()}")

        bricks, linkers = chop(mol)
        
        
        
        