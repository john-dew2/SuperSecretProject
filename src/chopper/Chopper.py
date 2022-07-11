from rdkit import Chem
import utilities
import constants
import DoubleBondRecoonector
import Fragmenter
import FragmentPartitioner    

#
# Main chop routine:
#     (1) Fragment with BRICS
#     (2) Reconnect C.2 = C.2 double bonds broken by BRICS
#     (3) Identify bricks and linkers
#     (4) Reconnect adjacent linkers that share a dummy bond
#
# @input: An rdkit molecule with all hydrogens removed
# @output: list of (rdkit) brick and linker objects
#
def chop(stripped_rdkit_mol):

    #
    # (1) Fragment
    #
    try: 
        rdkit_fragments = fragmenter.fragment(stripped_rdkit_mol)

    except RDKitError:
        return []

    #
    # (2) Reconnect C.2 = C.2 Double Bonds
    #
    fragCountBefore = len(rdkit_fragments)
    rdkit_fragments = DoubleBondReconnector.reconnect(molecule.GetRDKitObject(),
                                                      rdkit_fragments)

    if len(rdkit_fragments) != fragCountBefore:
        print(f'Reconnected {fragCountBefore - len(rdkit_fragments)} C.2-C.2 double bonds.')

    #
    # (3) Distribute the fragments into bricks / linkers
    #
    bricks, linkers = FragmentPartitioner.partition(fragments)

    #
    # (4) 'Adjacent' linkers will be combined into single linkers
    #
    linkers = LinkerCombiner.combineLinkers(molecule.GetRDKitObject(), bricks, linkers)
    
    return bricks, linkers

#
# Create and add local Brick/Linker objects to the database; report results
#
def processFragments(level, bricks, linkers, parent):

    results = database.addAll([Brick(fragment, parent) for fragment in bricks])
    
    emit(level, f'Added {result.count(True)} TC-unique bricks; \
         {result.count(False)} were TC-redundant')

    results = database.addAll([Linker(fragment, parent) for fragment in linkers])
    
    emit(level, f'Added {result.count(True)} TC-unique linkers; \
         {result.count(False)} were TC-redundant')

#
# Chop many molecules
#
# @input: list of Molecule objects (containing rdkit objects)
# @output: MoleculeDatabase containing all fragments
#
def chopall(level, molList):

    database = MoleculeDatabase()

    for parent_mol in molList:

        emit(level, f'Processing molecule{parent_mol.GetFileName()}')

        # Remove all hydrogens from our molecule for simplicity
        stripped_rdkit = parent_mol.GetRDKitObject().RemoveAllHs()

        # Actually chop and store the molecular fragments in the database
        bricks, linkers = chop(stripped_rdkit))
        processFragments(bricks, linkers, parent_mol)      

    return database