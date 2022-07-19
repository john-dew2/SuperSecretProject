from rdkit import Chem
from eMolFrag2.src.utilities import logging
import Preprocessor
import Deconstructor
import Connectivity

def chop(rdkit_mol):
    """
        0. We work on a copy of the input molecule with hydrogens removed and
           atom type information within the molecule 

        Chopping consists of the following algorithm:
            1. Build graph of molecule
            2. Find where BRICS would cleave (BRICS bonds) [snips in the code]
                 * We modify FindBRICSBonds code to ensure that L7 bonds are never
                   cleaved (according eMolFrag v1.0)
                 * Modified in BRICS_custom.py
            3. Break graph into fragments
                 Take the graph of molecule, substract proposed BRICS bonds
                 The resulting molecule is a disconnected set of subgraphs.
                 Each of these subgraphs is a fragment.
                 Those fragments with fewer than 4 atoms are the first phase of
                 linkers; all others are bricks
            4. Combine sequences of Linker-Linker fragments into single a Linker

            5. Compute connectivity of free radicals (from each snip calc atomtype of tuple pairs)

            6. Break each fragment into Rdkit.Mol objectPerform actual break of chop
            
        Input: An rdkit molecule (ideally, with AtomType info from mol2 format)
            
        Output: list of Rdkit.Mol Bricks, list of Rdkit.Mol Linkers
    """
    # (0)
    mol = preprocess(rdkit_mol)

    #
    # Steps (1) - (4)
    # Deconstruct our molecule into the known fragments
    # These are sets of atom indices (bricks, linkers); snips are a set of bonds.
    #    * All such information 'references' the molecule, but does not modify it
    #
    bricks, linkers, snips = Deconstructor.deconstruct(mol)

    #
    # (5): Compute connectivity among free radicals
    #      This computation modifies the molecule with property information
    Connectivity.compute(mol, snips)

    #
    # (6): Perform the actual chop
    #
    return Fragmenter.fragmentAll(mol, bricks, linkers)

def chopall(mols):
    """
        Chop many molecules
        
        @input: list of Molecule objects (containing rdkit objects)
        
        @output: MoleculeDatabase of brick fragments
        @output: MoleculeDatabase of linker fragments        
    """
    brick_db = MoleculeDatabase(constants.DEFAULT_TC_UNIQUENESS)
    linker_db = MoleculeDatabase(constants.DEFAULT_TC_LINKER_UNIQUENESS)

    for mol in mols:

        logging.logger.debug(f'Processing molecule{mol.GetFileName()}')

        #
        # Chop
        #
        bricks, linkers = chop(mol.getRDKitObject())

        #
        # Process the results
        #
        results = database.addAll([Brick(b, mol) for b in bricks])
    
        logging.logger.debug(f'Added {result.count(True)} TC-unique bricks; \
                             {result.count(False)} were TC-redundant')

        results = database.addAll([Linker(ell, mol) for ell in linkers])

        logging.logger.debug(f'Added {result.count(True)} TC-unique linkers; \
                             {result.count(False)} were TC-redundant')

    return brick_db, linker_db