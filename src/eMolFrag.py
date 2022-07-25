import sys
from pathlib import Path

from SuperSecretProject.src.input import MoleculeFileReader, MoleculeReader, Options
from SuperSecretProject.src.chopper import Chopper
from SuperSecretProject.src.representation import Molecule
from SuperSecretProject.src.output import writer
from SuperSecretProject.src.utilities import logging

def main():
    """
        eMolFrag

        (1) Parse input arguments
        (2) Read input files into Molecule objects
        (3) Fragment molecules
        (4) Output fragments as specified by the user

    """
    options = Options.Options()
    if not options.isRunnable():
        logging.logger.error(f'Command-line arguments failed to parse; execution of eMolFrag will stop.')
        return

    dataset = []
    
  	# Verify Tools and Parse Command Line

    # Get files
    mol_files = MoleculeFileReader.getFiles(options)
    logging.logger.info(f'{len(mol_files)} files to be processed.')
    
    # Get molecules
    molecules = MoleculeReader.getMolecules(mol_files)
    logging.logger.info(f'{len(molecules)} molecules to be chopped.')
 
    # CHOP
    brick_db, linker_db = Chopper.chopall(molecules)
    
    # Output fragments
    logging.logger.info(f'{brick_db.numUnique()} unique bricks among {brick_db.numAllMolecules()} bricks')
    logging.logger.info(f'{linker_db.numUnique()} unique linkers among {linker_db.numAllMolecules()} linkers')

    writer.write(options, brick_db, linker_db)


if __name__ == '__main__':
    main()