import shutil
from pathlib import Path

from eMolFrag2.src.representation import MoleculeDatabase
from eMolFrag2.src.utilities import logging
from eMolFrag2.src.utilities import constants

def prepareDirectory(out_path):
    """
        If the directory does not exist, create it.
        If the directory does exist, clean it.        
    """
    if out_path.is_file():
        logging.logger.error(f'Output path {str(out_path)} is a file, not a directory.')
        raise ValueError(f'Malformed output specification path {str(out_path)}')

    # Delete the directory, if needed
    if not out_path.exists():
        logging.logger.info(f'Output path {str(out_path)} does not exist; will be created.')        
    else:
        logging.logger.warning(f'Output path {str(out_path)} exists; all contents will be deleted.')
        shutil.rmtree(out_path)
    
    # Recreate the diretory
    out_path.mkdir()

def writeSingleFile(indicator, name, out_dir, mols, extension = constants.SDF_FORMAT_EXT):
    """
        indicator --  'u' --> unique or 'a' --> all
        name -- main part of the output file: 'bricks' or 'linkers'
        out_dir -- output directory path
        mols -- the actual Molecule objects to write
    """
    file_name = f'{indicator}{name}{extension}'

    logging.logger.debug(f'Writing file {file_name}')

    # Delimiter needed? Or does SDWriter put it there?
    text = constants.SDF_MOLECULE_DELIMITER.join([mol.toSDF() for mol in mols])

    out_path = out_dir / file_name

    out_path.touch()
    out_path.open('w').write(text)
    out_path.close()

def writeIndividualFiles(out_dir, mols, extension = constants.SDF_FORMAT_EXT):
    """
        indicator --  'u' --> unique or 'a' --> all
        name -- main part of the output file: 'bricks' or 'linkers'
        out_dir -- output directory path
        mols -- the actual Molecule objects to write
    """
    for mol in mols:    
        out_path = out_dir / mol.getFileName()
        out_path.touch() # Needed?
        out_path.open('w').write(mol.toSDF())
        out_path.close()

def write(options, brick_db, linker_db):
    """
        Main output routine
        The focus is what fragments (unique OR all) and format how to
        output it (many files OR a single file).
    """
    out_dir = Path(options.OUTPUT_PATH)
    prepareDirectory(out_dir)

    logging.logger.debug(f'Writing to directory {str(out_dir)}')

    bricks_to_write = []
    linkers_to_write = []
    indicator = ''

    # Only unique fragments wanted
    if options.UNIQUE:
        indicator = constants.FILE_OUTPUT_UNIQUE_INDICATOR
        bricks_to_write = brick_db.GetUniqueMolecules()
        linkers_to_write = linkers_db.GetUniqueMolecules()

    # All fragments wanted
    else:
        indicator = constants.FILE_OUTPUT_ALL_INDICATOR
        bricks_to_write = brick_db.GetAllMolecules()
        linkers_to_write = linker_db.GetAllMolecules()

    # Write all fragments to their own files
    if options.INDIVIDUAL:
        writeIndividualFiles(out_dir, bricks_to_write + linkers_to_write)

    # Write all fragments to a single brick and a signle linker file
    else:
        writeSingleFile(indicator, BRICK_SINGLE_FILE_OUTPUT_NAME, out_dir, bricks_to_write)
        writeSingleFile(indicator, LINKER_SINGLE_FILE_OUTPUT_NAME, out_dir, linkers_to_write)