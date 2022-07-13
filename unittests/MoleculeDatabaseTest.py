from pathlib import Path

import sys

from eMolFrag2.unittests import utilities
from eMolFrag2.src.representation import MoleculeDatabase


def run_TCEquiv(pathForMol1, pathForMol2, expec_result):

    mol1 = utilities.getRDKitMolecule(pathForMol1,pathForMol1.suffix)
    mol2 = utilities.getRDKitMolecule(pathForMol2,pathForMol2.suffix)
    md = MoleculeDatabase.MoleculeDatabase(given_tc = 0.5)

    assert md._TCEquiv(mol1, mol2) == expec_result


def run_TCEquivTests():

    mol1path = Path("/content/eMolFrag2/unittests/data/uniqueMol(SMI)/DB00415.smi")
    mol2path = Path("/content/eMolFrag2/unittests/data/uniqueMol(SMI)/DB01208.smi")
    mol3path = Path("/content/eMolFrag2/unittests/data/uniqueMol(SMI)/DB04626.smi")

    # Test tc unique molecules (antibiotics with tc < 1.0 )
    run_TCEquiv(mol1path, mol2path, False)
    run_TCEquiv(mol1path, mol3path, False)
    run_TCEquiv(mol2path, mol3path, False)

    # Test similar pair antibiotics 
    mol1 = Path("/content/eMolFrag2/unittests/data/similarPairMols/DB00452.mol2")
    mol2 = Path("/content/eMolFrag2/unittests/data/similarPairMols/DB01421.mol2")

    run_TCEquiv(mol1, mol2, True)
    


def run(func):
    try: 
      func()
      return True 

    except: 
      return False

def runtest(test_name, test_func, successful, failed):
  (successful if run(test_func) else failed).append(test_name)

def runtests(printlevel):

    utilities.emit(printlevel, f'Executing {__file__} unit tests.')

    #
    # Define all tests as a Dictionary: {str-name, <function-to-execute>}
    # 
    tests = { "_TCEquiv" : run_TCEquivTests
              # "add" : run_addTests, 
              # "addAll" : run_addAllTests,
              # "GetAllMolecules" : run_GetAllMoleculesTests
            }
    #
    # Run 
    #
    successful = []
    failed = []

    for(test_name, test_func) in tests.items():
      runtest(test_name, test_func, successful, failed)

    #
    # Report 
    #
    if not failed: 
      utilities.emit(printlevel, f'{__name__} unit tests are successful.')

    else:
        for test in failed:
            utilities.emit(printlevel+1, f'Failed {test}.')

if __name__ == "__main__":
    runtests(0)
