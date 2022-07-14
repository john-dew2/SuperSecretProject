from pathlib import Path

import sys
from eMolFrag2.unittests import utilities
from eMolFrag2.src.representation import Molecule, MoleculeDatabase


def run_TCEquiv(mol1path, mol2path, expec_result):

    mol1 = utilities.getRDKitMolecule(mol1path,Path(mol1path).suffix)
    mol2 = utilities.getRDKitMolecule(mol2path,Path(mol2path).suffix)
    
    md = MoleculeDatabase.MoleculeDatabase(given_tc = 1.0)

    assert md._TCEquiv(mol1, mol2) == expec_result

def run_TCEquivTests():

    cwd = Path.cwd().joinpath("eMolFrag2", "unittests", "data")

    # Test tc unique molecules (antibiotics with tc < 1.0 )
    mol1path = cwd.joinpath("uniqueMol(SMI)/DB00415.smi")
    mol2path = cwd.joinpath("uniqueMol(SMI)/DB01208.smi")
    mol3path = cwd.joinpath("uniqueMol(SMI)/DB04626.smi")

    run_TCEquiv(mol1path, mol2path, False)
    run_TCEquiv(mol1path, mol3path, False)
    run_TCEquiv(mol2path, mol3path, False)


    # Test similar pair antibiotics 
    # set 1
    mol1 = cwd.joinpath("similarPairMols/1/DB00452.mol2")
    mol2 = cwd.joinpath("similarPairMols/1/DB01421.mol2") 
    run_TCEquiv(mol1, mol2, True)    # kekulize / sanitise problem 

    # set 2 
    # kekulize / sanitise problem
    #mol3 = cwd.joinpath("similarPairMols/2/DB01137.mol2")
    #mol4 = cwd.joinpath("similarPairMols/2/DB01165.mol2")
    mol3 = cwd.joinpath("similarPairSMI/2/DB01137.smi")
    mol4 = cwd.joinpath("similarPairSMI/2/DB01165.smi")
    run_TCEquiv(mol3, mol4, True)
    
    # prove sets of pairs of molecules do not have tc_equiv
    run_TCEquiv(mol1, mol3, False) 
    run_TCEquiv(mol1, mol4, False) 
    run_TCEquiv(mol2, mol3, False) 
    run_TCEquiv(mol2, mol4, False) 

    # set 3 
    # kekulize / sanitise problem
    #mol5 = cwd.joinpath("similarPairMols/3/DB12447.mol2")
    #mol6 = cwd.joinpath("similarPairMols/3/DB16219.mol2")
    mol5 = cwd.joinpath("similarPairSMI/3/DB12447.smi")
    mol6 = cwd.joinpath("similarPairSMI/3/DB16219.smi")
    #run_TCEquiv(mol5, mol6, True)     


    # Test tc unique molecules (tc < 1.0)
    uniqueMol1 = cwd.joinpath("uniqueMol(SMI)/DB00415.smi")
    uniqueMol2 = cwd.joinpath("uniqueMol(SMI)/DB01208.smi")
    uniqueMol3 = cwd.joinpath("uniqueMol(SMI)/DB04626.smi")
    uniqueMol4 = cwd.joinpath("uniqueMol(SMI)/DB11774.smi")
    uniqueMol5 = cwd.joinpath("uniqueMol(SMI)/DB13499.smi")

    run_TCEquiv(uniqueMol1, uniqueMol2, False)
    run_TCEquiv(uniqueMol1, uniqueMol3, False)
    run_TCEquiv(uniqueMol1, uniqueMol4, False)
    run_TCEquiv(uniqueMol1, uniqueMol5, False)

    run_TCEquiv(uniqueMol2, uniqueMol3, False)
    run_TCEquiv(uniqueMol2, uniqueMol4, False) 
    run_TCEquiv(uniqueMol2, uniqueMol5, False)

    run_TCEquiv(uniqueMol3, uniqueMol4, False)
    run_TCEquiv(uniqueMol3, uniqueMol5, False)
    
    run_TCEquiv(uniqueMol4, uniqueMol5, False)



def run_add(md, molPath, expec_result):

    # Create rdkit object from file path           extension
    mol = utilities.getRDKitMolecule(molPath, Path(molPath).suffix)
    assert md.add(mol) == expec_result 

def run_addTests():

    cwd = Path.cwd().joinpath("eMolFrag2", "unittests", "data")
    md = MoleculeDatabase.MoleculeDatabase(given_tc = 1.0)

    # Molecules with tc = 1.0
    mol1 = cwd.joinpath("similarPairMols/1/DB00452.mol2")
    mol2 = cwd.joinpath("similarPairMols/1/DB01421.mol2")

    run_add(md, mol1, True)
    run_add(md, mol1, False)
    run_add(md, mol2, False)

    # Molecules with tc < 1.0
    uniqueMol1 = cwd.joinpath("uniqueMol(SMI)/DB00415.smi")
    uniqueMol2 = cwd.joinpath("uniqueMol(SMI)/DB01208.smi")
    uniqueMol3 = cwd.joinpath("uniqueMol(SMI)/DB04626.smi")
    uniqueMol4 = cwd.joinpath("uniqueMol(SMI)/DB11774.smi")
    uniqueMol5 = cwd.joinpath("uniqueMol(SMI)/DB13499.smi")

    run_add(md, uniqueMol1, True)
    run_add(md, uniqueMol2, True)
    run_add(md, uniqueMol3, True)
    run_add(md, uniqueMol4, True)
    run_add(md, uniqueMol5, True)
    

# return the list of UNIQUE database
def run_addAll(md, moleculesPath, length):
    molecules = []
    for molPath in moleculesPath: 
      mol = utilities.getRDKitMolecule(molPath, Path(molPath).suffix)
      molecules.append(mol)

    # frequency = md.addAll(molecules)
    # print("Size of list:\n", len(frequency))
    # print(frequency)

    # Test the size of the list 
    assert len(md.addAll(molecules)) == length

def run_addAllTests():
    cwd = Path.cwd().joinpath("eMolFrag2", "unittests", "data")
    md = MoleculeDatabase.MoleculeDatabase(given_tc = 1.0)
    moleculesPath = []

    # Test: adding 5 unique molecules 
    uniqueMol1 = cwd.joinpath("uniqueMol(SMI)/DB00415.smi")
    uniqueMol2 = cwd.joinpath("uniqueMol(SMI)/DB01208.smi")
    uniqueMol3 = cwd.joinpath("uniqueMol(SMI)/DB04626.smi")
    uniqueMol4 = cwd.joinpath("uniqueMol(SMI)/DB11774.smi")
    uniqueMol5 = cwd.joinpath("uniqueMol(SMI)/DB13499.smi")
    
    moleculesPath.extend([uniqueMol1, uniqueMol2, uniqueMol3, uniqueMol4, uniqueMol5])
    run_addAll(md, moleculesPath, 5)

    # Test: add pairs of similar molecules (tc = 1.0)
    mol1 = cwd.joinpath("similarPairSMI/1/DB00452.smi")
    mol2 = cwd.joinpath("similarPairSMI/1/DB01421.smi")
    moleculesPath.clear()
    moleculesPath.extend([mol1, mol2])
    run_addAll(md, moleculesPath, 1)

    # Test: add 2nd pair of similar molecules
    mol3 = cwd.joinpath("similarPairSMI/2/DB01137.smi")
    mol4 = cwd.joinpath("similarPairSMI/2/DB01165.smi")
    moleculesPath.append(mol3)
    moleculesPath.append(mol4)
    print("# of path: ", len(moleculesPath))
    run_addAll(md, moleculesPath, 2)


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
    tests = { "_TCEquiv" : run_TCEquivTests,
              "add" : run_addTests,
              "addAll" : run_addAllTests
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
