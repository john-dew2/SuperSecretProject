

from pathlib import Path

import utilities
import MoleculeDatabase

def run_TCEquiv(pathForMol1, pathForMol2):

    mol1 = utilities.getRDKitMolecule(Path(PathForMol1))
    mol2 = utilities.getRDKitMolecule(Path(PathForMol2))

    assert MoleculeDatabase._TCEquiv(mol1, mol2) == expec_result


def run_TCEquivTest():
    # Test 31738.mol2-001
    mol1Path = "C:\Users\tchen\Documents\GitHub\eMolFrag-2\unittests\data\linkersWithTC1.0\l-CHEMBL31738.mol2-001\l-CHEMBL31738.mol2-000.sdf"
    mol2Path = "C:\Users\tchen\Documents\GitHub\eMolFrag-2\unittests\data\linkersWithTC1.0\l-CHEMBL31738.mol2-001\l-CHEMBL427852.mol2-000.sdf"
    mol3Path = "C:\Users\tchen\Documents\GitHub\eMolFrag-2\unittests\data\linkersWithTC1.0\l-CHEMBL31738.mol2-001\l-CHEMBL427852.mol2-001.sdf"

    run_TCEquiv(mol1Path, mol2Path, True)
    run_TCEquiv(mol1Path, mol3Path, True)
    run_TCEquiv(mol2Path, mol3Path, True)

def runtests():

    successful = []
    try:
       successful.append("testName")
    except :
        print("test fail")
        
    print(f'{__name__} unit tests are successful.')

if __name__ == "__main__":
    runtests()
