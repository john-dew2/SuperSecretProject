import sys

sys.path.insert(0, 'C:\Users\tchen\Documents\GitHub\eMolFrag-2\src\representation\MoleculeDatabase.py')
from pathlib import Path


import MoleculeDatabase
import utilities


def runtests():

    successful = []
    try:
       successful.append("testName")
    except :
        print("test fail")
        
    print(f'{__name__} unit tests are successful.')

if __name__ == "__main__":
    runtests()
