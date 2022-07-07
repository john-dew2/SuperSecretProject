
import sys
from Molecule import Molecule
from Options import Options
from pathlib import Path
from rdkit import Chem
from AcquireFiles import acquireFiles
from AcquireMolecules import acquireMolecules

def main():
    dataset = []
    
	#Verify Tools and Parse Command Line
    args = sys.argv
    initializer = Options()
    initializer.parseCommandLine(args)
    
    #Input System
    files = acquireFiles(initializer)
    dataset = acquireMolecules(files)
    
    #Process
    
    #Post-Process
	
	





if __name__ == "__main__":
	main()