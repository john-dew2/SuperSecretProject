
import sys
import os
from Options import Options
from pathlib import Path
from rdkit import Chem

#takes the contents of a file and puts it in a string for processing
def fileToSting(file):
    retString = ""
    #get all the lines
    with open(file) as f:
      lines = f.readlines()
    f.close()
    #concatenate each line to one string
    for line in lines:
        retString += line
        
    return retString
   
#converts the contetns of a file to their respective molecule
def convertToRDkit(contents, extension):
    molecule = None
    
    if (extension == ".fasta"):
        molecule = Chem.MolFromFASTA(contents) 
    if (extension == ".yaml"):
        molecule = Chem.MolFromHELM(contents) 
    if (extension == ".mol2"):
        molecule = Chem.MolFromMol2Block(contents)
    if (extension == ".mol"):
        molecule = Chem.MolFromMolBlock(contents)
    if (extension == ".pbd"):
        molecule = Chem.MolFromPDBBlock(contents) 
    if (extension == ".sma"):
        molecule = Chem.MolFromSmarts(contents) 
    if (extension == ".smi"):
        molecule = Chem.MolFromSmiles(contents)
    if (extension == ".tpl"):
        molecule = Chem.MolFromTPLBlock(contents)
        
    return molecule

#Takes the input files and turns them into molecular objects
def handleInput(initializer):
    files = []
    original_filenames = []
    filenames = []
    data = []

    if os.path.exists(initializer.INPUT_PATH):
        folder =  initializer.INPUT_PATH
        
        #initialize two lists; one with the full path of a file and one with the file names
        for filename in os.listdir(folder):
            #file extension will help filter bad data
            extension = filename[filename.find("."):]

            #RDKit can use mol2, mol, smiles, smarts, FASTA, HELM, PDB, PNG
            if ((extension != ".mol2")and(extension != ".smi")):
                print("[Error] emolFrag 2.0 only accpets .mol2 or .smi formats")
            else:
                original_filenames.append(filename)
                files.append(folder + "/" + filename)
                
        #parse through each file to convert to RDKit molecule        
        for current_file in files:
            #get the contents of the file and the file type (extension) for processing
            file_contents = fileToSting(current_file)
            extension = current_file[current_file.find("."):]
            
            #get the molecule
            molecule = convertToRDkit(file_contents, extension)
            
            #if the molecule didnt process let the user know
            if (molecule == None):
                print("[Error] RDKit was unable to convert", original_filenames[files.index(current_file)],  "to a RDKit object")
            #otherwise add it to our dataset and update the filenames we have
            else:
                data.append(molecule)
                filenames.append(original_filenames[files.index(current_file)])
                
    return data

def main():
    dataset = []
	#Verify Tools and Parse Command Line
    args = sys.argv
    initializer = Options()
    
    initializer.parseCommandLine(args)
    
    #Input System
    dataset = handleInput(initializer)
    
    #Process
    
    #Post-Process
	
	





if __name__ == "__main__":
	main()