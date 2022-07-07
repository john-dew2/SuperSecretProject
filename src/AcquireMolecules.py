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
    #kekulize = false
    
    if (extension == ".fasta"):
        return Chem.MolFromFASTA(contents) 
    if (extension == ".yaml"):
        return Chem.MolFromHELM(contents) 
    if (extension == ".mol2"):
        return Chem.MolFromMol2Block(contents)
    if (extension == ".mol"):
        return Chem.MolFromMolBlock(contents)
    if (extension == ".pbd"):
        return Chem.MolFromPDBBlock(contents) 
    if (extension == ".sma"):
        return Chem.MolFromSmarts(contents) 
    if (extension == ".smi"):
        return Chem.MolFromSmiles(contents)
    if (extension == ".tpl"):
        return Chem.MolFromTPLBlock(contents)
        
    return None
    
def acquireMolecules(files):
    #parse through each file to convert to RDKit molecule        
    for current_file in files:
        #get the contents of the file and the file type (extension) for processing
        file_contents = fileToSting(current_file)
        extension = current_file.suffix
            
        #get the molecule
        try:
            molecule = convertToRDkit(file_contents, extension)
        except:
            #print(f"RDKit failed to read {}" current_file.name)
                
        #if the molecule didnt process let the user know
        if (molecule == None):
            print("[Error] RDKit was unable to convert", current_file.name,  "to a RDKit object")
        #otherwise add it to our dataset and update the filenames we have
        else:
            molObject = Molecule(molecule, current_file.name)
            data.append(molObject)
            
    return data