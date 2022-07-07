
from Options import Options

def acquireFiles(initializer):
    files = []
    filenames = []
    folderPath = Path(initializer.INPUT_PATH)

    if folderPath.exists():
        return
    
    #initialize two lists; one with the full path of a file and one with the file names
    for file in folderPath.iterdir():
        #file extension will help filter bad data
        extension = file.suffix

        #RDKit can use mol2, mol, smiles, smarts, FASTA, HELM, PDB, PNG
        formats = [".fasta",".yaml",".mol2",".mol",".pbd",".sma",".smi",".tpl"]
        if ((extension in formats) == False):
            print("[Error] emolFrag 2.0 only accpets .fasta,.yaml,.mol2,.mol,.pbd,.sma,.smi,.tpl")
        else:
            filenames.append(file.name)
            files.append(Path(initializer.INPUT_PATH + "/" + file.name))
            
    return files