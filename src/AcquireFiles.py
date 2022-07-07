
from Options import Options
from pathlib import Path

def acquireFiles(initializer):
  files = []
  folderPath = Path(initializer.INPUT_PATH)

  if not(folderPath.exists()):
    return
    
  
  for current_file in folderPath.iterdir():
    extension = current_file.suffix

    #RDKit can use mol2, mol, smiles, smarts, FASTA, HELM, PDB, PNG
    formats = [".fasta",".yaml",".mol2",".mol",".pbd",".sma",".smi",".tpl"]
    if not(extension in formats):
      print(f"[Error] emolFrag 2.0 only accpets fasta, yaml, mol2, mol, pbd, sma, smi, tpl files")
    else:
      files.append(Path(initializer.INPUT_PATH + "/" + current_file.name))
          
  return files
  
  #if __name__ == "__main__":
  