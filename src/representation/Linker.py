from Molecule import Molecule

class Linker(Molecule):

    def __init__(self, rdkit_object, file_name):
<<<<<<< Updated upstream:src/representation/Linker.py
    
        Molecule.__init__(self, rdkit_object, file_name):
=======
        Molecule.__init__(self, rdkit_object, file_name):

    def __init__(self, rdkit_object, rdkit_parent):
        Molecule.__init__(self, rdkit_object, rdkit_parent):
>>>>>>> Stashed changes:src/Linker.py
