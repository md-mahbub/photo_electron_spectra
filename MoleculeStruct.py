#test change
class MoleculeStruct:

    def __init__(self):
        self.atom_definition_map = None
        self.bond_definition_list = None

    def load(self, atom_definition_map, bond_definition_list):
        self.atom_definition_map = atom_definition_map
        self.bond_definition_list = bond_definition_list
        return self

    def get_atom_definition_map(self):
        return self.atom_definition_map

    def get_bond_definition_list(self):
        return self.bond_definition_list
