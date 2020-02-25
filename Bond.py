class Bond:

    def __init__(self):
        self.peerA = None
        self.peerB = None
        self.bond_count = 0
        self.length_range_low = 0.0
        self.length_range_high = 0.0

    def load(self, length_range_low, length_range_high):
        self.length_range_low = length_range_low
        self.length_range_high = length_range_high
        return self

    def load(self, bond_count, peerA, peerB, length_range_low, length_range_high):
        self.peerA = peerA
        self.peerB = peerB
        self.bond_count = bond_count
        self.length_range_low = length_range_low
        self.length_range_high = length_range_high
        return self

    @staticmethod
    def does_length_matches_to_bonds(bond_length, bond_definition_list):
        flag = False
        for bond in bond_definition_list:
            if bond.length_range_low <= bond_length <= bond.length_range_high:
                flag = True
                break
            else:
                flag = False
        return flag



