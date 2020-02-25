import math
from Bond import *

class Atom:

    def __init__(self):
        self.atom_name = None
        self.atom_index = -1
        self.x_coordinate = None
        self.y_coordinate = None
        self.z_coordinate = None
        self.coefficient_count = 0
        self.coefficient_start_index = 0
        self.coefficient_end_index = 0
        self.neighbors = None
        self.possible_inhome_neighbors = None

    def load(self, atom_name, coefficient_count, possible_inhome_neighbors):
        self.atom_name = atom_name
        self.coefficient_count = coefficient_count
        self.possible_inhome_neighbors = possible_inhome_neighbors
        return self

    def load_atoms_from_geometry(atom_list_str):
        atoms = []
        atom_geometry = atom_list_str.split()
        for i, atom in enumerate(atom_geometry):
            if atom == 'C' or atom == 'H' or atom == 'O' or atom == 'N':
                atom = Atom()
                atom.atom_name = atom_geometry[i]
                atom.atom_index = i
                atom.x_coordinate = float(atom_geometry[i+1])
                atom.y_coordinate = float(atom_geometry[i+2])
                atom.z_coordinate = float(atom_geometry[i+3])
                atoms.append(atom)
                i = i + 4

        return atoms

    @staticmethod
    def find_euclidian_dist(atom1, atom2):
        dist = round(math.sqrt((atom1.x_coordinate - atom2.x_coordinate)**2
                               + (atom1.y_coordinate - atom2.y_coordinate)**2
                               + (atom1.z_coordinate - atom2.z_coordinate)**2), 3)
        return dist

    @staticmethod
    def assign_neighbors_to_atoms(bond_definition_list, atoms, molecule_structure):
        start_idx = 0
        for i, atom in enumerate(atoms):
            neighbors = []
            for j in range(len(atoms)):
                if i != j:
                    euclidian_dist = Atom.find_euclidian_dist(atoms[i], atoms[j])
                    matches = Bond.does_length_matches_to_bonds(euclidian_dist, bond_definition_list)

                    if matches:
                        neighbors.append(j)
            atom.neighbors = neighbors

            # find coefficient_start_index and coefficient_end_index of each atom
            # get the coefficient count based on their atom type,
            coefficient_count = (molecule_structure.get_atom_definition_map()[atom.atom_name]).coefficient_count
            # get the index position
            end_idx = start_idx + coefficient_count - 1

            atom.coefficient_start_index = start_idx
            atom.coefficient_end_index = end_idx

            start_idx = end_idx + 1


    # Deprecated
    @staticmethod
    def print_matrix(matrix):
        for i in range(len(matrix)):
            print()
            for j in range(len(matrix)):
                print(matrix[i][j], end =" ")


    #Deprecated
    @staticmethod
    def find_neighbor_matrix(dist_matrix, bond_definition_list):
        # 2D array - initialized by all zeros
        rows, cols = (len(dist_matrix), len(dist_matrix))
        neighbor_matrix = [[0 for i in range(cols)] for j in range(rows)]

        for i in range(len(dist_matrix)):
            for j in range(len(dist_matrix)):
                neighbor_matrix[i][j] = Bond.does_length_matches_to_bonds(dist_matrix[i][j], bond_definition_list)

        return neighbor_matrix

    # Deprecated
    @staticmethod
    def assign_neighbors_to_atoms1(dist_matrix, bond_definition_list, atoms):
        # 2D array - initialized by all zeros
        rows, cols = (len(dist_matrix), len(dist_matrix))
        neighbor_matrix = [[0 for i in range(cols)] for j in range(rows)]


        for i, atom in enumerate(atoms):
            neighbors = []
            for j in range(len(atoms)):
                matches = Bond.does_length_matches_to_bonds(dist_matrix[i][j], bond_definition_list)

                if matches:
                    neighbors.append(j)
            atom.neighbors = neighbors
