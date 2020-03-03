import os
from xml.etree import ElementTree as ET
#from lxml import etre
from xml.dom import minidom


import os
import errno
import copy
from Atom import *
from Bond import *
from MoleculeStruct import *


def indent(elem, level=0):
    i = "\n" + level * "  "
    if len(elem):
        if not elem.text or not elem.text.strip():
            elem.text = i + "  "
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
        for elem in elem:
            indent(elem, level + 1)
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
    else:
        if level and (not elem.tail or not elem.tail.strip()):
            elem.tail = i

FILE_NAME = 'benzen.xml'
XML_FILE_PATH = './data/'
TRANSITION_NAME = '[ Reference -- EOM-IP-CCSD state  1/Ag   ]'

CARBON_COEFF_COUNT = 9
HYDROGEN_COEFF_COUNT = 2
OXYGEN_COEFF_COUNT = 0
NITROGEN_COEFF_COUNT = 0

'''
# not implemented in code
# unit is Angstrom(A), 1A = 10^(-10) meter
C_H_SINGLE_BOND_DIST = 1.09 # 1.04 ~ 1.14
C_C_SINGLE_BOND_DIST = 1.54 # 1.49 ~ 1.59
C_C_DOUBLE_BOND_DIST = 1.34 # 1.29 ~ 1.39
C_C_TRIPPLE_BOND_DIST = 1.20  # 1.15 ~ 1.25
# how about overlapping distance?
'''

# BOND_DEFINITION_LIST holds all possible bond types that belongs to the molecule.
# add/change the numbers in this list if your atoms have different bond lengths, count, etc.
BOND_DEFINITION_LIST = [
    Bond().load('C', 'H', 1, 1.040, 1.149),
    Bond().load('C', 'C', 1, 1.490, 1.599),
    Bond().load('C', 'C', 2, 1.290, 1.399),
    Bond().load('C', 'C', 3, 1.150, 1.259)]

# ATOM_DEFINITION_LIST holds all possible atom types that belongs to the molecule.
# add/change the list elements if your molecule have different atoms and possible inhome neighbors
#TODO: dynamically load the coefficient counts based on 'qchem_input->BASIS'
ATOM_DEFINITION_MAP = {
    'C': Atom().load('C', 9, ['H']),
    'H': Atom().load('H', 2, ['C'])
   }

#MOLECULE_STRUCT = MoleculeStruct().load(ATOM_DEFINITION_MAP, BOND_DEFINITION_LIST)
molecule_structure = MoleculeStruct()
molecule_structure.load(ATOM_DEFINITION_MAP, BOND_DEFINITION_LIST)

full_file = os.path.abspath(os.path.join('data', FILE_NAME))
tree = ET.parse(FILE_NAME)
root = tree.getroot()
geometry = tree.findall('geometry')
atom_list_str = (geometry[0]).get('text')
atoms = []
atoms = Atom.load_atoms_from_geometry(atom_list_str)

Atom.assign_neighbors_to_atoms(BOND_DEFINITION_LIST, atoms, molecule_structure)

# shadowing the original dom to make changes on the shadow and to save as output files
original_dyson_transition = tree.findall('dyson_molecular_orbitals/DMO[@transition="' + TRANSITION_NAME + '"]')
basis = tree.findall('basis')
n_of_basis_functions = int((basis[0]).get('n_of_basis_functions'))
offset = 0

for i, atom in enumerate(atoms):
    new_dom = copy.deepcopy(tree)
    original_dyson_transition = new_dom.findall('dyson_molecular_orbitals/DMO[@transition="' + TRANSITION_NAME + '"]')

    # for each atom get all the neighbors
    neighbors = atom.neighbors
    # from the neighbor list only keep inhome neighbors
    possible_inhome_neighbors = (molecule_structure.get_atom_definition_map()[atom.atom_name]).possible_inhome_neighbors
    inhome_neighbors = []
    outhome_neighbors = []

    for neighbor in neighbors:
        if atoms[neighbor].atom_name in possible_inhome_neighbors:
            inhome_neighbors.append(neighbor)
        else:
            outhome_neighbors.append(neighbor)

    # now we are left only inhome neighbors in the 'inhome_neighbors' list

    # find chunks of coefficients in the basis function list,
    # then keep those coefficients and make other's value 0
    for dyson in original_dyson_transition:
        if dyson.get('comment') == 'dyson right':
            basis_func_string = dyson.get('text')
            basis_func_list = basis_func_string.split()

            modified_basis_func_right = ['0'] * n_of_basis_functions  # taking list of zeros
            #for the atom itself avoid making zeros
            modified_basis_func_right[atom.coefficient_start_index: atom.coefficient_end_index] = basis_func_list[atom.coefficient_start_index: atom.coefficient_end_index]
            # now consider the inhome_neighbors
            for inhome_neighbor in inhome_neighbors:
                inhome_neighbor_start_idx = atoms[inhome_neighbor].coefficient_start_index
                inhome_neighbor_end_idx = atoms[inhome_neighbor].coefficient_end_index

                modified_basis_func_right[inhome_neighbor_start_idx: inhome_neighbor_end_idx] \
                    = basis_func_list[inhome_neighbor_start_idx: inhome_neighbor_end_idx + 1]
            dyson.set('text', ' '.join(modified_basis_func_right))

        elif dyson.get('comment') == 'dyson left':
            basis_func_string = dyson.get('text')
            basis_func_list = basis_func_string.split()

            modified_basis_func_left = ['0'] * n_of_basis_functions  # taking list of zeros
            # for the atom itself avoid making zeros
            modified_basis_func_left[atom.coefficient_start_index: atom.coefficient_end_index] = basis_func_list[
                                                                                                  atom.coefficient_start_index: atom.coefficient_end_index]
            # now the inhome_neighbors
            for inhome_neighbor in inhome_neighbors:
                inhome_neighbor_start_idx = atoms[inhome_neighbor].coefficient_start_index
                inhome_neighbor_end_idx = atoms[inhome_neighbor].coefficient_end_index

                modified_basis_func_left[inhome_neighbor_start_idx: inhome_neighbor_end_idx] \
                    = basis_func_list[inhome_neighbor_start_idx: inhome_neighbor_end_idx + 1]
            dyson.set('text', " ".join(modified_basis_func_left))

    directory = './data/' + atom.atom_name + str(i + 1)
    filepath = os.path.join(directory, FILE_NAME)

    try:
        os.makedirs(directory, exist_ok=True)
        new_dom.write(filepath, encoding="utf-8")
        #new_root = new_dom.getroot()
        #ET.tostring(new_dom).toprettyxml(indent="   ",  newl="",  encoding="").write(filepath, encoding="utf-8")
        #xmlstr = minidom.parseString(ET.tostring(new_dom)).toprettyxml(indent="   ",  newl="",  encoding="")
        #print(ET.tostring(new_dom))
        #print(ET.tostring(new_root, encoding='utf8').decode('utf8'))
    except Exception as e:
        print(e)


    # if atom doesn't have out-home neighbor then overlapping is invalid
    if outhome_neighbors:
        # now generate files for right overlap
        dom_overlap = copy.deepcopy(tree)
        dyson_transition_left_overlap = dom_overlap.find('dyson_molecular_orbitals/DMO[@transition="' + TRANSITION_NAME + '"][@comment="dyson left"]')
        dyson_transition_right_overlap = dom_overlap.find('dyson_molecular_orbitals/DMO[@transition="' + TRANSITION_NAME + '"][@comment="dyson right"]')

        basis_func_string = dyson_transition_right_overlap.get('text')
        right_basis_funcs = basis_func_string.split()

        modified_basis_func_right = ['0'] * n_of_basis_functions  # taking list of zeros
        #for the atom itself avoid making zeros
        modified_basis_func_right[atom.coefficient_start_index: atom.coefficient_end_index] = right_basis_funcs[atom.coefficient_start_index: atom.coefficient_end_index]
        # now consider the inhome_neighbors
        for inhome_neighbor in inhome_neighbors:
            inhome_neighbor_start_idx = atoms[inhome_neighbor].coefficient_start_index
            inhome_neighbor_end_idx = atoms[inhome_neighbor].coefficient_end_index

            modified_basis_func_right[inhome_neighbor_start_idx: inhome_neighbor_end_idx] \
                = right_basis_funcs[inhome_neighbor_start_idx: inhome_neighbor_end_idx + 1]
        dyson_transition_right_overlap.set('text', ' '.join(modified_basis_func_right))


        basis_func_string = dyson_transition_left_overlap.get('text')
        left_basis_funcs = basis_func_string.split()

        modified_basis_func_left = ['0'] * n_of_basis_functions  # taking list of zeros
        #for the atom's outhome neighbor itself and the neighbors residents avoid making zeros
        for outhome_neighbor_atom in outhome_neighbors:
            #for the atom's outhome neighbor itself

            modified_basis_func_left[atoms[outhome_neighbor_atom].coefficient_start_index: atoms[outhome_neighbor_atom].coefficient_end_index] \
                = left_basis_funcs[atoms[outhome_neighbor_atom].coefficient_start_index: atoms[outhome_neighbor_atom].coefficient_end_index]
            # now consider the neighbor's residents

            # for each atom get all the neighbors
            neighboring_neighbors = atoms[outhome_neighbor_atom].neighbors
            # from the neighbor list only keep inhome neighbors
            possible_inhome_neighbors = (molecule_structure.get_atom_definition_map()[atoms[outhome_neighbor_atom].atom_name]).possible_inhome_neighbors
            neighboring_inhome_neighbors = []

            for neighbor in neighboring_neighbors:
                if atoms[neighbor].atom_name in possible_inhome_neighbors:
                    neighboring_inhome_neighbors.append(neighbor)

            for neighbor in neighboring_inhome_neighbors:
                neighbor_start_idx = atoms[neighbor].coefficient_start_index
                neighbor_end_idx = atoms[neighbor].coefficient_end_index

                modified_basis_func_left[neighbor_start_idx: neighbor_end_idx] \
                    = left_basis_funcs[neighbor_start_idx: neighbor_end_idx + 1]
            dyson_transition_left_overlap.set('text', ' '.join(modified_basis_func_left))

        directory = './data/' + atom.atom_name + str(i + 1) + '_right'
        filepath = os.path.join(directory, FILE_NAME)

        try:
            os.makedirs(directory, exist_ok=True)
            dom_overlap.write(filepath, encoding="utf-8")

        except Exception as e:
            print(e)

        # now generate files for left overlap
        dom_overlap = copy.deepcopy(tree)
        dyson_transition_left_overlap = dom_overlap.find(
            'dyson_molecular_orbitals/DMO[@transition="' + TRANSITION_NAME + '"][@comment="dyson left"]')
        dyson_transition_right_overlap = dom_overlap.find(
            'dyson_molecular_orbitals/DMO[@transition="' + TRANSITION_NAME + '"][@comment="dyson right"]')

        basis_func_string = dyson_transition_left_overlap.get('text')
        left_basis_funcs = basis_func_string.split()

        modified_basis_func_left = ['0'] * n_of_basis_functions  # taking list of zeros
        # for the atom itself avoid making zeros
        modified_basis_func_left[atom.coefficient_start_index: atom.coefficient_end_index] = left_basis_funcs[
                                                                                              atom.coefficient_start_index: atom.coefficient_end_index]
        # now consider the inhome_neighbors
        for inhome_neighbor in inhome_neighbors:
            inhome_neighbor_start_idx = atoms[inhome_neighbor].coefficient_start_index
            inhome_neighbor_end_idx = atoms[inhome_neighbor].coefficient_end_index

            modified_basis_func_left[inhome_neighbor_start_idx: inhome_neighbor_end_idx] \
                = left_basis_funcs[inhome_neighbor_start_idx: inhome_neighbor_end_idx + 1]
        dyson_transition_left_overlap.set('text', ' '.join(modified_basis_func_left))

        basis_func_string = dyson_transition_right_overlap.get('text')
        right_basis_funcs = basis_func_string.split()

        modified_basis_func_right = ['0'] * n_of_basis_functions  # taking list of zeros
        # for the atom's outhome neighbor itself and the neighbors residents avoid making zeros
        for outhome_neighbor_atom in outhome_neighbors:
            # for the atom's outhome neighbor itself

            modified_basis_func_right[atoms[outhome_neighbor_atom].coefficient_start_index: atoms[
                outhome_neighbor_atom].coefficient_end_index] \
                = right_basis_funcs[atoms[outhome_neighbor_atom].coefficient_start_index: atoms[
                outhome_neighbor_atom].coefficient_end_index]
            # now consider the neighbor's residents

            # for each atom get all the neighbors
            neighboring_neighbors = atoms[outhome_neighbor_atom].neighbors
            # from the neighbor list only keep inhome neighbors
            possible_inhome_neighbors = (molecule_structure.get_atom_definition_map()[
                atoms[outhome_neighbor_atom].atom_name]).possible_inhome_neighbors
            neighboring_inhome_neighbors = []

            for neighbor in neighboring_neighbors:
                if atoms[neighbor].atom_name in possible_inhome_neighbors:
                    neighboring_inhome_neighbors.append(neighbor)

            for neighbor in neighboring_inhome_neighbors:
                neighbor_start_idx = atoms[neighbor].coefficient_start_index
                neighbor_end_idx = atoms[neighbor].coefficient_end_index

                modified_basis_func_right[neighbor_start_idx: neighbor_end_idx] \
                    = left_basis_funcs[neighbor_start_idx: neighbor_end_idx + 1]
            dyson_transition_right_overlap.set('text', ' '.join(modified_basis_func_right))

        directory = './data/' + atom.atom_name + str(i + 1) + '_left'
        filepath = os.path.join(directory, FILE_NAME)

        try:
            os.makedirs(directory, exist_ok=True)
            dom_overlap.write(filepath, encoding="utf-8")

        except Exception as e:
            print(e)
