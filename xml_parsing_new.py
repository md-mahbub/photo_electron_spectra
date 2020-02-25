import os
from xml.etree import ElementTree
import os
import copy

FILE_NAME = 'benzen.xml'
TRANSITION_NAME = '[ Reference -- EOM-IP-CCSD state  1/Ag   ]'

carbon_orbital_count = 9
hydrogen_orbital_count = 2
oxygen_orbital_count = 0
nitrogen_orbital_count = 0
n_of_atoms = 12
n_of_basis_functions = 66

full_file = os.path.abspath(os.path.join('data', FILE_NAME))

dom = ElementTree.parse(full_file)

new_dom = copy.deepcopy(dom)
children_of_dyson_orbital = dom.findall('dyson_molecular_orbitals/DMO')

#ref_transition = new_dom.findall('dyson_molecular_orbitals/DMO[@transition="' + TRANSITION_NAME + '"]')

geometry = dom.findall('geometry')
atom_str = (geometry[0]).get('text') # its expected to have only a single <geometry /> tag in benzen.xml
atom_geometry = atom_str.split()
geometric_atom_list = []

for atom in atom_geometry:
    if atom == 'C':
        geometric_atom_list.append(atom)
    elif atom == 'H':
        geometric_atom_list.append(atom)
    elif atom == 'O':
        geometric_atom_list.append(atom)

offset = 0

for i in range(n_of_atoms):
    new_dom = copy.deepcopy(dom)
    ref_transition = new_dom.findall('dyson_molecular_orbitals/DMO[@transition="' + TRANSITION_NAME + '"]')

    if geometric_atom_list[i] == 'C':
        orbital_count = carbon_orbital_count
    elif geometric_atom_list[i] == 'H':
        orbital_count = hydrogen_orbital_count
    elif geometric_atom_list[i] == 'O':
        orbital_count = oxygen_orbital_count
    elif geometric_atom_list[i] == 'N':
        orbital_count = nitrogen_orbital_count

    minIndex = offset
    maxIndex = offset + orbital_count - 1;
    offset = maxIndex + 1;

    modified_basis_func_right = []
    modified_basis_func_left = []


    for dyson in ref_transition:
        if dyson.get('comment') == 'dyson right':
            basis_func_string = dyson.get('text')
            basis_func_list = basis_func_string.split()

            modified_basis_func_right = ['0'] * n_of_basis_functions  # taking list of zeros
            modified_basis_func_right[minIndex: maxIndex + 1] = basis_func_list[minIndex: maxIndex + 1]
            dyson.set('text', ' '.join(modified_basis_func_right))

        elif dyson.get('comment') == 'dyson left':
            basis_func_string = dyson.get('text')
            basis_func_list = basis_func_string.split()

            modified_basis_func_left = ['0'] * n_of_basis_functions  # taking list of zeros
            modified_basis_func_left[minIndex: maxIndex + 1] = basis_func_list[minIndex: maxIndex + 1]
            dyson.set('text', ' '.join(modified_basis_func_left))

    directory = './Data/' + geometric_atom_list[i] + str(i+1)
    filepath = os.path.join(directory, FILE_NAME)

    try:
        os.makedirs(directory, exist_ok=True)
        new_dom.write(filepath)
    except FileExistsError:
        print('directory already exists')
    except Exception as e:
        print(e)





