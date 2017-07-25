#Celular Automata for Fluid Charge Adjustments
from atom import Atom
import numpy as np
import matplotlib.pyplot as plt
from optparse import OptionParser

#Opciones del script
parser=OptionParser()
parser.add_option("-i","--input", action="store", type="string", dest="InputFile", help="Path of the .mol2 input file.")
parser.add_option("-d","--distance", action="store", type="float", dest="Distance", default=10, help="Distance from atoms to simulation box (A).")

(options, args) = parser.parse_args()

inp_name=options.InputFile
dist=options.Distance

#importa la informacion de los atomos en el .mol2
inp_file=np.genfromtxt(inp_name, delimiter='\n', dtype='string')
for i in range(len(inp_file)):
    if "@<TRIPOS>ATOM" in inp_file[i]:
        first_line=i+1
    if "@<TRIPOS>BOND" in inp_file[i]:
        last_line=i
        break
N_atoms=last_line - first_line

#Define lista con elementos tipo atomo e inicializa sus atributos
def init_molecule(input_file):
    atoms_list=[]

    atom_info=input_file[first_line:last_line]
    atom_info_split=[]
    for i in range(N_atoms):
        atom_info_split.append(atom_info[i].split())

        name=atom_info_split[i][1]
        x=atom_info_split[i][2]
        y=atom_info_split[i][3]
        z=atom_info_split[i][4]
        tipo=atom_info_split[i][5]
        q=atom_info_split[i][-1]

        atoms_list.append(Atom(x,y,z,q,name,tipo))
    return atoms_list
#Escribe archivos con nombres, tipos, coordenadas y cargas
def save_molecule(atoms_list):
    names = open('names.outpy', "a")
    coords = open('coords.outpy', "a")
    types = open('types.outpy', "a")
    charges = open('charges.outpy', "a")

    for i in range(N_atoms-1):
        types.write(atoms_list[i].type+'\n')
        coords.write(str(atoms_list[i].x) + ' ' + str(atoms_list[i].y) + ' ' + str(atoms_list[i].z)+'\n')
        names.write(atoms_list[i].name+'\n')
        charges.write(str(atoms_list[i].q)+'\n')

    names.close()
    coords.close()
    types.close()
    charges.close()

molecule=init_molecule(inp_file)
save_molecule(molecule)
