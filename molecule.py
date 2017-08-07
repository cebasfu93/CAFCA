import init
import numpy as np
import re #Para get_element
import numpy as np
from atom import *

#importa la informacion de los atomos en el .mol2
inp_file=np.genfromtxt(init.inp_name, delimiter='\n', dtype='string')
for i in range(len(inp_file)):
    if "@<TRIPOS>ATOM" in inp_file[i]:
        first_line=i+1
    if "@<TRIPOS>BOND" in inp_file[i]:
        last_line=i
        break
N_atoms=last_line - first_line

Lx, Ly, Lz, Vx, Vy, Vz = 0, 0, 0, 0, 0, 0

#Convierte coordenadas a enteros poniendo el minimo en cero
def int_coordinates(at_info_split):
    temp_array=np.array(at_info_split, dtype='string')
    temp_array=temp_array[:,2:5].astype(float)
    for i in range(3):
        min_temp=np.min(temp_array[:,i])
        temp_array[:,i]=temp_array[:,i]-(min_temp-init.distx)

    global Lx
    global Ly
    global Lz
    global Vx
    global Vy
    global Vz

    Lx, Ly, Lz =int(np.round(np.max((temp_array[:,0])+init.distx)/init.dx,0)), int(np.round(np.max((temp_array[:,1])+init.distx)/init.dx,0)), int(np.round(np.max((temp_array[:,2])+init.distx)/init.dx,0))
    Vx, Vy, Vz = int(round(2*init.distv,0)), int(round(2*init.distv,0)), int(round(2*init.distv,0))

    temp_array=np.round(temp_array/init.dx,0).astype(int)

    return temp_array

#Define lista con elementos tipo atomo e inicializa sus atributos
def init_molecule(input_file):
    atoms_list=[]

    atom_info=input_file[first_line:last_line]
    atom_info_split=[]
    for i in range(N_atoms):
        atom_info_split.append(atom_info[i].split())

    integer_coords=int_coordinates(atom_info_split)

    for i in range(N_atoms):
        name=atom_info_split[i][1]
        element=get_element(name)
        x=integer_coords[i,0]
        y=integer_coords[i,1]
        z=integer_coords[i,2]
        vx=(0)
        vy=(0)
        vz=(0)
        tipo=atom_info_split[i][5]

        atoms_list.append(Atom(x,y,z, vx, vy, vz, name, tipo, element))
    return atoms_list

#Dado el nombre de un atomo, asigna su elemento
def get_element(nombre):
    match = re.match(r"([a-zA-Z]+)([0-9]+)", nombre, re.I)
    if match:
        items = match.groups()
    return items[0].title()

#Escribe archivo con nombres, tipos, coordenadas, velocidades y cargas
def save_molecule(at_list):
    atomos = open('atomos.outpy', "w")
    all_at=[]
    for i in range(N_atoms):
        all_at.append(at_list[i].info)
    col_width = max(len(word) for row in all_at for word in row) + 2
    for row in all_at:
        atomos.write("".join(word.ljust(col_width) for word in row) + "\n")
    atomos.close()

#Imprime constantes para el programa
def save_cons():
    cons=open('constants.outpy', "w")
    cons.write(str(Lx) + '\n')
    cons.write(str(Ly) + '\n')
    cons.write(str(Lz) + '\n')
    cons.write(str(Vx) + '\n')
    cons.write(str(Vy) + '\n')
    cons.write(str(Vz) + '\n')
    cons.write(str(N_atoms) + '\n')
    cons.close()

molecule=init_molecule(inp_file)
save_molecule(molecule)
save_cons()
