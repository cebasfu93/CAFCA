#Celular Automata for Fluid Charge Adjustments
import numpy as np
import sys
import os
import matplotlib.pyplot as plt
import shutil
from string import digits
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
atom_info=inp_file[first_line:last_line]
atom_info_split=[]

#Define arreglos para tipos, coordenadas, cargas y nombres, y los asigna
atom_types=[]
atom_coor=np.zeros((N_atoms,3))
atom_charge=np.zeros(N_atoms)
atom_name=[]

for i in range(N_atoms):
    atom_info_split.append(atom_info[i].split())

    atom_types.append(atom_info_split[i][1])
    atom_coor[i,0]=atom_info_split[i][2]
    atom_coor[i,1]=atom_info_split[i][3]
    atom_coor[i,2]=atom_info_split[i][4]
    atom_name.append(atom_info_split[i][5])
    atom_charge[i]=atom_info_split[i][-1]

#Centra coordenadas en cero. Tambien escribe constantes a archivo
cell_length=0.0

constants=open('constants.outpy', 'a')
for i in range(3):
    atom_coor[:,i]=atom_coor[:,i]-atom_coor[:,i].mean()
    var_max=atom_coor[:,i].max()
    var_min=atom_coor[:,i].min()
    constants.write(str(var_min)+ " " + str(var_max) + "\n")
    try_length=var_max-var_min
    if try_length > cell_length:
        cell_length=try_length

#atom_coor=atom_coor*10/cell_length

constants.write(str(cell_length)+ " empty \n")
constants.close()

#Escribe archivos con nombres, tipos, coordenadas y cargas
names = open('names.outpy', "a")
coords = open('coords.outpy', "a")
types = open('types.outpy', "a")
charges = open('charges.outpy', "a")

for i in range(N_atoms):
    types.write(atom_types[i]+'\n')
    coords.write(str(atom_coor[i,0]) + ' ' + str(atom_coor[i,1]) + ' ' + str(atom_coor[i,2])+'\n')
    names.write(atom_name[i]+'\n')
    charges.write(str(atom_charge[i])+'\n')

names.close()
coords.close()
types.close()
charges.close()
