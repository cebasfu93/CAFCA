#Celular Automata for Fluid Charge Adjustments
from atom import Atom
import numpy as np
import matplotlib.pyplot as plt
from operator import attrgetter
from optparse import OptionParser

#Opciones del script
parser=OptionParser()
parser.add_option("-i","--input", action="store", type="string", dest="InputFile", help="Path of the .mol2 input file.")
parser.add_option("-L","--realdistance", action="store", type="float", dest="Distance", default=5.0, help="Real-space distance from atoms to simulation box (A).")
parser.add_option("-V","--veldistance", action="store", type="float", dest="Distancev", default=5.0, help="Velocity distance from atoms to simulation box (A/ps).")
parser.add_option("-N","--resolution", action="store", type="int", dest="Resolution", default=512, help="Resolution of the phase space (pixels)")
parser.add_option("-Q","--charge", action="store", type="int", dest="ChargeDiscrete", default=0.04, help="Charge discretization (Q0)")

(options, args) = parser.parse_args()

inp_name=options.InputFile
distx=options.Distance
distv=options.Distancev
N_res=options.Resolution
q_fund=options.ChargeDiscrete

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
        x=float(atom_info_split[i][2])
        y=float(atom_info_split[i][3])
        z=float(atom_info_split[i][4])
        vx=float(0.0)
        vy=float(0.0)
        vz=float(0.0)
        tipo=atom_info_split[i][5]
        q=float(atom_info_split[i][-1])

        atoms_list.append(Atom(x,y,z, vx, vy, vz, q, name, tipo))
    return atoms_list
#Escribe archivos con nombres, tipos, coordenadas, velocidades y cargas
def save_molecule(atoms_list):
    names = open('names.outpy', "a")
    coords = open('coords.outpy', "a")
    speeds = open('speeds.outpy', "a")
    types = open('types.outpy', "a")
    charges = open('charges.outpy', "a")

    for i in range(N_atoms-1):
        types.write(atoms_list[i].type+'\n')
        coords.write(str(atoms_list[i].x) + ' ' + str(atoms_list[i].y) + ' ' + str(atoms_list[i].z)+'\n')
        speeds.write(str(atoms_list[i].vx) + ' ' + str(atoms_list[i].vy) + ' ' + str(atoms_list[i].vz)+'\n')
        names.write(atoms_list[i].name+'\n')
        charges.write(str(atoms_list[i].q)+'\n')

    names.close()
    coords.close()
    speeds.close()
    types.close()
    charges.close()
#Coge el objeto con el atributo minimo de una molecula
def take_min(atoms_list, atribute):
    return getattr(min(atoms_list, key=attrgetter(atribute)), atribute)
#Coge el objeto con el atributo maximo de una molecula
def take_max(atoms_list, atribute):
    return getattr(max(atoms_list, key=attrgetter(atribute)), atribute)
#Imprime constantes para el programa
def save_cons(atoms_list):
    Lx_min, Lx_max = take_min(atoms_list, 'x')-distx, take_max(atoms_list, 'x')+distx
    Ly_min, Ly_max = take_min(atoms_list, 'y')-distx, take_max(atoms_list, 'y')+distx
    Lz_min, Lz_max = take_min(atoms_list, 'z')-distx, take_max(atoms_list, 'z')+distx
    Vx_min, Vx_max = take_min(atoms_list, 'vx')-distv, take_max(atoms_list, 'vx')+distv
    Vy_min, Vy_max = take_min(atoms_list, 'vy')-distv, take_max(atoms_list, 'vy')+distv
    Vz_min, Vz_max = take_min(atoms_list, 'vz')-distv, take_max(atoms_list, 'vz')+distv
    cons=open('constants.outpy', "a")
    cons.write(str(Lx_min) + ' ' + str(Lx_max) + '\n')
    cons.write(str(Ly_min) + ' ' + str(Ly_max) + '\n')
    cons.write(str(Lz_min) + ' ' + str(Lz_max) + '\n')
    cons.write(str(Vx_min) + ' ' + str(Vx_max) + '\n')
    cons.write(str(Vy_min) + ' ' + str(Vy_max) + '\n')
    cons.write(str(Vz_min) + ' ' + str(Vz_max) + '\n')
    cons.write(str(N_res) + '\n')
    cons.write(str(q_fund) + '\n')
    cons.close()

molecule=init_molecule(inp_file)
save_molecule(molecule)
save_cons(molecule)
