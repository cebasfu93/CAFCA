import numpy as np
from init import dx

#Define lista con nombres de elementos admitidos y sus radios de van der waals
el_names=['H', 'He',
'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne',
'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar',
'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr']
radii=np.array([1.10, 1.43,
2.63, 2.23, 2.05, 1.96, 1.79, 1.71, 1.65, 1.58,
2.77, 2.42, 2.40, 2.26, 2.14, 2.06, 2.05, 1.94,
3.02, 2.78, 2.62, 2.44, 2.27, 2.23, 2.25, 2.27, 2.25, 2.23, 2.27, 2.24, 2.41, 2.32, 2.25, 2.18, 2.10, 2.07])
radii=np.round(radii/(10*dx),0).astype(int)

#Import archivo con numero de electrones para radios dados
file_elec=np.genfromtxt('SphericalGrid/N_elec.outpy').astype(int)

#Define diccionario que le da numero atomico, radio de vdw y color a cada elemento
atomic_numbers={}
el_colors={}
vdw_radii={}
N_bod={}
for i in range(len(el_names)):
    temp_el=el_names[i]
    atomic_numbers[temp_el]=float(i+1)
    vdw_radii[temp_el]=radii[i]
    el_colors[temp_el]='c'
    N_bod[temp_el]=file_elec[vdw_radii[temp_el]-1,1]
el_colors['H']='w'
el_colors['O']='r'
el_colors['N']='b'
el_colors['Cl']='g'

#Crea clase atomo con toda la info de un atomo
class Atom:
    def __init__(self, x , y, z, vx, vy, vz, name, tipo, el):
        self.name = name
        self.type = tipo
        self.el = el
        self.x = x
        self.y = y
        self.z = z
        self.vx = vx
        self.vy = vy
        self.vz = vz
        self.q = atomic_numbers[self.el]
        self.rvdw=vdw_radii[self.el]
        self.N_elec=N_bod[self.el]
        self.qpercube=-round(self.q/self.N_elec, 8)
        self.info = map(str,[self.el, self.name, self.type, self.x, self.y, self.z, self.vx, self.vy, self.vz, self.q, self.rvdw, self.N_elec, self.qpercube])
