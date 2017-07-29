import numpy as np

class Atom:
    def __init__(self, x , y, z, vx, vy, vz, q, name, tipo):
        self.x = x
        self.y = y
        self.z = z
        self.vx = vx
        self.vy = vy
        self.vz = vz
        self.q = q
        self.name = name
        self.type = tipo

atomic_numbers={}

el_names=['H', 'He',
'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne',
'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar',
'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr']
for i in range(len(el_names)):
    atomic_numbers[el_names[i]]=i+1

print atomic_numbers
