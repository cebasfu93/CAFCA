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
