import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

constantes=np.genfromtxt('constants.outpy')
cell_length=constantes[3,0]
x_min, x_max = constantes[0,0], constantes[0,1]
y_min, y_max = constantes[1,0], constantes[1,1]
z_min, z_max = constantes[2,0], constantes[2,1]
margin=1.20

#importa info de los atomos
atoms=np.genfromtxt('Atomos_c.outc', skip_header=1)

#re-escala las coordenadas
atoms[:,0:3]=atoms[:,0:3]*cell_length

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
colors=ax.scatter(atoms[:,0], atoms[:,1], atoms[:,2], c=atoms[:,3], cmap='seismic')
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ax.set_aspect('equal')
#ax.set_aspect('equal', 'datalim')
fig.colorbar(colors)
plt.savefig('Initial.pdf', format='pdf')
plt.show()
