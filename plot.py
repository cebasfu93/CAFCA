import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

#importa info de los atomos
atoms=np.genfromtxt('atomos.outc', skip_header=1)

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
