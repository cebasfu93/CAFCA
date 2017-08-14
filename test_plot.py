import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

cons=np.genfromtxt('constants.outpy')
inp=np.genfromtxt('rspace.outc')
Lx=cons[0]
Ly=cons[1]
Lz=cons[2]

nx=inp[:,0]
ny=inp[:,1]
nz=inp[:,2]
qs=inp[:,3]

def ndx2xyz(ndx):
    x = (ndx % (Lx*Ly)) % Lx
    y = (ndx % (Lx*Ly)) // Lx
    z = ndx // (Lx*Ly)
    return x, y, z

#xs, ys, zs = ndx2xyz(ind)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
colors=ax.scatter(nx, ny, nz, c=qs, cmap='seismic', vmin=min(qs), vmax=max(qs[qs<0.8]))
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
ax.set_aspect('equal')
fig.colorbar(colors)
plt.show()
#plt.close()
