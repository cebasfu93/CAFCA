from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import gridspec

Nx=100
Ny=150
Nz=100

"""data=np.transpose(np.genfromtxt('fftw_test.txt'))
data_out=np.transpose(np.genfromtxt('fftw_test_out.txt'))
print np.shape(data)
print np.shape(data_out)
xs=np.linspace(0.0, 1.0, 100)
ys=np.linspace(0.0, 2.5, 150)
gx, gy = np.meshgrid(xs, ys)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ax.plot_wireframe(gx, gy, data_out, color='r', rstride=10, cstride=10, label='out')
ax.plot_wireframe(gx, gy, data, rstride=10, cstride=10, label='init')
plt.legend()
plt.show()
"""

arrx=np.genfromtxt('arrx.outc')
arry=np.genfromtxt('arry.outc')
arrz=np.genfromtxt('arrz.outc')
x=np.linspace(0,Nx-1, Nx)
y=np.linspace(0,Ny-1, Ny)
z=np.linspace(0,Nz-1, Nz)

fig=plt.figure()
gs=gridspec.GridSpec(1,3)
ax0=fig.add_subplot(gs[0,0])
ax0.plot(x, arrx)
plt.ylabel(r'arrx')
plt.xlabel('x')

ax1=fig.add_subplot(gs[0,1])
ax1.plot(y, arry)
plt.ylabel(r'arry')
plt.xlabel('y')

ax2=fig.add_subplot(gs[0,2])
ax2.plot(z, arrz)
plt.ylabel(r'arrz')
plt.xlabel('z')

plt.savefig('arr.pdf', format='pdf')
plt.close()
