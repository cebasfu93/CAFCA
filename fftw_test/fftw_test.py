from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
import numpy as np

data=np.transpose(np.genfromtxt('fftw_test.txt'))
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
