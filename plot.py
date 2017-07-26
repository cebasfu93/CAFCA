import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

def scatter3d(xs, ys, zs, color, plot_name, xlabel='X', ylabel='Y', zlabel='Z'):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    colors=ax.scatter(xs, ys, zs, c=color, cmap='seismic')
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_zlabel(zlabel)
    ax.set_aspect('equal')
    #ax.set_aspect('equal', 'datalim')
    fig.colorbar(colors)
    plt.savefig(plot_name+'.pdf', format='pdf')
    #plt.show()
    plt.close()

#importa info de los atomos
atoms=np.genfromtxt('atomos.outc', skip_header=1)
norm_vel=np.linalg.norm(atoms[:,3:6], axis=1)

scatter3d(atoms[:,0], atoms[:,1], atoms[:,2], atoms[:,6], 'Initial_realq')
scatter3d(atoms[:,3], atoms[:,4], atoms[:,5], atoms[:,6], 'Initial_velsq', 'VX', 'VY', 'VZ')
scatter3d(atoms[:,0], atoms[:,1], atoms[:,2], norm_vel, 'Initial_realv')
