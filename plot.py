import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import pandas as pd
from atom import el_colors

def scatter3d(xs, ys, zs, color, plot_name, xlabel='X', ylabel='Y', zlabel='Z', barra=False, inter=False):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    colors=ax.scatter(xs, ys, zs, c=color, cmap='seismic')
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_zlabel(zlabel)
    ax.set_aspect('equal')
    #ax.set_aspect('equal', 'datalim')
    if barra:
        fig.colorbar(colors)
    plt.savefig(plot_name+'.pdf', format='pdf')
    if inter:
        plt.show()
    plt.close()

#importa info de los atomos
atoms=np.genfromtxt('atomos.outc', skip_header=1, dtype='str', delimiter=', ')
atoms_str=atoms[:,0:3].astype(str)
atoms_float=atoms[:,3:].astype(float)
colores=[el_colors[x] for x in atoms_str[:,0]]
norm_vel=np.linalg.norm(atoms_float[:,3:5], axis=1)

scatter3d(atoms_float[:,0], atoms_float[:,1], atoms_float[:,2], colores, 'Initial_realel', inter=True)
scatter3d(atoms_float[:,0], atoms_float[:,1], atoms_float[:,2], atoms_float[:,-1], 'Initial_realq', barra=True)
scatter3d(atoms_float[:,3], atoms_float[:,4], atoms_float[:,5], atoms_float[:,-1], 'Initial_velsq', 'VX', 'VY', 'VZ', barra=True)
scatter3d(atoms_float[:,0], atoms_float[:,1], atoms_float[:,2], norm_vel, 'Initial_realv', barra=True)
