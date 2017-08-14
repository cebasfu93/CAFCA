import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import imageio #Para hacer el gif
import os #Para crear directorio temporal
import shutil #Para eliminar directorio temporal
from matplotlib import gridspec
import numpy as np
from matplotlib import rc
rc('font',**{'family':'serif','serif':['Times New Roman']})
#https://matplotlib.org/users/usetex.html
#rc('text', usetex=True)

Z=18
cons=np.genfromtxt('constants.outc', dtype='int')
rspace=np.genfromtxt('rspace.outc')
densx=np.genfromtxt('densx.outc')
densy=np.genfromtxt('densy.outc')
densz=np.genfromtxt('densz.outc')
potx=np.genfromtxt('potx.outc')
poty=np.genfromtxt('poty.outc')
potz=np.genfromtxt('potz.outc')
accx=np.genfromtxt('accx.outc')
accy=np.genfromtxt('accy.outc')
accz=np.genfromtxt('accz.outc')
Lx=cons[0]
Ly=cons[1]
Lz=cons[2]
N_steps=cons[3]
x,y,z= np.linspace(0,Lx-1,Lx), np.linspace(0,Ly-1,Ly), np.linspace(0,Lz-1,Lz)

os.mkdir("temp")

with imageio.get_writer('./CAFCA.gif', mode='I') as writer:
    for i in range(N_steps):
        nx=rspace[i*N_steps:(i+1)*N_steps,0]
        ny=rspace[i*N_steps:(i+1)*N_steps,1]
        nz=rspace[i*N_steps:(i+1)*N_steps,2]
        qs=rspace[i*N_steps:(i+1)*N_steps,3]

        fig=plt.figure(figsize=(24,12))
        gs=gridspec.GridSpec(3,6)

        plt.suptitle("CAFCA", fontsize=Z+4)

        ax0=fig.add_subplot(gs[0,0])
        ax0.plot(x, densx[i*Lx:(i+1)*Lx])
        plt.ylabel(r'Density', fontsize=Z)
        plt.ylim((np.min(densx),np.max(densx)))
        plt.xlim((np.min(x),np.max(x)))
        plt.xticks(fontsize=Z)
        plt.yticks(fontsize=Z)

        ax1=fig.add_subplot(gs[0,1])
        ax1.plot(y, densy[i*Ly:(i+1)*Ly])
        plt.ylim((np.min(densy),np.max(densy)))
        plt.xlim((np.min(y),np.max(y)))
        plt.xticks(fontsize=Z)
        plt.yticks(fontsize=Z)

        ax2=fig.add_subplot(gs[0,2])
        ax2.plot(z, densz[i*Lz:(i+1)*Lz])
        plt.ylim((np.min(densz),np.max(densz)))
        plt.xlim((np.min(z),np.max(z)))
        plt.xticks(fontsize=Z)
        plt.yticks(fontsize=Z)

        ax3=fig.add_subplot(gs[1,0])
        ax3.plot(x, potx[i*Lx:(i+1)*Lx])
        plt.ylabel(r'Potential', fontsize=Z)
        plt.ylim((np.min(potx),np.max(potx)))
        plt.xlim((np.min(x),np.max(x)))
        plt.xticks(fontsize=Z)
        plt.yticks(fontsize=Z)

        ax4=fig.add_subplot(gs[1,1])
        ax4.plot(y, poty[i*Ly:(i+1)*Ly])
        plt.ylim((np.min(poty),np.max(poty)))
        plt.xlim((np.min(y),np.max(y)))
        plt.xticks(fontsize=Z)
        plt.yticks(fontsize=Z)

        ax5=fig.add_subplot(gs[1,2])
        ax5.plot(z, potz[i*Lz:(i+1)*Lz])
        plt.ylim((np.min(potz),np.max(potz)))
        plt.xlim((np.min(z),np.max(z)))
        plt.xticks(fontsize=Z)
        plt.yticks(fontsize=Z)

        ax6=fig.add_subplot(gs[2,0])
        ax6.plot(x, accx[i*Lx:(i+1)*Lx])
        plt.ylabel(r'Acceleration', fontsize=Z)
        plt.xlabel('x', fontsize=Z)
        plt.ylim((np.min(accx),np.max(accx)))
        plt.xlim((np.min(x),np.max(x)))
        plt.xticks(fontsize=Z)
        plt.yticks(fontsize=Z)

        ax6=fig.add_subplot(gs[2,1])
        ax6.plot(y, accy[i*Ly:(i+1)*Ly])
        plt.xlabel('y', fontsize=Z)
        plt.ylim((np.min(accy),np.max(accy)))
        plt.xlim((np.min(y),np.max(y)))
        plt.xticks(fontsize=Z)
        plt.yticks(fontsize=Z)

        ax6=fig.add_subplot(gs[2,2])
        ax6.plot(z, accz[i*Lz:(i+1)*Lz])
        plt.xlabel('z', fontsize=Z)
        plt.ylim((np.min(accz),np.max(accz)))
        plt.xlim((np.min(z),np.max(z)))
        plt.xticks(fontsize=Z)
        plt.yticks(fontsize=Z)

        ax7=fig.add_subplot(gs[:,3:])
        colors=ax7.scatter(nx, ny, nz, c=qs, cmap='seismic', vmin=min(qs), vmax=max(qs[qs<0.8]))
        plt.xlabel('x', fontsize=Z)
        plt.ylabel('y', fontsize=Z)
        plt.xticks(fontsize=Z)
        plt.yticks(fontsize=Z)
        #plt.zlabel('z')
        ax7.set_aspect('equal')
        cbar=fig.colorbar(colors)
        cbar.ax.tick_params(labelsize=Z)

        gs.update(wspace=0.5, hspace=0.5)
        fig = plt.gcf()
        plt.savefig('./temp/test'+str(i)+'.png', format='png')
        plt.close()

        image=imageio.imread('./temp/test'+str(i)+'.png')
        writer.append_data(image)

def ndx2xyz(ndx):
    x = (ndx % (Lx*Ly)) % Lx
    y = (ndx % (Lx*Ly)) // Lx
    z = ndx // (Lx*Ly)
    return x, y, z
