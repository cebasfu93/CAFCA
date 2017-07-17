#Celular Automata for Fluid Charge Adjustments
import numpy as np
import sys
import os
import matplotlib.pyplot as plt
import imageio
import shutil
from string import digits

#Hace diccionario con radios vdw (angs)
vdw={}
vdw['C']=1.7 #Angs
vdw['CL']=1.75
vdw['H']=1.2

#importa la informacion de los atomos en el .mol2
inp=np.genfromtxt(sys.argv[1], delimiter='\n', dtype='string')
for i in range(len(inp)):
    if "@<TRIPOS>ATOM" in inp[i]:
        ini=i+1
    if "@<TRIPOS>BOND" in inp[i]:
        fin=i
        break

g=inp[ini:fin]
N=len(g) #numero de atomos
L=0.02 #ancho de un slot (Angs)
dis_caja=3 #distancia (angs) de la molecula a la caja
decay=2 #angs

elementos=[]
coordx=np.zeros(N)
coordy=np.zeros(N)
coordz=np.zeros(N)
cargas=np.zeros(N)

#Hace diccionario con radios vdw (slots)
qvdw={}
for key, value in vdw.items():
    qvdw[key]=value/L

#Hace arrays con todas las coordenadas x, y y z (angs)
for i in range(N):
    g_temp=g[i].split()
    elementos.append(g_temp[1].translate(None, digits).upper())
    coordx[i]=g_temp[2]
    coordy[i]=g_temp[3]
    coordz[i]=g_temp[4]
    cargas[i]=g_temp[8]

#centra la molecula en 0A
x_new=coordx-coordx.mean()
y_new=coordy-coordy.mean()
z_new=coordz-coordz.mean()

#Hace las dimensiones de la caja (angs)
longx=(np.abs(x_new).max()+dis_caja)*2
longy=(np.abs(y_new).max()+dis_caja)*2
longz=(np.abs(z_new).max()+dis_caja)*2

#pone un vertice de la caja en (0,0,0)
x_fin=x_new+longx/2
y_fin=y_new+longy/2
z_fin=z_new+longz/2

#Hace dimensiones de la caja (slots)
dimx=int(longx/L)
dimy=int(longy/L)
dimz=int(longz/L)

#crea el espacio
space=np.zeros((dimx, dimy, dimz))

for i in range(N):
    R=vdw[elementos[i]]

    nucx=int(x_fin[i]/L)
    nucy=int(y_fin[i]/L)
    nucz=int(z_fin[i]/L)

    decslot=int(decay/L)+1
    radslot=np.zeros((decslot,decslot,decslot))
    radangs=np.zeros((decslot,decslot,decslot))

    for p in range(0,decslot):
        for q in range(0,decslot):
            for r in range(0,decslot):
                rad_temp=np.sqrt(p**2+q**2+r**2)
                radslot[p,q,r]=rad_temp
                radangs[p,q,r]=L*np.round(radslot[p,q,r])
    norma_temp=0
    for j in range(0,len(radangs[:,0,0])-1):
        delt=radangs[j+1,0,0]-radangs[j,0,0]
        norma_temp+=np.exp(-radangs[j,0,0]**2/(2*R**2))*radangs[j,0,0]**2*delt
    norma=1/(norma_temp*4*np.pi)

    for p in range(-decslot+1,decslot):
        for q in range(-decslot+1,decslot):
            for r in range(-decslot+1,decslot):
                if radangs[np.abs(p),np.abs(q),np.abs(r)]<=2.0:
                    space[nucx+p,nucy+q,nucz+r]+=norma*np.exp(-radangs[np.abs(p),np.abs(q),np.abs(r)]**2/(2*R**2))*cargas[i]

#Crea un gif con planos de diferentes valores de z
minimo=space.min()
maximo=space.max()
with imageio.get_writer('./Ejex.gif', mode='I') as writer:
    os.mkdir('temp')
    for i in range(dimx):
        if i%(dimx/60)==0:
            fig=plt.figure()
            ax=plt.axes()
            plt.title(str(i))
            im = plt.imshow(space[i,:,:], vmin=minimo,vmax=maximo)
            fig.colorbar(im)
            plt.savefig('./temp/espacio'+str(i)+'.png', format='png')
            plt.close()

            image=imageio.imread('./temp/espacio'+str(i)+'.png')
            writer.append_data(image)
shutil.rmtree('temp')
with imageio.get_writer('./Ejez.gif', mode='I') as writer:
    os.mkdir('temp')
    for i in range(dimz):
        if i%(dimz/60)==0:
            fig=plt.figure()
            ax=plt.axes()
            plt.title(str(i))
            im = plt.imshow(space[:,:,i], vmin=minimo,vmax=maximo)
            fig.colorbar(im)
            plt.savefig('./temp/espacio'+str(i)+'.png', format='png')
            plt.close()

            image=imageio.imread('./temp/espacio'+str(i)+'.png')
            writer.append_data(image)
shutil.rmtree('temp')
