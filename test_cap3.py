import numpy as np
from ipdc import *

G = 6.674e-11  # Unidades N*m^2/kg^2
eps = 0.05      #softening 
N = 2

e = 0.017     # excentricidad de la orbita
a = 1.496e11        # semieje mayor
p = a*(1-e)     #distancia al pericentro

m = np.array([2e30, 6e24]) #masas

mu = G *(m[0]+m[1])
nn = np.sqrt(mu/a**3)

vp = nn * a * np.sqrt((1+e)/(1-e)) #velocidad en el pericentro

apo = a*(1+e)
b = a * np.sqrt(1-e**2)

T = 2*np.pi/nn #periodo


#posiciones y velocidades iniciales
pos = np.array([[0.,0.,0.],[p,0.,0.]])
vel = np.array([[0.,0.,0.],[0.,vp,0.]])


#Aca calculamos las aceleraciones con fortran
ax,ay,az = integradores.aceleracion(eps,m,pos[:,0],pos[:,1],pos[:,2])


#Aca definimos el paso y el nro de pasos
tiempo = 1.
#dt = [10**-6.5, 1e-6, 10**-5.5, 1e-5, 10**-4.5, 1e-4, 1e-3, 10**-2.5, 1e-2]  #pasos
dt = [1e-1, 1e-2]
archivo = open('/home/ornela/Materias/numerico/mn2/_data/cap3/salida.dat','a')

for i in range(0,2):
    
    nit = np.int_(tiempo/dt[i]) #numero de pasos
    x,y,z,vx,vy,vz,ax,ay,az = integradores.euler(eps,dt[i],nit,m,pos[:,0],pos[:,1],pos[:,2],vel[:,0],vel[:,1],vel[:,2],ax,ay,az)

    archivo.write(str('%12.6f'%x)  + '\t'+
                  str('%12.6f'%y)  + '\t'+
                  str('%12.6f'%z)  + '\t'+
                  str('%12.6f'%vx) + '\t'+
                  str('%12.6f'%vy) + '\t'+
                  str('%12.6f'%vz) + '\t'+
                  str('%12.6f'%ax) + '\t'+
                  str('%12.6f'%ay) + '\t'+
                  str('%12.6f'%az) + '\n')
    archivo.close()