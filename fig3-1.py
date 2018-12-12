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

T = 2*pi/nn #periodo


#posiciones y velocidades iniciales
pos = np.array([[0.,0.,0.],[p,0.,0.]])
vel = np.array([[0.,0.,0.],[0.,vp,0.]])


#Aca calculamos las aceleraciones con fortran
ax,ay,az = integradores.aceleracion(eps,m,pos[:,0],pos[:,1],pos[:,2])


#Aca definimos el paso y el nro de pasos
tiempo = 100.
dt = [10**-6.5, 1e-6, 10**-5.5, 1e-5, 10**-4.5, 1e-4, 1e-3, 10**-2.5, 1e-2]  #pasos

for i in range(0,len(dt)):
    
    nit = np.int_(tiempo/dt[i]) #numero de pasos
    integradores.euler(eps,dt[i],nit,m,pos[:,0],pos[:,1],pos[:,2],vel[:,0],vel[:,1],vel[:,2],ax,ay,az)