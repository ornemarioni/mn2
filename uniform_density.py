import numpy as np
import random

# Este programa nos devuelve una distribucion uniforme de particulas
# dentro de una esfera de radio 1

def dens(N):     
        
    range_r3  = np.linspace(0, 1/3., 2*N)
    range_ct  = np.linspace(-1, 1, 2*N)
    range_phi = np.linspace(0, 2*np.pi, 2*N)
    
    r3  = np.random.choice(range_r3,  N, replace = True) #r**3
    ct  = np.random.choice(range_ct,  N, replace = True) #cos(t)
    phi = np.random.choice(range_phi, N, replace = True) #phi
    
    st = np.sqrt(1. - ct**2) #sen(t)
        
    x = (r3*3)**(1/3.) * st * np.cos(phi)
    y = (r3*3)**(1/3.) * st * np.sin(phi)
    z = (r3*3)**(1/3.) * ct
    
    return x,y,z