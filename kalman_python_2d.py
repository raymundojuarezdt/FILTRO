# -*- coding: utf-8 -*-
"""
@author: gerardo
"""

## Creado por Juan Gerardo Castrejon
## Kalman extendido 
## 06/09/2019
#%%
import numpy as np
from numpy.random import randn
from numpy.linalg import inv
import matplotlib.pyplot as plt

#%%
F = np.array([[1, 1],[0, 1]])
H = np.array([1, 0])
G = np.array([0.5, 1])
g = 1

Q = np.array([[0, 0],[0, 0]])
R = 1

#P0 = np.array([])
x0 = np.array([100, 0])

n = 6
u = -g

x = np.zeros((2,n))
z = np.zeros((n,1))
w = 0.5*randn(n,1)

#%%
def xs(xa): 
    return F@xa + G*u

#%%
x[:,0] = x0
for k in range(1,n):
    x[:,k] = xs(x[:,k-1])
    z[k] = H@x[:,k] + w[k]

#%%
xp = np.zeros((2,n))
xc = np.zeros((2,n))
Pp = np.zeros((2,2,n))
Pc = np.zeros((2,2,n))

xc0 = np.array([95, 1])
Pc0 = np.array([[10, 0], [0, 1]])

#%%
xc[:,0] = xc0
Pc[:,:,0] = Pc0

for k in range(1,n):
    # Prediccion
    xp[:,k] = F@xc[:,k-1] + G*u
    Pp[:,:,k] = F@Pc[:,:,k-1]@F.T + Q
    
    # Correccion
    S = H@Pp[:,:,k]@H.T + R
    K = Pp[:,:,k]@H.T/S #OBSERVE: en python inv no puede aplicarse a escalares
    xc[:,k] = xp[:,k] + K*(z[k] - H@xp[:,k])
    Pc[:,:,k] = Pp[:,:,k] - K*S@K.T #OBSERVE: en python @ no puede aplicarse a escalares


#%%
t = range(0,6)
plt.plot(t, x[0,:], t, z, t, xc[0,:])
plt.axis([1,6,85,101])
plt.grid(True)
plt.legend(['True', 'Measurement','Estimate'])
plt.show()
