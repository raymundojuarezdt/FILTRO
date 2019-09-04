#"""***************************************************************************
#*                        Filtro de Kalman Extendido                          *
#* Versión: 1.0                                                               *
#* Autor: Arturo Enrique Gil García                                           *
#* Fecha: 07/Marzo/2018                                                       *
#* Notas: Hasta el momento solo sirve mostrar el sistema con y sin ruido
#****************************************************************************"""
# Librerías útiles
#import scipy # scientific library
import numpy as np # math library se usa para usar arreglos de matrices
import matplotlib.pyplot as plt # Se usa para graficar
import matplotlib.patches as mpatches #No se bien que hace, pero sin el no hace

#%% Declaro mis variables para el sistema real con/sin ruido
T = 0.01
k = 36
t = np.linspace(0,(k-1)*T,k)

# Creo el vector del sistema X y el vector del sistema con ruido Z con c.i=x0
x0 = 0.99
X = np.zeros((1, k))
Z = np.zeros((1, k))
X[0,0] = x0
Z[0,0] = x0

# Creo mi vector de ruido blanco
W = 0.5*np.random.standard_normal((k,))

#%% Obtengo la respuesta a los sistemas real X y con ruido blanco Z y grafico
for n in range(1, k):
    xs = T+(1+T)*X[0,n-1] + (T/2)*X[0,n-1]**2+(T/6)*X[0,n-1]**3
    X[0,n] = xs
    Z[0,n] = X[0,n] + W[n-1]
    
# Despliego resultados
plt.figure(1)
plt.plot(t,X.transpose())
plt.xlabel('time (s)')
plt.ylabel('function')
plt.grid(True)
plt.title('Real System - Pyton :D')

red_patch = mpatches.Patch(color='blue', label='x(t)')
plt.legend(handles=[red_patch],loc=1)
plt.savefig('SistemaSinRuido.png', format='png',dpi = 700) 
#plt.show()

plt.figure(2)
plt.plot(t,Z.transpose())
plt.grid(True)
plt.title('System + white noise :(')
plt.savefig('SistemaRuidoso.png', format='png',dpi = 700) 

#%% Filtro de Kalman Extendido