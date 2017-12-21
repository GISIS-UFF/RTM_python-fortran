
import numpy as np
import matplotlib.pylab as pl
import time

start_time = time.time()



#import op_nsg


# Variaveis 

Nx = 383                   # Numero de pontos no Grid (x)
Nz = 141                   # Numero de pontos no Grid (z)
#Nx = 2                      # Numero de pontos no Grid (x)
#Nz = 3                      # Numero de pontos no Grid (z)
h  = 10                     # Espacamento do Grid
dt = 1.0e-04                # Incremento de tempo
Pc = np.zeros((Nz,Nx))      # Matriz do Tempo Atual (t)
Pf = np.zeros((Nz,Nx))      # Matriz do Tempo Futuro (t+1)
C  = np.zeros((Nz,Nx))      # Matriz do Modelo de Velocidade
f_corte = 30

# Abrindo o Marmousi

filename = '../data/marmousi_vp_383x141.bin'

with open(filename, 'rb') as f:
    
      data = np.fromfile(f, dtype=np.float32, count= Nz*Nx)
      C = np.reshape(data, [Nz, Nx], order='F')

pl.imshow(C,cmap='jet')
pl.colorbar()
pl.show()

#%%
# Matrizes Simples Aleatorias 

# C = np.ones((Nz,Nx))
# Pc = np.random.rand(Nz,Nx)
# Pf = np.random.rand(Nz,Nx)

# Loop do Tempo

#* Colocar depois

# Operador de Diferencas Finitas de Quarta Ordem

#op_nsg.operador_quarta_ordem(h,dt,C,Pc,Pf)

#op_nsg.wavelet(1,dt,1,f_corte)

#W  = 'wavelet_ricker.dat'


