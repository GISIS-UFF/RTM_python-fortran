"""
Read Parameters
"""
from numpy import zeros

filename = '../data/marmousi_vp_383x141.bin'

Nx       = 383               # Numero de pontos no Grid (x)
Nz       = 141               # Numero de pontos no Grid (z)
h        = 10                # Espacamento do Grid
dt       = 1.0e-04           # Incremento de tempo
Pc       = zeros((Nz,Nx)) # Matriz do Tempo Atual (t)
Pf       = zeros((Nz,Nx)) # Matriz do Tempo Futuro (t+1)
C        = zeros((Nz,Nx)) # Matriz do Modelo de Velocidade
f_corte  = 30
