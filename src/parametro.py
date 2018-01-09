"""
Read Parameters
"""
from numpy import zeros

filename = '../data/marmousi_vp_383x141.bin'

Nx       = 383               # Numero de pontos no Grid (x)
Nz       = 141               # Numero de pontos no Grid (z)
h        = 10                # Espacamento do Grid
dt       = 1.0e-04           # Incremento de tempo
Pc       = zeros((Nz,Nx)) 	  # Matriz do Tempo Atual (t)
Pf       = zeros((Nz,Nx))    # Matriz do Tempo Futuro (t+1)
C        = zeros((Nz,Nx))    # Matriz do Modelo de Velocidade
f_corte  = 30                # Frequencia de corte
fat = 1.0e-03                # Fator de amortecimento
nat = 80                     # Numero de pontos do Grid que farao parte da camada de amortecimento
Fx = int(Nx/2)               # Posicao da Fonte (x)
Fz = 2;                      # Pos5icao da Fonte (z)
T = 4                        # Tempo final
nt = int(T/dt)
beta = 4                     # Parametro de estabilidade
alfa = 5                     # Parametro de dispersao