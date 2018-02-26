# -*- coding: utf-8 -*-
"""
Read Parameters
"""
from numpy import zeros

modeloreal = '../modelo_real/marmousi_vp_383x141.bin'
modelosuavizado = '../modelo_suavizado/Suave_v15_marmousi_vp_383x141.bin'
modelohomogeneo = '../modelo_homogeneo/velocitymodel_Homo_383x141.bin'



Nx         = 383               # Numero de pontos no Grid (x)
Nz         = 141               # Numero de pontos no Grid (z)
h          = 10                # Espacamento do Grid
dt         = 1.0e-04           # Incremento de tempo
gera_pos_fonte = 1             # 
N_shot     = 2                 # Numero total de tiros
Fx0        = 50
Fz0        = 2
SpaFonte   = 100
f_corte    = 30                # Frequencia de corte
fat	   = 1.5e-03	       # Fator de amortecimento
nat	   = 80 	       # Numero de pontos do Grid que farao parte da camada de amortecimento
T   	   = 1		       # Tempo final
Nt	   = int(T/dt)         # Numero de amostras
beta	   = 4		       # Parametro de estabilidade
alfa	   = 5		       # Parametro de dispersao
shotshow   = 0                 # Tiro para snapshot (0 nao gera snapshot)
Nsnap      = 21                # Numero de Snapshots 
zr         = 3                 # Profundidade dos receptores
