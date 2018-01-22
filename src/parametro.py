"""
Read Parameters
"""
from numpy import zeros

modeloreal = '../modelo_real/marmousi_vp_383x141.bin'

Nx         = 383               # Numero de pontos no Grid (x)
Nz         = 141               # Numero de pontos no Grid (z)
h          = 10                # Espacamento do Grid
dt         = 1.0e-04           # Incremento de tempo
shot       = 1                 # Numero do tiro 
f_corte    = 30                # Frequencia de corte
fat	   = 1.5e-03	       # Fator de amortecimento
nat	   = 80 	       # Numero de pontos do Grid que farao parte da camada de amortecimento
T   	   = 4		       # Tempo final
Nt	   = int(T/dt)         # Numero de amostras
beta	   = 4		       # Parametro de estabilidade
alfa	   = 5		       # Parametro de dispersao
shotshow   = 1                 # Tiro para snapshot (0 nao gera snapshot)
Nsnap      = 21                # Numero de Snapshots 
