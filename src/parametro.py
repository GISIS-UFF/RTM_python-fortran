# -*- coding: utf-8 -*-
"""
Read Parameters
"""
from numpy import zeros

modeloreal = '../modelos_utilizados/marmousi_vp_383x141.bin'#'../modelos_utilizados/modelo_v2_VR_275x701.bin'
modelosuavizado = '../modelos_utilizados/Suave_v15_marmousi_vp_383x141.bin'
modelohomogeneo = '../modelo_utilizados/velocitymodel_Homo_383x141.bin'



Nx         = 383  # 701             # Numero de pontos no Grid (x)
Nz         = 141  #275              # Numero de pontos no Grid (z)
h          = 10   #0.15             # Espacamento do Grid
dt         = 1.0e-04 #5.0e-04       # Incremento de tempo
gera_pos_fonte =1           	    # Gera ou não um arquivo automatico com as posicoes da fonte(0 nao gera, 1 gera)
N_shot     = 10                      # Numero total de tiros
Fx0        = 10 #6                  # Posicao inicial da fonte em X
Fz0        = 2                      # Posicao inicial da fonte em Z
SpaFonte   = 10                     # Espacamento entre as posicoes da fonte
f_corte    = 30 #350                # Frequencia de corte
fat	   = 1.5e-03	            # Fator de amortecimento
nat	   = 80 	            # Numero de pontos do Grid que farao parte da camada de amortecimento
T   	   = 3 #0.10		    # Tempo final
Nt	   = int(T/dt)              # Numero de amostras
beta	   = 4		            # Parametro de estabilidade
alfa	   = 5		            # Parametro de dispersao
shotshow   = 1                      # Tiro para snapshot (0 nao gera snapshot)
Nsnap      = 20                     # Numero de Snapshots 
zr         = 2                      # Profundidade dos receptores
