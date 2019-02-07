# -*- coding: utf-8 -*-

# This script is a module with the parameters used in the modelling and migration scripts
"""
Read Parameters
"""
from numpy import zeros

# Variables with the folders that will be used by Fortran
modeloreal              = '../modelos_utilizados/marmousi_vp_383x141.bin'
modelosuavizado         = '../modelos_utilizados/Suave_v15_marmousi_vp_383x141.bin'
modelocamadadeagua      = '../modelos_utilizados/velocitymodel_Hmgns_wtrly.bin'

sismogramaobservado     = '../sismograma/'
sismogramacamadadeagua  = '../sismograma_modelo_camada_de_agua/'
sismogramasemondadireta = '../sismograma_sem_onda_direta/'

caminho_TTM = '../matriz_tempo_transito/'
caminho_migracao = '../Imagem/' 

# Parameters 
nome_prin               = 'Marmousi'

Nx         = 383 #214    # Numero de pontos no Grid (x)
Nz         = 141 #121    # Numero de pontos no Grid (z)
h          = 10          # Espacamento do Grid
dt         = 1.0e-04     # Incremento de tempo
gera_pos_fonte = True  	 # Criar arquivo com posições da fonte? (False or True)
N_shot     = 18          # Numero total de tiros
processors = 3           # Numero de processadores usados simultaneamente
Fx0        = 10          # Posicao inicial da fonte em X
Fz0        = 2           # Posicao inicial da fonte em Z
SpaFonte   = 20          # Espacamento entre as posicoes da fonte
f_corte    = 30          # Frequencia de corte
fat	       = 1.5e-03	 # Fator de amortecimento
nat	       = 80 	     # Numero de pontos do Grid que farao parte da camada de amortecimento
T   	   = 2       	 # Tempo final
Nt	       = int(T/dt)   # Numero de amostras
beta	   = 4	         # Parametro de estabilidade
alfa	   = 5	#3       # Parametro de dispersao
shotshow   = 0           # Tiro para snapshot (0 nao gera snapshot)
Nsnap      = 20          # Numero de Snapshots 
zr         = 2           # Profundidade dos receptores