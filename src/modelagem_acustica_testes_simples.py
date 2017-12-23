# -*- coding: utf-8 -*-
"""
Created on Sat Dec 23 14:06:01 2017

@author: Cintiaq
"""

import numpy as np
#import matplotlib.pylab as pl
import time
import fortransubroutines

start_time = time.time()

Nx = 2                      # Numero de pontos no Grid (x)
Nz = 3                      # Numero de pontos no Grid (z)
h  = 10                     # Espacamento do Grid
dt = 1.0e-04                # Incremento de tempo
Pc = np.zeros((Nz,Nx))      # Matriz do Tempo Atual (t)
Pf = np.zeros((Nz,Nx))      # Matriz do Tempo Futuro (t+1)
C  = np.zeros((Nz,Nx))      # Matriz do Modelo de Velocidade
f_corte = 30                # Frequencia de corte

#%%Matrizes Simples Aleatorias 

C = np.ones((Nz,Nx))
Pc = np.random.rand(Nz,Nx)
Pf = np.random.rand(Nz,Nx)

# Loop do Tempo

#* Colocar depois

#Operador de Diferencas Finitas de Quarta Ordem

fortransubroutines.operador_quarta_ordem(h,dt,C,Pc,Pf)