#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 24 20:29:37 2019

@author: cintia
"""

import math 
import numpy as np
import matplotlib.pyplot as pl
import parametro
import auxfunctionsmodule as aux

def function_taper(input,x,wx,n_taper):
    taper = np.zeros(n_taper)
    for i in np.arange(0,n_taper):
        taper[i] = i*(math.pi/n_taper)
        taper[i] = math.cos(taper[i])/2 + 0.5
             
    
    [m, n] = input.shape
    mute = np.zeros(n)
           
    mute[x - wx/2:x+wx/2+1] = 1
    mute[x-wx/2-1:x-wx/2-1 - n_taper:-1] = taper
    mute[x+wx/2+1:x+wx/2+1 + n_taper] = taper
        
    test = input * mute    
    return mute


#filename = "../Imagem/Marmousi_v2_shot050.bin"
#N_shot = 50
#matrix = aux.readbinaryfile(parametro.Nz,parametro.Nx) 

#teste = function_taper(matrix,N_shot, 100, 100)

#pl.imshow(matrix)
#pl.imshow(teste)

