#!/usr/bin/python
#  Plot 2D binary file.          
#  This script import a 2D binary file
#  and plot as image
#
#  INPUT:  
#  filename      = path with the 2D binary file.
#  dim1          = Number of elements in first dimension
#  dim2          = Number of elements in second dimension
# 
#  OUTPUT: 
#  file          = None;
#  
#  Code Written by Felipe Timoteo and Cintia Queiroz
#                  Last update: 10 May, 2019
# 
#  Copyright (C) 2019 Grupo de Imageamento Sísmico e Inversão Sísmica (GISIS)
#                     Departamento de Geologia e Geofísica
#                     Universidade Federal Fluminense
###############################################################################

def readbinaryfile(dim1,dim2,filename):
      """
      readbinaryfile - Functions that read a binary file.
      Usage
      Input:
      dim1     = Number of sample of 1st Dimension
      dim2     = Number of sample of 2nd Dimension
      filename = path of binary file     
      """      
      import numpy as np
      with open(filename, 'rb') as f:    
            data   = np.fromfile(filename, dtype= np.float32, count= dim1*dim2)
            matrix = np.reshape(data, [dim1,dim2], order='F')
      return matrix


def plotmatrix(matrix,colormap):
      """
      Plot a 2D matrix read from binary file.

      """      
      import matplotlib.pyplot as pl
      from mpl_toolkits.axes_grid1 import make_axes_locatable

      pl.figure()
      ax = pl.gca()
      im = pl.imshow(matrix,cmap= colormap)

      # create an axes on the right side of ax. The width of cax will be 5%
      # # of ax and the padding between cax and ax will be fixed at 0.05 inch.
      divider = make_axes_locatable(ax)
      cax = divider.append_axes("right", size="5%", pad=0.05)

      pl.colorbar(im, cax=cax)
      #pl.draw() # drawing figure to be plotted later 
      pl.show()


if __name__ =="__main__":
      '''
      Ambiente para teste das funções criadas      
      '''

      import parametro

      
      filename  = parametro.modelosuavizado
      #filename  = '../sismograma/Marmousi_v2_sismograma002.bin'
      dim1=parametro.Nz
      dim2=parametro.Nx

      
      print('Plotting: '+ filename)
      print('Vertical dimension   = ', dim1)
      print('Horizontal dimension = ', dim2)
      matrix = readbinaryfile(dim1,dim2,filename)
      plotmatrix(matrix,'jet')
      