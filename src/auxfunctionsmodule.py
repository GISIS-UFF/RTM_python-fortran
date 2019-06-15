# Modules from python

import numpy as np
import matplotlib.pylab as pl
from mpl_toolkits.axes_grid1 import make_axes_locatable

# Modules created 

import parametro
import fortransubroutines as fortran

def readbinaryfile(dim1,dim2,filename):
      """
      readbinaryfile - Functions that read a binary file.
      Usage
      Input:
      dim1     = Number of sample of 1st Dimension
      dim2     = Number of sample of 2nd Dimension
      filename = path of binary file     
      """      
      with open(filename, 'rb') as f:    
            data   = np.fromfile(filename, dtype= np.float32, count= dim1*dim2)
            data = np.reshape(data, [dim1,dim2], order='F')
      return data

def plotmodel(matrix,colormap):
      """
      Plot a 2D velocity model read from binary file.

      """      

      pl.figure()
      ax = pl.gca()
      im = pl.imshow(matrix,cmap= colormap)

    # create an axes on the right side of ax. The width of cax will be 5%
    # # of ax and the padding between cax and ax will be fixed at 0.05 inch.
      divider = make_axes_locatable(ax)
      cax = divider.append_axes("right", size="5%", pad=0.05)

      pl.colorbar(im, cax=cax)
      pl.draw() # drawing figure to be plotted later 
      #pl.figure()
      #pl.imshow(matrix,cmap= colormap)

      #pl.colorbar()
      #pl.draw() # drawing figure to be plotted later   

def plotgraphics(ID,filename,color):  

      """    
      Plotgraphics -> Function that read a file and plot a graphic.

      REMEMBER: The function 'np.loadtxt' only can be used with files written
            with dot separation

            ID = number of columns of the file    
      """      
      if ID == 2:
            
            X,Y = np.loadtxt(filename, unpack = True)

            pl.figure()
            pl.plot(X,Y, c = color)
            pl.draw()

      if ID == 1:
            
            X = np.loadtxt(filename, unpack = True)

            pl.figure()
            pl.plot(X, c = color)
            pl.draw()
            
      return

def amort(fat_amort,n_grid):
	"""
	amort -> Creates a damping function that will be used in the Cerjan condition.
               The function created is saved in a txt file. 
      
	"""	
	w = np.zeros(n_grid)
  
	for i in np.arange(0,(n_grid)):
		w[i] = np.exp(-(fat_amort * (n_grid-i)) ** 2)  
	  
	np.savetxt('f_amort.dat',w, delimiter='.')

def posicao_fonte(Nz,Nx,N_shot,Fx0,Fz0,SpaFonte):

      """
      posicao_fonte -> Creates a txt file with the position of the shots.               
      """
      
      posicao = np.zeros((N_shot,2))

      Fz = np.arange(0,N_shot)*0 + Fz0
      Fx = np.arange(0,N_shot)*SpaFonte + Fx0


      # if (N_shot*SpaFonte + Fx0) > (Nx -1):
      #    raise ValueError ("Fonte fora do modelo valido (horizontal). Modifique ou o espacamento da fonte ou a posicao inicial Fx0")

      # if (N_shot*SpaFonte + Fz0) > (Nz -1):
      #    raise ValueError ("Fonte fora do modelo valido (vertical). Modifique a posicao inicial Fz0")
     
      posicao[0:N_shot,0] = Fx
      posicao[0:N_shot,1] = Fz

      np.savetxt("posicoes_fonte.dat",posicao,fmt = '%i')

def modelagemparalela(shot,\
                      Fx,\
                      Fz,\
                      fonte,\
                      regTTM,\
                      sismograma,\
                      modelo,\
                      nome_prin): 

      """
      This function is responsible for the parallelization of the shots in the modelling script.
      """

      print("Fx =", Fx, "Fz =", Fz, "shot", shot)
      fortran.nucleomodelagem(parametro.Nz,parametro.Nx,parametro.Nt,\
                                    parametro.h,parametro.dt,parametro.nat,\
                                    shot,parametro.shotshow,\
                                    Fx,Fz,fonte,parametro.Nsnap,regTTM,\
                                    modelo,sismograma,\
                                    nome_prin,\
                                    parametro.zr,)
      print(" shot= ",shot," Finalizado.")

def migracao_rtm(shot,\
                 Fx,\
                 Fz,\
                 modelo,\
                 nome_prin):

      """
      This function is responsible for the parallelization of the shots in the migration script.
      """

      print("Fx =", Fx, "Fz =", Fz, "shot", shot)

      fortran.migracao(parametro.Nz,parametro.Nx,parametro.Nt,\
                                    parametro.h,parametro.dt,parametro.nat,parametro.zr,\
                                    shot,parametro.shotshow,\
                                    parametro.Nsnap,\
                                    modelo,nome_prin,)
      print(" shot= ",shot," Finalizado.")

def migracao_crosscorrelation(\
            shot,
            Fx,\
            Fz,\
            fonte,\
            modelo,\
            nome_prin):

      """
      This function is responsible for the parallelization of the shots in the migration script.
      """
      print("Fx =", Fx, "Fz =", Fz, "shot",shot)
      fortran.migracaocrosscorrelation(\
                        parametro.Nz,\
                        parametro.Nx,\
                        parametro.Nt,\
                        parametro.h,\
                        parametro.dt,\
                        parametro.nat,\
                        parametro.zr,\
                        shot,\
                        parametro.shotshow,\
                        Fx,\
                        Fz,\
                        fonte,\
                        parametro.Nsnap,\
                        parametro.modelosuavizado,\
                        parametro.nome_prin,)
                        
      print(" shot= ",parametro.N_shot," Finalizado.")


def remove_onda_direta(shot,Fx,Fz,nome_prin):

      """
      This function is responsible for the parallelization of the direct wave's removal.
      """
      print("Fx =", Fx, "Fz =", Fz, "shot", shot)
      fortran.removeondadireta(parametro.Nt,parametro.Nx,shot,nome_prin)
      print(" shot= ",shot," Finalizado.")
