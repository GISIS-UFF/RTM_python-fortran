import numpy as np
import matplotlib.pylab as pl

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
      Output:      
      """      
      with open(filename, 'rb') as f:    
            data   = np.fromfile(f, dtype=np.float32, count= dim1*dim2)
            matrix = np.reshape(data, [dim1,dim2], order='F')
      return matrix

def plotmodel(matrix,colormap):
      """
      Plot 2D velocity model read from binary file.
      """      
      pl.figure()
      pl.imshow(matrix,cmap= colormap)

      pl.colorbar()
      pl.draw() # drawing figure to be plotted later   

def plotgraphics(ID,filename,color):        
      """    
      plotgraphics - Function that read a file and plot a graphic.

      REMEMBER: The function 'np.loadtxt' only can be used with files written
            with dot               
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
	amort - Cria a funcao de amortecimento que vai ser usada nas bordas do modelo e a salva em um arquivo de texto
	
	"""	
	w = np.zeros(n_grid)
  
	for i in np.arange(0,(n_grid)):
		w[i] = np.exp(-(fat_amort * (n_grid-i)) ** 2)  
	  
	np.savetxt('f_amort.dat',w, delimiter='.')

def posicao_fonte(Nz,Nx,N_shot,Fx0,Fz0,SpaFonte):
      
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
                      sismogramaobservado,\
                      modeloreal,\
                      nome_prin):      
      print("Fx =", Fx, "Fz =", Fz, "shot", shot)
      fortran.nucleomodelagem(parametro.Nz,parametro.Nx,parametro.Nt,\
                                    parametro.h,parametro.dt,parametro.nat,\
                                    shot,parametro.shotshow,\
                                    Fx,Fz,fonte,parametro.Nsnap,regTTM,\
                                    modeloreal,sismogramaobservado,\
                                    nome_prin,\
                                    parametro.zr,)
      print(" shot= ",shot," Finalizado.")



def remove_onda_direta(shot,Fx,Fz):
      print("Fx =", Fx, "Fz =", Fz, "shot", shot)
      fortran.removeondadireta(parametro.Nt,parametro.Nx)
      print(" shot= ",shot," Finalizado.")

def square(x,numbers):
 
    for x in numbers:
        print('%s squared  is  %s' % (x, x**2))

def doubler(number):
    """
    A doubling function that can be used by a process
    """
    import os
    result = number * 2
    proc = os.getpid()
    print('{0} doubled to {1} by process id: {2}'.format(
        number, result, proc))        
