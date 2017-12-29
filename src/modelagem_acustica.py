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
      pl.imshow(matrix,cmap='jet')
      pl.colorbar()
      pl.draw() # drawing figure to be plotted later
      
      return matrix
  
    
def plotgraphics(filename,color):
    
    """
    
    plotgraphics - Function that read a file and plot a graphic.
    
    REMEMBER: The function 'np.loadtxt' only can be used with files written
              with dot 
    
    """
    
    X,Y = np.loadtxt(filename, unpack = True)

    pl.figure()
    pl.plot(X,Y, c = color)
    pl.draw()
    
    return 
    
    
        
      
def main():      
      '''
      Main program is here      
      '''
      
      import parametro      
      #from fortransubroutines import wavelet

      C = readbinaryfile(parametro.Nz,parametro.Nx,parametro.filename)

      #wavelet(1,parametro.dt,1,parametro.f_corte)
      
      fonte = plotgraphics('wavelet_ricker.dat', 'k')

      # Loop do Tempo
      
      #* Colocar depois
      
      # Operador de Diferencas Finitas de Quarta Ordem      
      #op_nsg.operador_quarta_ordem(h,dt,C,Pc,Pf)


      
if __name__ == '__main__':
    
      """
      Structure prepared for object-oriented programming
      """
      
      import numpy as np      
      import matplotlib.pylab as pl      
      import time
      
      start_time = time.time()


      main() # Call main function


      elapsed_time_python = time.time() - start_time
      print ("Tempo de processamento python = ", elapsed_time_python,"s")


      pl.show() # Showing all figures draw

