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
  

def plotgraphics(ID,filename,color):
        
        """
    
            plotgraphics - Function that read a file and plot a graphic.
    
            REMEMBER: The function 'np.loadtxt' only can be used with files written
              with dot 
    
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
	amort - Cria a funcao de amortecimento que vai ser usada nas bordas do modelo e a salva em umarquivo de texto
	
	"""
	
	w = np.zeros(n_grid-1)
  
	for i in np.arange(0,(n_grid -1)):
		w[i] = np.exp(-(fat_amort * (n_grid-i)) ** 2)
  
	  
	np.savetxt('f_amort.dat',w, delimiter='.')  
 

        
def main():      
      '''
      Main program is here      
      '''
      
      import parametro      
      #from fortransubroutines import wavelet
      #from fortransubroutines import modelagem

      C = readbinaryfile(parametro.Nz,parametro.Nx,parametro.filename)

      #wavelet(1,parametro.dt,1,parametro.f_corte)
      
      plotgraphics(2,'wavelet_ricker.dat', 'k')
      
      amort(parametro.fat,parametro.nat)
      
      plotgraphics(1,'f_amort.dat', 'k')
      
      ## Subrotina de Modelagem
      
      snapshots = raw_input("Deseja salvar snapshots(y or n): ")
      
      if snapshots == 'y':
          print 'deu bom!'
          #from fortransubroutines import snapshots
          
      #modelagem(C, ....)
      
      
      
      
      
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
      print ("Tempo de processamento python = ", elapsed_time_python, "s")


      pl.show() # Showing all figures draw

