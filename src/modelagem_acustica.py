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

 
def estabilidade(C,h,beta,dt):
    
    if dt > h / beta * np.max(np.max(C)):

       raise ValueError ("Erro de Estabilidade")

def dispersao(C,h,alfa,f_corte):
    
    if h > np.min(np.min(C)) / (alfa * f_corte):

        raise ValueError ("Erro de Dispersao Numerica")   

   
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
	
def plotmodel(matrix):

      pl.figure()
      pl.imshow(matrix,cmap='jet')
      pl.colorbar()
      pl.draw() # drawing figure to be plotted later
      
def plotseism(Sismograma,Nt,Nx):

      pl.figure()
      pl.imshow(Sismograma,extent=[0,Nt,1,Nx],aspect=100)    
      pl.imshow(Sismograma,aspect="auto",cmap = pl.cm.gray)    
      pl.draw()


def plotsnaps(filesnap):
    
    with open(filesnap, 'rb') as f:
        pl.figure()
        pl.imshow(f, 'jet')
        pl.draw()
        
	
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
      from fortransubroutines import wavelet
      from fortransubroutines import modelagem
      
      
      # Modelo de Velocidade

      C = readbinaryfile(parametro.Nz,parametro.Nx,parametro.modeloreal)
      plotmodel(C)    
      
      # Condicao de estabilidade e dispersao
      
      estabilidade(C,parametro.h,parametro.beta,parametro.dt)
      dispersao(C,parametro.h,parametro.alfa,parametro.f_corte)
      
      # Fonte Sismica
      
      wavelet(1,parametro.dt,1,parametro.f_corte)
      
      plotgraphics(2,'wavelet_ricker.dat', 'k')
      
      lixo, fonte = np.loadtxt('wavelet_ricker.dat', unpack = True)
      Nfonte      = np.size(fonte)

      # Funcao Amortecedora
      
      amort(parametro.fat,parametro.nat)
      
      plotgraphics(1,'f_amort.dat', 'k')
      
      # Modelagem
      
      Fx = int(parametro.Nx/2)               # Posicao da Fonte (x)
      Fz = int(parametro.Nz/2)               # Posicao da Fonte (z)
      shot = 1                               # Numero do tiro 

      modelagem(parametro.Nz,parametro.Nx,parametro.Nt,parametro.h,parametro.dt,\
                shot,Fx,Fz,fonte,Nfonte,parametro.snapshot,parametro.Nsnap)
               
      Sismograma = readbinaryfile(parametro.Nt,parametro.Nx,"../sismograma/Marmousi_sismograma001.bin")
 
      plotseism(Sismograma,parametro.Nt,parametro.Nx)


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

