<<<<<<< HEAD
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
      
def main():      
      '''
      Main program is here      
      '''
      
      import parametro      
      from fortransubroutines import wavelet

      C = readbinaryfile(parametro.Nz,parametro.Nx,parametro.filename)

      wavelet(1,parametro.dt,1,parametro.f_corte)

      # Matrizes Simples Aleatorias 

      # C = np.ones((Nz,Nx))
      # Pc = np.random.rand(Nz,Nx)
      # Pf = np.random.rand(Nz,Nx)

      # Loop do Tempo
      
      #* Colocar depois
      # Operador de Diferencas Finitas de Quarta Ordem      
      #op_nsg.operador_quarta_ordem(h,dt,C,Pc,Pf)



<<<<<<< HEAD
#W  = 'wavelet_ricker.dat'

      
if __name__ == '__main__':
      """
      Structure prepared for object-oriented programming
      """
      
      import numpy as np      
      import matplotlib.pylab as pl      
      import time
      
      start_time = time.time()
=======
=======
import numpy as np
import matplotlib.pylab as pl
import time
#import fortransubroutines

start_time = time.time()

# Variaveis 

>>>>>>> Plot da Wavelet adicionado
Nx = 383                    # Numero de pontos no Grid (x)
Nz = 141                    # Numero de pontos no Grid (z)
h  = 10                     # Espacamento do Grid
dt = 1.0e-04                # Incremento de tempo
Pc = np.zeros((Nz,Nx))      # Matriz do Tempo Atual (t)
Pf = np.zeros((Nz,Nx))      # Matriz do Tempo Futuro (t+1)
C  = np.zeros((Nz,Nx))      # Matriz do Modelo de Velocidade
f_corte = 30

# Abrindo o Marmousi

filename = '../data/marmousi_vp_383x141.bin'

with open(filename, 'rb') as f:
    
      data = np.fromfile(f, dtype=np.float32, count= Nz*Nx)
      C = np.reshape(data, [Nz, Nx], order='F')

pl.figure(1)	  
pl.imshow(C,cmap='jet')
pl.colorbar()
pl.show()

#%% Loop do Tempo

#* Colocar depois

# Operador de Diferencas Finitas de Quarta Ordem

#fortransubroutines.operador_quarta_ordem(h,dt,C,Pc,Pf)
>>>>>>> Tentativa de abrir o arquivo wavelet_ricker.dat adicionada; Não deu para testar devido a problemas da f2py no Windows

<<<<<<< HEAD
      main() # Call main function

<<<<<<< HEAD
      elapsed_time_python = time.time() - start_time
      print ("Tempo de processamento python = ", elapsed_time_python,"s")
=======
# Abrindo o arquivo wavelet.dat criado 
=======
#fortransubroutines.wavelet(1,dt,1,f_corte)

# Lendo o arquivo wavelet.dat criado e armazenando em dois vetores: X,Y
>>>>>>> Plot da Wavelet adicionado

X,Y = np.loadtxt('wavelet_ricker.dat', unpack = True)

pl.figure(2)
pl.plot(X,Y, color = 'k')
pl.show()


>>>>>>> Tentativa de abrir o arquivo wavelet_ricker.dat adicionada; Não deu para testar devido a problemas da f2py no Windows

      pl.show() # Showing all figures draw

