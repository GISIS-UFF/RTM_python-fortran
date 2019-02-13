# Program to read the migration images, sums they and then saves the result in a binary file

# Modules from Python

import time
import matplotlib.pyplot as pl
import numpy as np
import multiprocessing as mp

# Modules created 
import parametro
import auxfunctionsmodule as aux
import fortransubroutines as fortran

start_time = time.time()

StackImage = np.zeros((parametro.Nz,parametro.Nx)) # variable responsible for storing the migrated images

# Loads the source position
Fx, Fz = np.loadtxt('posicoes_fonte.dat',dtype = 'int',unpack = True)
N_shot = np.size(Fx)

# Generate binary with final image
fortran.savefinalimage(parametro.Nz,\
                       parametro.Nx,\
                       parametro.N_shot,\
                       parametro.caminho_migracao,\
                       parametro.nome_prin)

# plot Stack image generated in fortran
filename_imagem = "../Imagem/"'%s'%(parametro.nome_prin)+"_FinalImage.bin"
StackImage  =  aux.readbinaryfile(parametro.Nz,parametro.Nx,filename_imagem)     

aux.plotmodel(StackImage,'gray')
pl.show()

elapsed_time_python = time.time() - start_time
print ("Tempo de processamento python = ", elapsed_time_python, "s")

