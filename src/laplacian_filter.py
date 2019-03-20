# Modules from python
import time
import matplotlib.pyplot as pl
import numpy as np
import multiprocessing as mp

# Modules created 
import parametro
import auxfunctionsmodule as aux
import fortransubroutines as fortran


start_time = time.time()
# plot Stack image generated in fortran
filename_imagem = parametro.caminho_migracao + parametro.nome_prin + "_FinalImage.bin"
StackImage  =  aux.readbinaryfile(parametro.Nz,parametro.Nx,filename_imagem)     

aux.plotmodel(StackImage,'gray')
#pl.show()

fortran.laplacian(parametro.Nz,\
                    parametro.Nx,\
                    parametro.h,\
                    parametro.caminho_migracao,\
                    parametro.nome_prin+"_FinalImage")

# plot Stack image generated in fortran
filename_imagem = parametro.caminho_migracao + parametro.nome_prin + "_FinalImage_Laplacian.bin"
Laplacian  =  aux.readbinaryfile(parametro.Nz,parametro.Nx,filename_imagem)         

aux.plotmodel(Laplacian,'gray')
pl.show()

elapsed_time_python = time.time() - start_time
print ("Tempo de processamento python = ", elapsed_time_python, "s")