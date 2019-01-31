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


Mig_image = aux.readbinaryfile(parametro.Nz,parametro.Nx,"../Imagem/Plano_Paralelo_shot001.bin")
aux.plotmodel(Mig_image,'gray')
pl.show()



