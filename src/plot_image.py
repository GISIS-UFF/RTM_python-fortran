
# Modules from Python

import time
import matplotlib.pyplot as pl
import numpy as np
import multiprocessing as mp

# Modules created 
import parametro
import auxfunctionsmodule as aux
import fortransubroutines as fortran

filename_imagem = "../imagens_artigo/MarmousiRTMTT_Seismogram020.bin"
Image  =  aux.readbinaryfile(parametro.Nt,parametro.Nx,filename_imagem)  
#Image[0:36,:] = 0.0
#Image = Image/np.max(np.abs(Image))
aux.plotmodel(Image,'gray')
pl.show()