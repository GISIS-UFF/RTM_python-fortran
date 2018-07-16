## Modules that will be used in this code

# Modules from python
import time
import matplotlib.pyplot as pl
import numpy as np
import multiprocessing as mp

# Modules created 
import parametro
import auxfunctionsmodule as aux
import fortransubroutines as fortran

# This variable will tell fortran tha we want to do the modelling of the seismograms
regTTM = 0 

start_time = time.time()

# Define the source's samples
lixo, fonte = np.loadtxt('wavelet_ricker.dat', unpack = True)
Nfonte      = np.size(fonte)

# Loads the source position
Fx, Fz = np.loadtxt('posicoes_fonte.dat',dtype = 'int',unpack = True)
N_shot = np.size(Fx)

# This part of the script is responsable for the homogeneous seismogram creation
print("Modelagem com modelo Homogeneo para remover onda direta")

# Case 1: Only one shot
if N_shot == 1:
    print("Fx =", Fx, "Fz =", Fz, "shot",N_shot)
    fortran.nucleomodelagem(parametro.Nz,parametro.Nx,parametro.Nt,\
                    parametro.h,parametro.dt,parametro.nat,\
                    N_shot,parametro.shotshow,\
                    Fx,Fz,fonte,parametro.Nsnap,regTTM,\
                    parametro.modelocamadadeagua,parametro.sismogramacamadadeagua,\
                    parametro.nome_prin,\
                    parametro.zr,)

    print("shot=",shot,"Finalizado")
            
# Case 2: More than one shot -> use parallelization
else:
    procs = []    
    for shot in np.arange(0,N_shot):
        proc = mp.Process(target=aux.modelagemparalela, \
        args=(shot+1,\
        Fx[shot],\
        Fz[shot],\
        fonte,\
        regTTM,\
        parametro.sismogramacamadadeagua,\
        parametro.modelocamadadeagua,\
        parametro.nome_prin))
        procs.append(proc)
        proc.start()
    
    for proc in procs:
        proc.join()

# This part of the script is responsable for the direct wave removal
print("Removendo a onda direta")

# Case 1: Only one shot
if N_shot == 1:
     print("Fx =", Fx, "Fz =", Fz, "shot",N_shot)
     fortran.removeondadireta(parametro.Nt,\
                              parametro.Nx,\
                              N_shot,\
                              parametro.nome_prin)

# Case 2: More than one shot -> use parallelization
else: 
     procs = []    
     for shot in np.arange(0,N_shot):
         proc = mp.Process(target=aux.remove_onda_direta,\
         args=(shot+1,\
         Fx[shot],\
         Fz[shot],\
         parametro.nome_prin))
         procs.append(proc)
         proc.start()
    
     for proc in procs:
         proc.join()


elapsed_time_python = time.time() - start_time
print ("Tempo de processamento python = ", elapsed_time_python, "s")
