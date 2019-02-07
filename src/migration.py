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

start_time = time.time()

# Loads the source position
Fx, Fz = np.loadtxt('posicoes_fonte.dat',dtype = 'int',unpack = True)

print("Number of shots=",parametro.N_shot)

# Case 1: Only one shot
if parametro.N_shot == 1:
    print("Fx =", Fx, "Fz =", Fz, "shot",parametro.N_shot)
    fortran.migracao(parametro.Nz,\
                     parametro.Nx,\
                     parametro.Nt,\
                     parametro.h,\
                     parametro.dt,\
                     parametro.nat,\
                     parametro.zr,\
                     parametro.N_shot,\
                     parametro.shotshow,\
                     parametro.Nsnap,\
                     parametro.modelosuavizado,\
                     parametro.nome_prin,)

    print(" shot= ",parametro.N_shot," Finalizado.")

# Case 2: More than one shot -> use parallelization   
else: 
      # list with shot indices
    shots= np.arange(0,parametro.N_shot)
    num_processor=parametro.processors  
    # separete the list of shots according to the number
    # of processors    
    for shot_index in np.arange(0,parametro.N_shot,num_processor):
        
        # Run multiprocessing modeling    
        procs = []    
        for shot in shots[shot_index:shot_index+num_processor]:   
            proc = mp.Process(target=aux.migracao_rtm, \
            args=(shot+1,\
            Fx[shot],\
            Fz[shot],\
            parametro.modelosuavizado,\
            parametro.nome_prin))

            procs.append(proc)
            proc.start()
        
        for proc in procs:
            proc.join()

elapsed_time_python = time.time() - start_time
print ("Tempo de processamento python = ", elapsed_time_python, "s")
