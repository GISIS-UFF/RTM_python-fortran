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

# This variable will tell fortran that 
# we want to do the modelling of the seismograms
regTTM = 0

start_time = time.time()

# Velocity Model used 
C = aux.readbinaryfile(parametro.Nz,parametro.Nx,parametro.modeloreal)
aux.plotmodel(C,'jet')

# Creates the seismic source
fortran.wavelet(1,parametro.dt,1,parametro.f_corte) 

# Shows the seismic pulse
aux.plotgraphics(2,'wavelet_ricker.dat', 'k')
#pl.show()

# Define the source's samples
lixo, fonte = np.loadtxt('wavelet_ricker.dat', unpack = True)

# Creates the damping layer function
func_amort = aux.amort(parametro.fat,parametro.nat)
aux.plotgraphics(1,'f_amort.dat','k')
#pl.show()

# Creates a file with the source positions
if parametro.gera_pos_fonte: 
    aux.posicao_fonte(parametro.Nz,\
                      parametro.Nx,\
                      parametro.N_shot,\
                      parametro.Fx0,\
                      parametro.Fz0,\
                      parametro.SpaFonte)

# Loads the source position
Fx, Fz = np.loadtxt('posicoes_fonte.dat',dtype = 'int',unpack = True)

print("Number of shots=",parametro.N_shot)

# Case 1: Only one shot
if parametro.N_shot == 1:
    print("Fx =", Fx, "Fz =", Fz, "shot",parametro.N_shot)
    fortran.nucleomodelagem(parametro.Nz,\
                            parametro.Nx,parametro.Nt,\
                            parametro.h,parametro.dt,\
                            parametro.nat,\
                            parametro.N_shot,\
                            parametro.shotshow,\
                            Fx,\
                            Fz,\
                            fonte,\
                            parametro.Nsnap,\
                            regTTM,\
                            parametro.modeloreal,\
                            parametro.sismogramaobservado,\
                            parametro.nome_prin,\
                            parametro.zr,)
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
            proc = mp.Process(target=aux.modelagemparalela,\
            args=(shot+1,\
            Fx[shot],\
            Fz[shot],\
            fonte,\
            regTTM,\
            parametro.sismogramaobservado,\
            parametro.modeloreal,\
            parametro.nome_prin))

            procs.append(proc)
            proc.start()
        
        for proc in procs:
            proc.join()

elapsed_time_python = time.time() - start_time
print ("Tempo de processamento python = ", elapsed_time_python, "s")
