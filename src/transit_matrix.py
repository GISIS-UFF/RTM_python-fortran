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

# This variable will tell fortran tha we want to do the modelling of the transit time matrix
regTTM = 1 

start_time = time.time()

# Velocity Model used 

C = aux.readbinaryfile(parametro.Nz,parametro.Nx,parametro.modelosuavizado)
aux.plotmodel(C, 'jet')

# Create the seismic source

fortran.wavelet(1,parametro.dt,1,parametro.f_corte)

# Shows the seismic pulse

aux.plotgraphics(2,'wavelet_ricker.dat','k')
#pl.show()

# Define the source's samples

lixo, fonte = np.loadtxt('wavelet_ricker.dat', unpack = True)
Nfonte      = np.size(fonte)

# Create damping layer

func_amort = aux.amort(parametro.fat,parametro.nat)
aux.plotgraphics(1,'f_amort.dat','k')
#pl.show()

if parametro.gera_pos_fonte:
    aux.posicao_fonte(parametro.Nz,parametro.Nx,parametro.N_shot, parametro.Fx0,parametro.Fz0,parametro.SpaFonte)

# Loads the source position

Fx, Fz = np.loadtxt('posicoes_fonte.dat',dtype = 'int',unpack = True)
N_shot = np.size(Fx)

print(N_shot)

# Case 1: Only one shot

if N_shot == 1:
     print("Fx =", Fx, "Fz =", Fz, "shot",N_shot)
     fortran.nucleomodelagem(parametro.Nz,parametro.Nx,parametro.Nt,\
                                parametro.h,parametro.dt,parametro.nat,\
                                N_shot,parametro.shotshow,\
                                Fx,Fz,fonte,parametro.Nsnap,regTTM,\
                                parametro.modelosuavizado,parametro.caminho_TTM,\
                                parametro.nome_prin,\
                                parametro.zr,)
     print(" shot= ",shot," Finalizado.")

# Case 2: More than one shot -> use parallelization

else:
     procs = []    
     for shot in np.arange(0,N_shot):
        proc = mp.Process(target=aux.modelagemparalela, args=(shot+1,Fx[shot],Fz[shot],fonte,regTTM,\
                          parametro.caminho_TTM,parametro.modelosuavizado,parametro.nome_prin))
        procs.append(proc)
        proc.start()
    
     for proc in procs:
        proc.join()

TTM = aux.readbinaryfile(parametro.Nz,parametro.Nz,'../matriz_tempo_transito/Marmousi_TTM_shot001.bin')
aux.plotmodel(TTM,'jet')
pl.show()

elapsed_time_python = time.time() - start_time

print ("Tempo de processamento python = ", elapsed_time_python, "s")
