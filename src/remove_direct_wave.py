import parametro
import time
import matplotlib.pyplot as pl
import numpy as np
import multiprocessing as mp

import auxfunctionsmodule as aux
import fortransubroutines as fortran

regTTM = 0

start_time = time.time()

# Define o numero de amostas da fonte
lixo, fonte = np.loadtxt('wavelet_ricker.dat', unpack = True)
Nfonte      = np.size(fonte)

# Carrega posicao da fonte
Fx, Fz = np.loadtxt('posicoes_fonte.dat',dtype = 'int',unpack = True)
N_shot = np.size(Fx)

print("Modelagem com modelo Homogeneo para remover onda direta")
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
            
else: # Se numeros de tiros e maior que 1 use a paralelizacao
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

print("Removendo a onda direta")
if N_shot == 1:
     print("Fx =", Fx, "Fz =", Fz, "shot",N_shot)
     fortran.removeondadireta(parametro.Nt,\
                              parametro.Nx,\
                              N_shot,\
                              parametro.nome_prin)
            
else: # Se numeros de tiros e maior que 1 use a paralelizacao
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
