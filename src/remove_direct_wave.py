import parametro
import time
import matplotlib.pyplot as pl
import numpy as np
import multiprocessing as mp

import auxfunctionsmodule as aux
import fortransubroutines as fortran

regTTM = 0

start_time = time.time()

# Modelo de Velocidade Usado

C = aux.readbinaryfile(parametro.Nz,parametro.Nx,parametro.modelocamadadeagua)
aux.plotmodel(C,'jet')

# Gera fonte sismica - não precisa, gerado com os dados observados
#fortran.wavelet(1,parametro.dt,1,parametro.f_corte) 

# Visualiza pulso sismico
#aux.plotgraphics(2,'wavelet_ricker.dat', 'k')
#pl.show()

# Define o numero de amostas da fonte
lixo, fonte = np.loadtxt('wavelet_ricker.dat', unpack = True)
Nfonte      = np.size(fonte)

# Cria camada de amortecimento
func_amort = aux.amort(parametro.fat,parametro.nat)
aux.plotgraphics(1,'f_amort.dat','k')
#pl.show()

# Carrega posicao da fonte
Fx, Fz = np.loadtxt('posicoes_fonte.dat',dtype = 'int',unpack = True)
N_shot = np.size(Fx)

# Modelagem com modelo homogeneo - usado para remover a onda direta
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

# Removendo a onda direta

# if N_shot == 1:
#     print("Fx =", Fx, "Fz =", Fz, "shot",N_shot)
#     fortran.removeondadireta(parametro.Nt,parametro.Nx,N_shot)
            
# else: # Se numeros de tiros e maior que 1 use a paralelizacao
#     procs = []    
#     for shot in np.arange(0,N_shot):
#         proc = mp.Process(target=aux.remove_onda_direta, args=(shot+1,Fx[shot],Fz[shot]))
#         procs.append(proc)
#         proc.start()
    
#     for proc in procs:
#         proc.join()


elapsed_time_python = time.time() - start_time
print ("Tempo de processamento python = ", elapsed_time_python, "s")
