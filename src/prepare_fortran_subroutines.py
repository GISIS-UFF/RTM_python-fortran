import os

os.system("mkdir ../sismograma")
os.system("mkdir ../snapshot")
os.system("mkdir ../sismograma_modelo_camada_de_agua")
os.system("mkdir ../sismograma_sem_onda_direta")
os.system("mkdir ../snapshot_migracao_rtm")
os.system("mkdir ../matriz_tempo_transito")
os.system("mkdir ../Imagem")

os.system("rm *.so")
os.system("f2py -c -m  fortransubroutines fortransubroutines.f90")
#os.system("clear")
