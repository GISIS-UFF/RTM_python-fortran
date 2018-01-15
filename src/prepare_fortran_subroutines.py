import os

os.system("mkdir ../sismograma")
os.system("mkdir ../snapshot")
os.system("mkdir ../modeloreal")
os.system("rm *.so")
os.system("f2py -c -m  fortransubroutines fortransubroutines.f90")
#os.system("clear")
