import os

os.system("rm *.so")
os.system("f2py -c -m  fortran_subroutines fortran_subroutines.f90")
os.system("clear")
