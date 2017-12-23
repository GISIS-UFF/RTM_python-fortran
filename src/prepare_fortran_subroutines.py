import os

os.system("rm *.so")
os.system("f2py -c -m  fortransubroutines fortransubroutines.f90")
os.system("clear")
