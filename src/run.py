from os import system as os




#Clear data directory

#os("rm $RayTracingDir/*.txt")





os("gfortran -fopenmp -ffree-form -fdefault-real-8 -O3 *.f -J mod/") 
os("gfortran -fopenmp -ffree-form -fdefault-real-8 -O3 *.f -J mod/") 

os("./a.out") 
