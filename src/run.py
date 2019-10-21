from os import system as os




#Clear data directory

os("rm $DarkMatterDir/*.txt")

os("gfortran -fopenmp -ffree-form -fdefault-real-8 -O3 parameters.f constants.f tensors.f initial_conditions.f integration.f performance.f main.f -J mod/") 

os("./a.out") 
