from os import system as os

#Clear data directory

os("rm $DarkMatterDir*.txt")
os("rm $DarkMatterDir/RT/*.txt")

os("gfortran -fopenmp -Werror -ffree-form -fdefault-real-8 -Ofast parameters.f constants.f IO.f metric.f derivatives.f NumericalMethods.f OrbitalDynamics.f RayTracing.f main.f -J mod/") 

os("./a.out")
