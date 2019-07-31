module constants

use parameters

implicit none

real(kind=dp), parameter :: Newton_g = 6.67408d-11 
real(kind=dp), parameter :: Msolar = 1.989d30 
real(kind=dp), parameter :: mu = Newton_g*MBH*Msolar
real(kind=dp), parameter :: light_c = 3.0d8
real(kind=dp), parameter :: convert_m = light_c**2/(Newton_g*MBH*Msolar) !Multiply by this to go TO Natural units
real(kind=dp), parameter :: convert_s = light_c**3/(Newton_g*MBH*Msolar) !Multiply by this to go TO Natural units
real(kind=dp), parameter :: convert_spin= light_c/(Newton_g*(MBH*Msolar)**2.0_dp) !Multiply by this to go TO Natural units
real(kind=dp), parameter :: RPSR_M = RPSR * 1.0d3 * convert_m
real(kind=dp), parameter :: m = sqrt(rCOM**2 + a**2)




!Global constants declared later
real(kind=dp) :: E, Lz, kappa

end module constants
