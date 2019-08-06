module parameters
implicit none



!Define float precision
integer, parameter :: dp = selected_real_kind(33,4931)


real(kind=dp), parameter :: PI = 4.D0*ATAN(1.D0) 
real(kind=dp), parameter :: MBH = 4.310d6 !BH mass in solar  !BH mass
real(kind=dp), parameter :: RPSR = 10.0_dp !PSR radius in random units
real(kind=dp), parameter :: a = 0.60_dp !BH spin parameter

real(kind=dp), parameter :: stheta = 0.0_dp, sphi = 0.0_dp !Minus since rotation matrix defined counterclock
real(kind=dp), parameter :: psi = PI/2.0_dp, chi = 0.0_dp


real(kind=dp), parameter :: rCOM = 20.0_dp !R coordinate of centre of mass
real(kind=dp), parameter :: thetaCOM = PI/2.0_dp !theta coordinate of centre of mass
real(kind=dp), parameter :: phiCOM = PI - PI/6.0_dp !theta coordinate of centre of mass


real(kind=dp), parameter :: N = 1d8 !plasma density normalisation

integer(kind=dp), parameter :: plot = 1 !turn on/off (1/0) numerical accuracy evaluation

end module parameters
