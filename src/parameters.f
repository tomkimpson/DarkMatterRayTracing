module parameters
implicit none



!Define float precision
integer, parameter :: dp = selected_real_kind(33,4931)


real(kind=dp), parameter :: PI = 4.D0*ATAN(1.D0) 
real(kind=dp), parameter :: MBH = 4.310d6 !BH mass in solar  !BH mass
real(kind=dp), parameter :: RPSR = 10.0_dp !PSR radius in random units
real(kind=dp), parameter :: a = 0.60_dp !BH spin parameter

real(kind=dp), parameter :: stheta = -PI/4.0_dp, sphi = PI/2.0_dp !Minus since rotation matrix defined counterclock
real(kind=dp), parameter :: psi = PI/8.0_dp, chi = PI/6.0_dp


real(kind=dp), parameter :: rCOM = 48.0_dp !R coordinate of centre of mass
real(kind=dp), parameter :: thetaCOM = PI/2.0_dp !theta coordinate of centre of mass
real(kind=dp), parameter :: phiCOM = 0.0_dp !theta coordinate of centre of mass

end module parameters
