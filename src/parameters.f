module parameters
implicit none



!Define float precision
integer, parameter :: dp = selected_real_kind(33,4931)


real(kind=dp), parameter :: PI = 4.D0*ATAN(1.D0) 
real(kind=dp), parameter :: MBH = 4.310d6 !BH mass in solar  !BH mass
real(kind=dp), parameter :: RPSR = 10.0_dp !PSR radius in km
real(kind=dp), parameter :: a = 0.60_dp !BH spin parameter

real(kind=dp), parameter :: stheta = 0.0_dp, sphi = 0.0_dp
real(kind=dp), parameter :: psi = PI/2.0_dp, chi = 0.0_dp




end module parameters
