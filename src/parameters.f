module parameters
implicit none



!Define float precision
integer, parameter :: dp = selected_real_kind(33,4931)

!Universal constants
real(kind=dp), parameter :: PI = 4.D0*ATAN(1.D0) 


!BH parameters
real(kind=dp), parameter :: MBH = 4.310d6 !BH mass in solar  !BH mass
real(kind=dp), parameter :: a = 0.9980_dp !BH spin parameter


!PSR parameters
real(kind=dp), parameter :: RPSR = 10.0_dp !PSR radius in random units
real(kind=dp), parameter :: stheta = 0.0_dp, sphi = 0.0_dp !Minus since rotation matrix defined counterclock
real(kind=dp), parameter :: psi = PI/2.0_dp, chi = 0.0_dp


real(kind=dp), parameter :: xCOM = -1000.0_dp
real(kind=dp), parameter :: yCOM = +7.0_dp
real(kind=dp), parameter :: zCOM = 00.0_dp

real(kind=dp), parameter :: rC = 10.0_dp !R coordinate of centre of mass
real(kind=dp), parameter :: thetaC = PI/2.0_dp !theta coordinate of centre of mass
real(kind=dp), parameter :: phiC = PI/6.0_dp !theta coordinate of centre of mass

!Use Cartesian/BL coords for setting initial PSR location
integer(kind=dp), parameter :: PSR_Location = 1




!Plasma paramters
real(kind=dp), parameter :: N = 0.00_dp !plasma density normalisation


!I/O options
integer(kind=dp), parameter :: plot = 1 !turn on/off (1/0) numerical accuracy evaluation

end module parameters
