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
real(kind=dp), parameter :: Rhor = 1.0_dp+sqrt(1.0_dp-a**2) + 1.0d-2 !Horizon + eps
real(kind=dp), PARAMETER :: electron_charge = 4.80320425D-10 !CGS
real(kind=dp), PARAMETER :: electron_mass = 9.10938356D-28 !CGS


!Cash-Karp parameters and integration constants
real(kind = dp) :: B21=1.0_dp/5.0_dp
real(kind = dp) :: B31 = 3.0_dp/40.0_dp , B32 = 9.0_dp/40.0_dp
real(kind = dp) :: B41 = 3.0_dp/10.0_dp, B42 = -9.0_dp/10.0_dp, B43 = 6.0_dp/5.0_dp 
real(kind = dp) :: B51 = -11.0_dp/54.0_dp, B52 = 5.0_dp/2.0_dp, B53 = -70.0_dp/27.0_dp
real(kind = dp) :: B54 = 35.0_dp/27.0_dp
real(kind = dp) :: B61 = 1631.0_dp/55296.0_dp, B62 = 175.0_dp/512.0_dp, B63 = 575.0_dp/13824.0_dp
real(kind = dp) :: B64 = 44275.0_dp/110592.0_dp, B65 = 253.0_dp/4096.0_dp


real(kind = dp) :: c1 = 37.0_dp/378.0_dp, c3 = 250.0_dp/621.0_dp, c4 = 125.0_dp/594.0_dp
real(kind = dp) :: c6=512.0_dp/1771.0_dp
real(kind = dp) :: cbar1 = 2825.0_dp/27648.0_dp, cbar3 = 18575.0_dp/48384.0_dp
real(kind = dp) :: cbar4=13525.0_dp/55296.0_dp, cbar5 = 277.0_dp/14336.0_dp, cbar6 = 1.0_dp/4.0_dp

real(kind = dp) :: escal = 1.0d19
real(kind = dp) :: PSHRINK = -0.25_dp,PGROW=-0.20_dp, S=1.0_dp






!Global constants declared later
character(len=300) :: path
!real(kind=dp) :: dh = 1.0d-6
!real(kind=dp) :: E, Lz, kappa
real(kind=dp), dimension(6) :: delta = 0.0_dp
real(kind=dp) :: rCOM, thetaCOM, phiCOM

end module constants
