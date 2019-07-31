program main


use parameters
use constants
use tensors
use initial_conditions

implicit none 

real(kind=dp),dimension(3) :: ki !Comoving beam direction vector. Magnitude 1.
real(kind=dp), dimension(3) :: xi !local vector location
real(kind=dp), dimension(3) :: xi_global !global vector location
real(kind=dp), dimension(6) :: v !Variables (r,theta,phi,t,pr,ptheta)

real(kind=dp) :: xglobal, yglobal, zglobal, w

real(kind=dp) :: rmag, mag_check1, mag_check2

!Define tangent vector in cartesian components 
ki(1) = 1.0_dp*sin(psi)*cos(chi) !xdot
ki(2) = 1.0_dp*sin(psi)*sin(chi) !ydot
ki(3) = 1.0_dp*cos(psi) !zdot

!Define location of vector in local coordinates
rmag = RPSR_M
xi(1) = rmag*sin(psi)*cos(chi) !xdot
xi(2) = rmag*sin(psi)*sin(chi) !ydot
xi(3) = rmag*cos(psi) !zdot


!For now we take the 'location' of the vector as just being the COM.
!It is an open question 

xi_global(1) = rCOM
xi_global(2) = thetaCOM
xi_global(3) = phiCOM


call set_initial_conditions(ki,xi,xi_global,v)



end program main
