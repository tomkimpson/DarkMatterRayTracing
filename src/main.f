program main


use parameters
use constants
use tensors
use initial_conditions

implicit none 

real(kind=dp),dimension(3) :: ki !Comoving beam direction vector. Magnitude 1.
real(kind=dp), dimension(3) :: xi !local vector location
real(kind=dp), dimension(3) :: xi_global !global vector location

real(kind=dp) :: xglobal, yglobal, zglobal, w

real(kind=dp) :: rmag

!Define tangent vector in cartesian components 
ki(1) = 1.0_dp*sin(psi)*cos(chi) !xdot
ki(2) = 1.0_dp*sin(psi)*sin(chi) !ydot
ki(3) = 1.0_dp*cos(psi) !zdot

!Define location of vector in local coordinates
rmag = RPSR_M
xi(1) = rmag*sin(psi)*cos(chi) !xdot
xi(2) = rmag*sin(psi)*sin(chi) !ydot
xi(3) = rmag*cos(psi) !zdot


!Define location of vector in global coordinates (COM + flat ST.)
xglobal = m*sin(thetaCOM)*cos(phiCOM)! + xi(1)
yglobal = m*sin(thetaCOM)*sin(phiCOM)! + xi(2)
zglobal = rCOM*cos(thetaCOM)! + xi(3)


w = xglobal**2 + yglobal**2 + zglobal**2 - a**2

xi_global(1) = sqrt(0.50_dp*(w + sqrt(w**2 + 4*a**2*zglobal**2)))
xi_global(2) = acos(zglobal/xi_global(1))
xi_global(3) = atan2(yglobal,xglobal)

xi_global(1) = rCOM
xi_global(2) = thetaCOM
xi_global(3) = phiCOM

call set_initial_conditions(ki,xi,xi_global)

end program main
