module derivatives

use parameters
use constants
use metric
implicit none

public geodesic, DM_gr

private DM_gr_deriv

contains


subroutine geodesic(v,constants,dv)


        !Arguments
real(kind=dp), dimension(6) :: v,dv
real(kind=dp), dimension(4),intent(in) :: constants
!Other
real(kind=dp) :: r,theta,phi,t,pr,ptheta
real(kind=dp) :: sigma,delta, SD,csc
real(kind=dp) :: Lz, kappa,B2, plasma_correction_pr, plasma_correction_pt
real(kind=dp) :: fr, fr_prime, ft, ft_prime
real(kind=dp) :: Gr, upsilon, d_gr
real(kind=dp) :: delta_prime, psi

!Read in coordinate variables
r= v(1)
theta= v(2)
phi= v(3)
t= v(4)
pr= v(5)
ptheta= v(6)

!Get the constants
Lz = constants(1)
kappa = constants(2)
B2 = constants(3)

!Define some useful quantities


sigma = r**2 + a**2 * cos(theta)**2


if (DM .EQ. 1) then
call DM_gr(r,Gr)
call DM_gr_deriv(r,d_gr)
delta = r**2 * Gr + a**2
upsilon = r*(1.0_dp - Gr)/2.0_dp

delta_prime = r**2 * d_gr + 2.0_dp*r*Gr
psi = r*(2.0_dp - 2.0_dp*Gr - r*d_gr)/2.0_dp
else
Gr = 1.0_dp - 2.0_dp/r
d_gr = 2.0_dp*r**(-2.0_dp)
delta = r**2 -2.0_dp*r + a**2
upsilon = 1.0_dp
delta_prime = 2.0_dp*(r - 1.0_dp)
psi = 1.0_dp
endif


csc = 1.0_dp /sin(theta)
SD = sigma*delta




dv(1) = pr*delta/sigma
dv(2) = ptheta/sigma
dv(3) = (2.0_dp*a*r*upsilon + (sigma - 2.0_dp*r*upsilon)*Lz*csc**2 ) / SD
dv(4) = 1.0_dp + (2.0_dp*r*upsilon*(r**2+a**2)-2.0_dp*a*r*upsilon*Lz)/SD
dv(5) = (-kappa*delta_prime/2.0_dp + 2.0_dp*r*(r**2 +a**2) - 2.0_dp*a*Lz*psi)/SD - pr**2*delta_prime/sigma
dv(6) = sin(theta)*cos(theta)*(Lz**2 * csc**4 -a**2) / sigma 

 


end subroutine geodesic





subroutine DM_gr(r,gr)
!Arguments
real(kind=dp), intent(in) :: r
real(kind=dp), intent(out) :: gr
!Other
real(kind=dp) :: BigK, rr,AA



rr = r/r0
BigK = -2.0_dp*kDM*PI/r
AA = 1.0_dp + rr



gr = (1.0_dp + rr**2)**(BigK*(1.0_dp-rr)) * &
     AA**(2.0_dp*BigK*AA) * &
     exp(-2.0_dp*BigK*AA*atan(rr)) -&
     2.0_dp/r





print *, 'GR'
print *, r,BigK, (1.0_dp-rr), r0, r0_SI
stop

!print *, 'GR:'
!print *, (1.0_dp + rr**2)**(BigK*(1.0_dp-rr))
!print *, AA**(2.0_dp*BigK*AA)
!print *, exp(-2.0_dp*BigK*AA*atan(rr)) 
!print *, '-------'



end subroutine DM_gr

subroutine DM_gr_deriv(r,d_gr)
!Argument
real(kind=dp),intent(in) :: r 
real(kind=dp),intent(out) :: d_gr
!Other
real(kind=dp) :: rr, BigK, AA

rr = r/r0
BigK = -2.0_dp*kDM*PI/r
AA = 1.0_dp + rr



d_gr = 2.0_dp*(1.0_dp - exp(-2.0_dp*BigK*AA*atan(rr)) * &
              kDM*PI*(1+rr**2)**(BigK*(1-rr)) * &
              AA**(2.0_dp*BigK*AA)* (2.0_dp*atan(rr) - log(1+rr**2) - 2.0_dp*log(AA)) &
        )/r**2 







end subroutine DM_gr_deriv




end module derivatives
