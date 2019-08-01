module numerical_performance


use parameters
use constants
use tensors

implicit none

private

public eval_performance




contains




subroutine eval_performance(v,c,dk)
!Arguments
real(kind=dp),intent(in), dimension(6) :: v !Variables (r,theta,phi,t,pr,ptheta)
real(kind=dp),intent(in), dimension(3) :: c

real(kind=dp),intent(out) :: dk
!Other
real(kind=dp) :: Lz, kappa
real(kind=dp) :: theta, ptheta, kappa_check, r,pr
real(kind=dp), dimension(4,4) :: metric
r = v(1)
theta = v(2)
pr = v(5)
ptheta = v(6)

Lz = c(1)
kappa = c(2)

kappa_check = ptheta**2 + Lz**2/sin(theta)**2 +a**2*sin(theta)**2

dk = (kappa_check - kappa) / kappa
end subroutine eval_performance




end module numerical_performance
