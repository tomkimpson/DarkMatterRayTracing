module tensors

use parameters
use constants

implicit none


private

public calculate_covariant_metric, calculate_contravariant_metric, magnitude, &
       calculate_minkowski_metric, magnitude_minkowski, mag_3space


contains


subroutine mag_3space(vector,mag)
!arguments
real(kind=dp), intent(in), dimension(3) :: vector
real(kind=dp), intent(out) :: mag


mag = sqrt(vector(1)**2 + vector(2)**2 + vector(3)**2)

end subroutine mag_3space

subroutine magnitude_minkowski(vector,mag)
!Argiuments
real(kind=dp), intent(in), dimension(4) :: vector
real(kind=dp), intent(out) :: mag


mag = -1.0_dp *vector(1)**2.0_dp + &
       1.0_dp*vector(2)**2.0_dp + &
       1.0_dp*vector(3)**2.0_dp + &
       1.0_dp*vector(4)**2.0_dp 



 end subroutine magnitude_minkowski


subroutine magnitude(metric,vector,mag)
!Arguments
real(kind=dp), intent(IN), dimension(4,4) :: metric
real(kind=dp), intent(IN), dimension(4) :: vector
real(kind=dp), intent(out) :: mag



mag = metric(1,1)*vector(1)**2.0_dp + metric(2,2)*vector(2)**2.0_dp + metric(3,3)*vector(3)**2.0_dp+metric(4,4)*vector(4)**2.0_dp&
      +2.0_dp*(metric(1,2)*vector(1)*vector(2) + & 
               metric(1,3)*vector(1)*vector(3) + &
               metric(1,4)*vector(1)*vector(4) + &
               metric(2,3)*vector(2)*vector(3) + &
               metric(2,4)*vector(2)*vector(4) + &
               metric(3,4)*vector(3)*vector(4)   &
               )





end subroutine magnitude





subroutine calculate_contravariant_metric(r,theta,metric)
!Arguments
real(kind=dp), intent(IN) :: r, theta
real(kind=dp), intent(out), dimension(4,4) :: metric


!Internals
real(kind=dp) :: sigma, delta


sigma = r**2.0_dp + a**2.0_dp*cos(theta)**2.0_dp
delta = r**2.0_dp - 2.0_dp*r + a**2.0_dp


metric(1,1) = - ((r**2.0_dp + a**2.0_dp) + 2.0_dp*r*a**2.0_dp*sin(theta)**2.0_dp/sigma)/delta
metric(2,2) = delta/sigma
metric(3,3) = 1.0_dp/sigma
metric(4,4) = (1.0_dp-2.0_dp*r/sigma)/(delta*sin(theta)**2.0_dp)

metric(1,4) = -2.0_dp*r*a/(sigma*delta)
metric(4,1) = metric(1,4)

!All other terms are zero
metric(1,2) = 0.0_dp
metric(1,3) = 0.0_dp

metric(2,1) = 0.0_dp
metric(2,3) = 0.0_dp
metric(2,4) = 0.0_dp

metric(3,1) = 0.0_dp
metric(3,2) = 0.0_dp
metric(3,4) = 0.0_dp

metric(4,2) = 0.0_dp
metric(4,3) = 0.0_dp


end subroutine calculate_contravariant_metric


subroutine calculate_covariant_metric(r,theta,metric)
!Arguments
real(kind=dp), intent(IN) :: r, theta
real(kind=dp), intent(out), dimension(4,4) :: metric


!Internals
real(kind=dp) :: sigma, delta


sigma = r**2.0_dp + a**2.0_dp*cos(theta)**2.0_dp
delta = r**2.0_dp - 2.0_dp*r + a**2.0_dp



metric(1,1) = -(1.0_dp-2.0_dp*r/sigma)
metric(2,2) = sigma/delta
metric(3,3) = sigma
metric(4,4) = (r**2.0_dp + a**2.0_dp + 2.0_dp*r*a**2.0_dp*sin(theta)**2.0_dp/sigma)*sin(theta)**2.0_dp


metric(1,4) = -2.0_dp*r*a*sin(theta)**2.0_dp/sigma
metric(4,1) = metric(1,4)

!All other terms are zero
metric(1,2) = 0.0_dp
metric(1,3) = 0.0_dp

metric(2,1) = 0.0_dp
metric(2,3) = 0.0_dp
metric(2,4) = 0.0_dp

metric(3,1) = 0.0_dp
metric(3,2) = 0.0_dp
metric(3,4) = 0.0_dp

metric(4,2) = 0.0_dp
metric(4,3) = 0.0_dp


end subroutine calculate_covariant_metric


subroutine calculate_minkowski_metric(r,theta,metric)
!Note this is the covariant form
!Contravariant form is just metric**-1
        
!Arguments
real(kind=dp), intent(IN) :: r, theta
real(kind=dp), intent(out), dimension(4,4) :: metric


metric = 0.0_dp

metric(1,1) = -1.0_dp
metric(2,2) = 1.0_dp
metric(3,3) = r**2
metric(4,4) = r**2 * sin(theta)**2


end subroutine calculate_minkowski_metric



end module tensors


