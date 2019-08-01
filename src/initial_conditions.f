module initial_conditions


use parameters
use constants
use tensors

implicit none

private rotate_vector, transform_to_global

public set_initial_conditions


contains








subroutine set_initial_conditions(ki,xi,xi_global,v)

!Arguments
real(kind=dp), dimension(3), intent(inout) :: ki
real(kind=dp), dimension(3), intent(inout) :: xi, xi_global
real(kind=dp), dimension(6), intent(out) :: v
!Other
real(kind=dp) :: r_dot, theta_dot, phi_dot
real(kind=dp) :: r,theta, phi, t, sigma, delta
real(kind=dp) :: pr, ptheta
real(kind=dp) :: E2


!play
real(kind=dp) :: xdot, ydot, zdot
!Rotate vector components and vector location


call rotate_vector(ki)
call rotate_vector(xi)


!Get vector in global frame
call transform_to_global(ki,xi,xi_global)


!Define start point as COM
r = xi_global(1)
theta = xi_global(2)
phi = xi_global(3)
t = 0.0_dp


sigma = r**2 + a**2 * cos(theta)**2
delta = r**2 -2.0_dp*r +a**2



r_dot = ki(1)
theta_dot = ki(2)
phi_dot = ki(3)




xdot = 1.0_dp
ydot = 0.0_dp
zdot = 0.0_dp


r_dot = xdot*sin(theta)*cos(phi) + ydot*sin(theta)*sin(phi) + zdot*cos(theta)

theta_dot = -(-r_dot*cos(theta)/(r*sin(theta)) + zdot/(r*sin(theta)))

phi_dot = (-xdot*sin(phi) + ydot*cos(phi)) / (r*sin(theta))



!r_dot = cos(phi)
!theta_dot = 0.0_dp
!phi_dot = -sin(phi)



print *, 'BL dots:', r_dot, theta_dot, phi_dot


!xdot = r_dot*cos(phi) - r*phi_dot*sin(phi)
!ydot = r_dot*sin(phi) + r*phi_dot*cos(phi)
!zdot = -theta_dot


print *, 'Cartesian dots', xdot, ydot,zdot



!print *, r_dot, theta_dot, phi_dot


pr = r_dot * sigma/delta
ptheta = sigma*theta_dot



!Compute the Energy and angular momentum (i.e. pt, phi)
E2 = (sigma-2.0_dp*r)*(r_dot**2/delta + theta_dot**2) + delta*(sin(theta)*phi_dot)**2
E = sqrt(E2)
Lz = (sigma*delta*phi_dot - 2.0_dp*a*r*E)*sin(theta)**2 / (sigma-2.0_dp*r)

!Normalise to E = 1
pr = pr/E
ptheta = ptheta/E
Lz = Lz/E

!And define a normalized kappa
kappa = ptheta**2 + Lz**2/sin(theta)**2 + a**2*sin(theta)**2


v(1) = r
v(2) = theta
v(3) = phi
v(4) = t
v(5) = pr
v(6) = ptheta





!print *, r_dot, theta_dot, phi_dot
!print *, v(1:3)
!print *, v(4:6)

!print *, E, Lz, kappa

!stop





end subroutine set_initial_conditions



subroutine rotate_vector(ki)
!Rotate the direction vector to account for spin-axis alignment

!Arguments
real(kind=dp), dimension(3), intent(inout) :: ki

!Other
real(kind=dp), dimension(3,3) :: Rz, Ry
real(kind=dp), dimension(3) :: kprime

!Can also do rotation matrixes here.
!I prefer the algebra, just to make the coord transform explicit.
!Rz = 0.0_dp
!Ry= 0.0_dp

!Rz(1,1) = cos(sphi)
!Rz(1,2) = -sin(sphi)
!Rz(2,1) = sin(sphi)
!Rz(2,2) = cos(sphi)
!Rz(3,3) = 1.0_dp

!Ry(1,1) = cos(stheta)
!Ry(1,3) = sin(stheta)
!Ry(3,1) = -sin(stheta)
!Ry(3,3) = cos(stheta)
!Ry(2,2) = 1.0_dp

!kprime = MATMUL(Rz,MATMUL(Ry,ki))





kprime(1) = ki(1) * (cos(stheta)*cos(sphi)) + ki(2)*(-sin(sphi)) + ki(3)*(sin(stheta)*cos(sphi))
kprime(2) = ki(1) * (cos(stheta)*sin(sphi)) + ki(2)*(cos(sphi)) + ki(3)*(sin(stheta)*sin(sphi))
kprime(3) = ki(1)*(-sin(stheta)) + ki(2)*(0.0_dp) + ki(3)*(cos(stheta))


ki = kprime

end subroutine rotate_vector



subroutine transform_to_global(ki,xi,xi_global)

!Arguments        
real(kind=dp), dimension(3), intent(inout) :: ki
real(kind=dp), dimension(3), intent(in) :: xi, xi_global
!Other
real(kind=dp) :: r, theta, phi
real(kind=dp), dimension(4,4) :: metric_contra, metric_covar, transform_matrix, metric_minkowski
real(kind=dp), dimension(4) :: u_covar, u_contra, k_contra, k_covar
real(kind=dp) :: cpt, magnitudeTetrad,magnitudeGlobal
real(kind=dp) :: delta, N1, N2, N3, mag_4_vel
real(kind=dp), dimension(4) :: k_contra_global
real(kind=dp), dimension(3,3) :: jacobian
integer(kind=dp) :: i

!Play
real(kind=dp) :: mag1, rr,tt, pp

!Transform to spherical polar basis
rr = sqrt(xi(1)**2 + xi(2)**2 + xi(3)**2)
tt = acos(xi(3)/rr)
pp = atan(xi(2)/xi(1))


call mag_3space(ki,mag1)
!print *, 'Ori Magnitude = ', mag1



jacobian(1,1) = sin(tt)*cos(pp)
jacobian(1,2) = sin(tt)*sin(pp)
jacobian(1,3) = cos(tt)

jacobian(2,1) = cos(tt)*cos(pp)
jacobian(2,2) = cos(tt)*sin(pp)
jacobian(2,3) = - sin(tt)

jacobian(3,1) = -sin(pp)
jacobian(3,2) = cos(pp)
jacobian(3,3) = 0.0_dp

!Transform vector from cartesian to spherical coords.
k_contra(1) = 0.0_dp
k_contra(2:4) = MATMUL(jacobian, ki)


print *, ki

call mag_3space(k_contra(2:4),mag1)

print *, k_contra(2:4)


!Define metric
r = xi_global(1)
theta = xi_global(2)
phi = xi_global(3)



call calculate_contravariant_metric(r,theta,metric_contra)
call calculate_covariant_metric(r,theta,metric_covar)



!Define the 4-velocity - ultimately this will be an argument. For now just define
cpt = metric_contra(1,1) + 2.0_dp*metric_contra(1,4) + metric_contra(4,4) 
cpt = sqrt(-1.0_dp/cpt)

u_covar(1) = cpt
u_covar(2) = 0.0_dp
u_covar(3) = 0.0_dp
u_covar(4) = cpt
u_contra = MATMUL(metric_contra, u_covar)


mag_4_vel =u_covar(1)*u_contra(1) + &
           u_covar(2)*u_contra(2) + &
           u_covar(3)*u_contra(3) + &
           u_covar(4)*u_contra(4) 



print *, 'Magnitude of 4-velocity should be negative unity:', mag_4_vel







!Knowing the magntiude of the vector will be a nice check after the transform



!Now construct the transform matrix

delta = r**2.0_dp + a**2.0_dp - 2.0_dp*r
N1 = sqrt(- metric_covar(2,2) * (u_covar(1) * u_contra(1) + u_covar(4)*u_contra(4)) * (1.0_dp + u_covar(3)*u_contra(3)) )
N2 = sqrt(metric_covar(3,3) * (1.0_dp + u_covar(3) * u_contra(3)) )
N3 = sqrt(-(u_covar(1) * u_contra(1) + u_covar(4)*u_contra(4))*delta*sin(theta)**2)


transform_matrix(1,:) = u_contra

transform_matrix(2,1) = u_covar(2)*u_contra(1)/N1 
transform_matrix(2,2) = -(u_covar(1) * u_contra(1) + u_covar(4)*u_contra(4))/N1
transform_matrix(2,3) = 0.0_dp
transform_matrix(2,4) = u_covar(2)*u_contra(4)/N1


transform_matrix(3,1) = u_covar(3)*u_contra(1)/N2
transform_matrix(3,2) = u_covar(3)*u_contra(2) / N2
transform_matrix(3,3) = (1.0_dp + u_covar(3)*u_contra(3))/N2
transform_matrix(3,4) = u_covar(3)*u_contra(4)/N2



transform_matrix(4,1) = u_covar(4)/N3
transform_matrix(4,2) = 0.0_dp
transform_matrix(4,3) = 0.0_dp
transform_matrix(4,4) = -u_covar(1)/N3



!Note this is not just a MatMul. See Kulkarni et al.

do i = 1,4

k_contra_global(i) = transform_matrix(1,i)*k_contra(1) + &
                     transform_matrix(2,i)*k_contra(2) + &
                     transform_matrix(3,i)*k_contra(3) + &
                     transform_matrix(4,i)*k_contra(4) 
enddo



!print *, k_contra
!print *, k_contra_global


!Now get magntiude
call magnitude(metric_covar, k_contra_global,magnitudeGlobal)



!print *, magnitudeGlobal



ki = k_contra_global(2:4)


end subroutine transform_to_global


end module initial_conditions
