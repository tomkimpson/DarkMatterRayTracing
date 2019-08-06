module initial_conditions


use parameters
use constants
use tensors

implicit none

private rotate_vector, transform_to_global, inverse_transform

public set_initial_conditions


contains








subroutine set_initial_conditions(ki,v,c)

!Arguments
real(kind=dp), dimension(4), intent(inout) :: ki
!real(kind=dp), dimension(3), intent(inout) :: xi, xi_global
real(kind=dp), dimension(6), intent(out) :: v
real(kind=dp), dimension(4), intent(out) :: c

!Other
real(kind=dp) :: r_dot, theta_dot, phi_dot
real(kind=dp) :: r,theta, phi, t, sigma, delta
real(kind=dp) :: pr, ptheta
real(kind=dp) :: E2
real(kind=dp), dimension(4,4) :: metric
real(kind=dp) :: E, kappa, Lz
!play
real(kind=dp) :: xdot, ydot, zdot, dot_mag, mag_4_vel
real(kind=dp) :: Tmag,grr,gtt,gpp, N3,N2,N1, scalar
real(kind=dp) :: AA, w, dr_x, dr_y,dr_z, xp,yp,zp, m, cpt, Eprime
real(kind=dp), dimension(4) :: u_covar, u_contra
real(kind=dp), dimension(4,4) :: metric_contra, metric_covar, transform_matrix
real(kind=dp), dimension(4) :: p_tetrad, p_coordinate
integer(kind=dp) :: i
real(kind=dp) :: tbar, pbar, mm, Enorm, B2,Eprime2
real(kind=dp) :: rdot_prime, thetadot_prime, phidot_prime
real(kind=dp) :: fr, ft, omega2



!Define start point as COM
r = rCOM
theta = thetaCOM
phi = phiCOM
t = 0.0_dp


!Define ray direction
xdot = ki(1)
ydot = ki(2)
zdot = ki(3)


!Define some useful quantities
sigma = r**2 + a**2 * cos(theta)**2
delta = r**2 -2.0_dp*r +a**2
mm = sqrt(r**2 + a**2)




r_dot = mm*r*sin(theta)*cos(phi)*xdot/sigma &
       +mm*r*sin(theta)*sin(phi)*ydot/sigma &
       +mm**2*cos(theta)*zdot/sigma



theta_dot = (mm*cos(theta)*cos(phi) * xdot &
           +mm*cos(theta)*sin(phi) * ydot &
           -r*sin(theta)* zdot&
           )/sigma


phi_dot = (-sin(phi)*xdot + cos(phi)*ydot)/(mm*sin(theta))



!Define plasma frequency

B2 =  N*4.0_dp*PI*electron_charge**2 / electron_mass
call plasma_fr(r,fr)
call plasma_ft(r,ft)
omega2 = B2 * (fr+ft)/sigma



!Compute the energy normalization
!E = 2.0_dp*PI*1.0d9 !1GHz
E = ki(4)
Enorm = (sigma-2.0_dp*r)*(r_dot**2/delta + theta_dot**2) + delta*(sin(theta)*phi_dot)**2
Eprime2 = (E**2 - (sigma-2.0_dp*r)*omega2/sigma)/Enorm 
Eprime = sqrt(Eprime2)




!E2 = (sigma-2.0_dp*r)*(r_dot**2/delta + theta_dot**2 + omega2/sigma) + delta*(sin(theta)*phi_dot)**2
!Enorm = sqrt(E2)




!Eprime = E/Enorm

!Correct for the normalisation
r_dot = r_dot * Eprime
theta_dot = theta_dot*Eprime
phi_dot = phi_dot*Eprime



!Check energies are equal
E2 = (sigma-2.0_dp*r)*(r_dot**2/delta + theta_dot**2 + omega2/sigma) + delta*(sin(theta)*phi_dot)**2

!print *, 'Energy check:', E/ sqrt(E2)





!Define everything else

pr = r_dot * sigma/delta
ptheta = sigma*theta_dot
Lz = (sigma*delta*phi_dot - 2.0_dp*a*r*E)*sin(theta)**2 / (sigma-2.0_dp*r)



!Finish off. Normalise and define kappa



!Normalise to E = 1
pr = pr/E
ptheta = ptheta/E
Lz = Lz/E
B2 = B2 / E**2





!And define a normalized kappa
kappa = ptheta**2 + Lz**2/sin(theta)**2 + a**2*sin(theta)**2


v(1) = r
v(2) = theta
v(3) = phi
v(4) = t
v(5) = pr
v(6) = ptheta


!Dont declare these globally as need to be careful when running in parallel
c(1) = Lz
c(2) = kappa
c(3) = 1.0d-6 !initial stepsize
c(4) = B2



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








subroutine inverse_transform(r,theta,u_covar, u_contra,E,v1,v2)
!Arguments
real(kind=dp), intent(in) :: r,theta,E
real(kind=dp), dimension(4), intent(in) :: u_covar, u_contra
real(kind=dp), dimension(4), intent(in) :: v1
real(kind=dp), dimension(4), intent(out) :: v2

!Other
real(kind=dp), dimension(4,4) :: transform_matrix, metric_contra, metric_covar
real(kind=dp) :: sigma, delta,N1,N2,N3, scalar, Eprime
real(kind=dp) :: p_tetrad(4)
integer(kind=dp) :: i




!Read in p vector in the tetrad frame, calculate Eprime



p_tetrad = v1



call calculate_contravariant_metric(r,theta,metric_contra)
call calculate_covariant_metric(r,theta,metric_covar)


sigma = r**2 + a**2 * cos(theta)**2
delta = r**2 -2.0_dp*r +a**2





N1 = sqrt(- metric_covar(2,2) * (u_covar(1) * u_contra(1) + u_covar(4)*u_contra(4)) * (1.0_dp + u_covar(3)*u_contra(3)) )
N2 = sqrt(metric_covar(3,3) * (1.0_dp + u_covar(3) * u_contra(3)) )
N3 = sqrt(-(u_covar(1) * u_contra(1) + u_covar(4)*u_contra(4))*delta*sin(theta)**2)




scalar = u_covar(2)*u_covar(1)*p_tetrad(2)/N1 &
        +u_covar(3)*u_covar(1)*p_tetrad(3)/N2 &
        - delta*sin(theta)**2*u_contra(4)*p_tetrad(4)/N3







Eprime = (E + scalar)/u_covar(1)
p_tetrad(1) = Eprime




!Construct transformation matrix
call calculate_contravariant_metric(r,theta,metric_contra)
call calculate_covariant_metric(r,theta,metric_covar)


transform_matrix(1,:) = -u_covar

transform_matrix(2,1) = u_covar(2)*u_covar(1)/N1
transform_matrix(2,2) = -metric_covar(2,2)*(u_covar(1)*u_contra(1) + u_covar(3)*u_contra(3))/N1
transform_matrix(2,3) = 0.0_dp
transform_matrix(2,4) = u_covar(2)*u_covar(4)/N1


transform_matrix(3,1) = u_covar(3)*u_covar(1)/N2
transform_matrix(3,2) = u_covar(3)*u_covar(2)/N2
transform_matrix(3,3) = metric_covar(3,3)*(1.0_dp + u_covar(3)*u_contra(3))/N2
transform_matrix(3,4) = u_covar(3)*u_covar(4) / N2


transform_matrix(4,1) = -delta*sin(theta)**2 * u_contra(4) / N3
transform_matrix(4,2) = 0.0_dp
transform_matrix(4,3) = 0.0_dp
transform_matrix(4,4) = delta*sin(theta)**2 * u_contra(1)/N3


do i = 1,4
v2(i)=               transform_matrix(1,i)*p_tetrad(1) + &
                     transform_matrix(2,i)*p_tetrad(2) + &
                     transform_matrix(3,i)*p_tetrad(3) + &
                     transform_matrix(4,i)*p_tetrad(4) 

enddo


end subroutine inverse_transform



end module initial_conditions
