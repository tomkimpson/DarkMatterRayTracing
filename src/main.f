program main


use parameters
use constants
use tensors


implicit none 
real(kind=dp) r,theta, phi

real(kind=dp), dimension(3) :: xi !radiation point w.r.t COM
real(kind=dp), dimension(3) :: xR !radiation point w.r.t BH
real(kind=dp) x, y, z, t, mm
real(kind=dp) Ry_x, Ry_y, Ry_z
real(kind=dp) Rz_x, Rz_y, Rz_z

real(kind=dp), dimension(4) :: u_covar, u_contra !emitter velocity
real(kind=dp), dimension(4,4) :: metric_contra, transform
real(kind=dp) :: N1, N2, N3, delta
!Define emiatter location
!BL coordinates of PSR COM
r = 50.0_dp
theta = PI/2.0_dp
phi = 0.0_dp


!Caresian coords
mm = sqrt(r**2 + a**2)
x = mm * sin(theta)*cos(phi)
y = mm * sin(theta)*sin(phi)
z = mm * cos(theta)

!Location of radiation point w.r.t. COM

xi(1) = RPSR_M*sin(psi)*cos(chi)
xi(2) = RPSR_M*sin(psi)*sin(chi)
xi(3) = RPSR_M*cos(psi)


Ry_x = xi(1)*cos(stheta) + xi(3)*sin(stheta)
Ry_y = xi(2)
Ry_z = xi(1)*(-sin(stheta)) + xi(3)*cos(stheta)

Rz_x = Ry_x*cos(sphi) - Ry_y*sin(sphi)
Rz_y = Ry_x*sin(sphi) + Ry_y*cos(sphi)
Rz_z = Ry_z


xR(1) = Rz_x + x
xR(2) = Rz_y + y
xR(3) = Rz_z + z

!Define emitter 4 velocity
u_covar(1) = 1.0_dp
u_covar(2) = 0.0_dp
u_covar(3) = 0.0_dp
u_covar(4) = 0.0_dp


call calculate_contravariant_metric(r,theta,metric_contra)
call calculate_covariant_metric(r,theta,metric_covar)

u_contra = MATMUL(metric_contra, u_covar)


!Construct matrix to transform to the global (BL) coordinate basis
N1 = - metric


transform(1,:) = u_contra

!Define tangent ray vector in comoving tetrad frame


!Transform to coodinate frame

end program main
