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
real(kind=dp), dimension(4,4) :: metric_contra, transform, metric_covar, transform_to
real(kind=dp) :: N1, N2, N3, delta
real(kind=dp),dimension(4) :: k, kp, k_contra, k_covar
real(kind=dp), dimension(4) :: k_contra_global, k_covar_global
real(kind=dp) :: mag1, mag2, cpt, sigma, summation


integer(kind=dp)  :: ii, jj, kk
real(kind=dp) :: u0,u1,u2,u3, grr, gthth
real(kind=dp) :: u0_covar,u1_covar,u2_covar,u3_covar
real(kind=dp), dimension(4,4) :: Bmatrix

!Define emiatter location
!BL coordinates of PSR COM


!Caresian coords
!mm = sqrt(r**2 + a**2)
!x = mm * sin(theta)*cos(phi)
!y = mm * sin(theta)*sin(phi)
!z = mm * cos(theta)

!Location of radiation point w.r.t. COM

!xi(1) = RPSR_M*sin(psi)*cos(chi)
!xi(2) = RPSR_M*sin(psi)*sin(chi)
!xi(3) = RPSR_M*cos(psi)


!Ry_x = xi(1)*cos(stheta) + xi(3)*sin(stheta)
!Ry_y = xi(2)
!Ry_z = xi(1)*(-sin(stheta)) + xi(3)*cos(stheta)

!Rz_x = Ry_x*cos(sphi) - Ry_y*sin(sphi)
!Rz_y = Ry_x*sin(sphi) + Ry_y*cos(sphi)
!Rz_z = Ry_z


!xR(1) = Rz_x + x
!xR(2) = Rz_y + y
!xR(3) = Rz_z + z




!Define emitter location (centre of mass )
r = 50.0_dp
theta = PI/2.0_dp
phi = 0.0_dp


!Define metric
call calculate_contravariant_metric(r,theta,metric_contra)
call calculate_covariant_metric(r,theta,metric_covar)






!Check the metrics

!ii = 1
!print *, metric_contra(1,ii) * metric_covar(ii,1) + &
 !        metric_contra(2,ii) * metric_covar(ii,2) + &
  !       metric_contra(3,ii) * metric_covar(ii,3) + &
   !      metric_contra(4,ii) * metric_covar(ii,4) 

! stop



!Define emitter 4 velocity
cpt = metric_contra(1,1) + 2.0_dp*metric_contra(1,4) + metric_contra(4,4) 
cpt = sqrt(-1.0_dp/cpt)

!Ensures the normalization is correct
u_covar(1) = cpt
u_covar(2) = 0.0_dp
u_covar(3) = 0.0_dp
u_covar(4) = cpt


call magnitude(metric_contra, u_covar, mag1)
print *, '4-velocity magnitude:',  mag1
u_contra = MATMUL(metric_contra, u_covar)
call magnitude(metric_covar, u_contra, mag1)
print *, '4-velocity magnitude:',  mag1


!Construct some random vector. 
!This vector is in the tetrad frame and so is lowered/raised by Minkowski
k_contra(1) = 1.0_dp
k_contra(2) = 2.0_dp
k_contra(3) = 3.0_dp
k_contra(4) = 4.0_dp


k_covar(1) = -k_contra(1)
k_covar(2) = k_contra(2)
k_covar(3) = k_contra(3)
k_covar(4) = k_contra(4)


!The size of this vector is
call magnitude_minkowski(k_contra, mag1)
print *, 'Vector size in comoving frame', mag1


call magnitude_minkowski(k_covar, mag1)
print *, 'Vector size in comoving frame', mag1



!call magnitude(metric_covar, k_contra, mag1)
!print *, mag1
!stop

!Construct matrix to transform to the global (BL) coordinate basis
delta = r**2.0_dp + a**2.0_dp - 2.0_dp*r
N1 = sqrt(- metric_covar(2,2) * (u_covar(1) * u_contra(1) + u_covar(4)*u_contra(4)) * (1.0_dp + u_covar(3)*u_contra(3)) )
N2 = sqrt(metric_covar(3,3) * (1.0_dp + u_covar(3) * u_contra(3)) )
N3 = sqrt(-(u_covar(1) * u_contra(1) + u_covar(4)*u_contra(4))*delta*sin(theta)**2)

transform(1,:) = u_contra

transform(2,1) = u_covar(2)*u_contra(1)/N1 
transform(2,2) = -(u_covar(1) * u_contra(1) + u_covar(4)*u_contra(4))/N1
transform(2,3) = 0.0_dp
transform(2,4) = u_covar(2)*u_contra(4)/N1


transform(3,1) = u_covar(3)*u_contra(1)/N2
transform(3,2) = u_covar(3)*u_contra(2) / N2
transform(3,3) = (1.0_dp + u_covar(3)*u_contra(3))/N2
transform(3,4) = u_covar(3)*u_contra(4)/N2



transform(4,1) = u_covar(4)/N3
transform(4,2) = 0.0_dp
transform(4,3) = 0.0_dp
transform(4,4) = -u_covar(1)/N3








Bmatrix(1,:) = -u_covar

Bmatrix(2,1) = u_covar(2)*u_covar(1)/N1
Bmatrix(2,2) = -metric_covar(2,2)*(u_covar(1)*u_contra(1) + u_covar(4)*u_contra(4))/N1
Bmatrix(2,3) = 0.0_dp
Bmatrix(2,4) = u_covar(2)*u_covar(4)/N1


Bmatrix(3,1) = u_covar(3)*u_covar(1)/N2
Bmatrix(3,2) = u_covar(3)*u_covar(2)/N2
Bmatrix(3,3) = metric_covar(3,3)*(1.0_dp + u_covar(3)*u_contra(3))/N2
Bmatrix(3,4) = u_covar(3)*u_covar(4)/N2


Bmatrix(4,1) = u_contra(4)
Bmatrix(4,2) = 0.0_dp
Bmatrix(4,3) = 0.0_dp
Bmatrix(4,4) = -u_contra(1)

Bmatrix(4,:) = -Bmatrix(4,:) * delta*sin(theta)**2/N3




!print *, 'Check transform matrix'
!do ii = 1,4
!print *, transform(ii,1)*Bmatrix(ii,1) + &
 !       transform(ii,2)*Bmatrix(ii,2) + &
  !       transform(ii,3)*Bmatrix(ii,3) + &
   !     transform(ii,4)*Bmatrix(ii,4) 


!enddo









k_contra_global(1) = transform(1,1) * k_contra(1) + &
                     transform(2,1) * k_contra(2) + &
                     transform(3,1) * k_contra(3) + &
                     transform(4,1) * k_contra(4)


print *, k_contra_global(1)


k_covar_global(1) = Bmatrix(1,1) * k_covar(1) + &
                    Bmatrix(2,1) * k_covar(2) + &
                    Bmatrix(3,1) * k_covar(3) + &
                    Bmatrix(4,1) * k_covar(4) 

print *, k_covar_global(1)*metric_contra(1,1)


stop




!call magnitude_minkowski(k_contra_global, mag1)




!Transform to global coordinate frame

k_contra_global = MATMUL(transform, k_contra)
k_covar_global = MATMUL(Bmatrix, k_covar)






print *, 'Checks'
print *, 'Global Vectors:'
print *, k_covar_global
print *, k_contra_global






!print *, 'covar metric'
!print *, metric_covar(1,:)
!print *, metric_covar(2,:)
!print *, metric_covar(3,:)
!print *, metric_covar(4,:)




print *, 'compare'
ii = 3
print *, k_covar_global(ii)
print *, k_contra_global(ii)*metric_covar(ii,ii)



stop



print *, 'sum magnitude'



print *, k_contra_global(1)*k_covar_global(1) &
         + k_contra_global(2)*k_covar_global(2) &
         + k_contra_global(3)*k_covar_global(3) &
         + k_contra_global(4)*k_covar_global(4) 








print *, 'Check Kronekar'
do ii = 1,4
do kk = 1,4

 summation = 0
 do jj = 1,4
 summation = summation + metric_contra(ii,jj)*metric_covar(jj,kk)
 enddo
 print *, ii,kk,summation





 enddo


 enddo








!print *, k_contra_global(ii)
!print *, k_covar_global(1)*metric_contra(1,ii) + & 
 !        k_covar_global(2)*metric_contra(2,ii) + &
  !       k_covar_global(3)*metric_contra(3,ii) + &
   !      k_covar_global(4)*metric_contra(4,ii) 

 
 



!print *, k_contra_global(1) 
!print *, metric_contra(1,1) * k_covar_global(1) + metric_contra(1,4) * k_covar_global(4)


print *, k_contra_global(1)*k_covar_global(1) &
         + k_contra_global(2)*k_covar_global(2) &
         + k_contra_global(3)*k_covar_global(3) &
         + k_contra_global(4)*k_covar_global(4) 



!The magnitude is now given by the full metric
!call magnitude_minkowski(k_contra_global, mag1)
call magnitude(metric_covar, k_contra_global, mag1)
print *, 'Vector size in global frame', mag1


end program main
