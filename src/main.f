program main


use parameters
use constants
use tensors
use initial_conditions
use integrate
use numerical_performance




implicit none 

real(kind=dp),dimension(3) :: ki !Comoving beam direction vector. Magnitude 1.
real(kind=dp), dimension(3) :: xi !local vector location
real(kind=dp), dimension(3) :: xi_global !global vector location
real(kind=dp), dimension(6) :: v !Variables (r,theta,phi,t,pr,ptheta)

real(kind=dp) :: xglobal, yglobal, zglobal, w

real(kind=dp) :: rmag, mag_check1, mag_check2, dk

integer(kind=dp) :: counter, j,i
real(kind=dp), dimension(int(1e6),6) :: AllData !Big array to save all data
real(kind=dp) :: mm, xC, yC, zC
real(kind=dp) :: dummy_angle,res
character(len=200) :: ID
!Set up save location
call get_environment_variable("RayTracingDir", path)


!print *, real(1e6)*real(6)*real(dp)/1.0d9, ' GB'


 !!!!!! Temporarily excluded -------



!Define tangent vector in cartesian components 
!ki(1) = 1.0_dp*sin(psi)*cos(chi) !xdot
!ki(2) = 1.0_dp*sin(psi)*sin(chi) !ydot
!ki(3) = 1.0_dp*cos(psi) !zdot




!Define location of vector in local coordinates
!rmag = RPSR_M
!xi(1) = rmag*sin(psi)*cos(chi) !xdot
!xi(2) = rmag*sin(psi)*sin(chi) !ydot
!xi(3) = rmag*cos(psi) !zdot


 !!!!!! Temporarily excluded -------




!!!!!$OMP PARALLEL DO PRIVATE(ki, v, &
!!!!!!$OMP& counter,AllData,xC,yC,zC    )






res = 20.0_dp
do j = 0,int(res)
dummy_angle = (PI/res) * real(j,kind=dp)


ki(1) = sin(dummy_angle) !xdot
ki(2) = cos(dummy_angle) !ydot
ki(3) = 0.0_dp !zdot


call set_initial_conditions(ki,v)


counter = 1


do while (v(1) .GT. Rhor .and. v(1) .LT. 40.0_dp)



call rk(v)


if (plot .EQ. 1) then
!Save data at each integration point for plotting - probably not necessary for PSR timing (where we only care about endpoints)
AllData(counter,:) = v

endif
counter = counter + 1

enddo




if (plot .EQ. 1) then


write( ID, '(f10.2)' )  dummy_angle

print *, ID

open(unit = 10, file = trim(adjustl(path))//'RT_'//trim(adjustl(ID))//'.txt', status = 'replace', form='formatted')



do i=1,counter-1
    xC = sqrt(AllData(i,1)**2 + a**2) * sin(AllData(i,2))*cos(AllData(i,3))
    yC = sqrt(AllData(i,1)**2 + a**2) * sin(AllData(i,2))*sin(AllData(i,3))
    zC = AllData(i,1) * cos(AllData(i,2))
    
 !   print *, xC, yC, zC

    write(10,*) xC,yC,zC

enddo


close(10)



endif









enddo
!!!!!!$OMP END PARALLEL DO








stop
print *, 'Ray Tracing completed succesfully'
call eval_performance(v,dk)
print *, 'The relative error in kappa was:', dk
print *, 'The total accumulated relative error in the variables was:'
print *, delta(1:3) 
print *, delta(4:6) 




end program main
