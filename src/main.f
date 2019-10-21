program main


use parameters
use constants
use tensors
use initial_conditions
use integrate
use numerical_performance




implicit none 
real(kind=dp), dimension(4) :: ki !variables
real(kind=dp), dimension(3) :: xi !local vector location
real(kind=dp), dimension(3) :: xi_global !global vector location
real(kind=dp), dimension(6) :: v !Variables (r,theta,phi,t,pr,ptheta)
real(kind=dp), dimension(4) :: c
real(kind=dp) :: xglobal, yglobal, zglobal, w

real(kind=dp) :: rmag, mag_check1, mag_check2, dk

integer(kind=dp) :: counter, j,i, save_unit
real(kind=dp),dimension(:,:),allocatable :: AllData !Big array to save all data
real(kind=dp) :: mm, xC, yC, zC,r_start
real(kind=dp) :: dummy_angle,res, E, ww, AA
character(len=200) :: ID
!Set up save location
call get_environment_variable("DarkMatterDir", path)

!Set up initial location

if (PSR_Location .EQ. 1) then
ww = xCOM**2 + yCOM**2 + zCOM**2 -a**2
AA = sqrt(ww**2 + 4.0_dp*a**2*zCOM**2)

rCOM = sqrt(0.50_dp * (ww+AA))
thetaCOM = acos(zCOM/rCOM)
phiCOM = atan2(yCOM,xCOM)


elseif (PSR_Location .EQ. 0) then
rCOM = rC
thetaCOM = thetaC
phiCOM = phiC
endif



!Some info for the user
print *, 'The intial location of the ray is'




res = 10.0_dp
allocate(AllData(int(1e6),6))
!$OMP PARALLEL DO PRIVATE(ki, v, &
!$OMP& counter,AllData,xC,yC,zC, &
!$OMP& dummy_angle, ID,i,c,save_unit, r_start,E) 

do j = 1,int(res)


E = 0.10_dp + real(j,kind=dp) * (2.0_dp - 0.10_dp)/res


print *, 'f = ', E, 'GHz'

E = 2.0_dp*PI*E*1d9 !GHz
!dummy_angle = (PI/res) * real(j,kind=dp)
!dummy_angle = PI/2.0_dp




ki(1) = 1.0_dp !xdot
ki(2) = 0.0_dp !cos(dummy_angle) !ydot
ki(3) = 0.0_dp !zdot
ki(4) = E



call set_initial_conditions(ki,v,c)
r_start = v(1)

!print *, j, ki, v(5:6)
counter = 1



!do while (counter .LT. 500)
do while (v(1) .GT. Rhor .and. v(1) .LT. r_start*10.0_dp)



!print *, j
call rk(v,c)




if (plot .EQ. 1) then
!Save data at each integration point for plotting - probably not necessary for PSR timing (where we only care about endpoints)
AllData(counter,:) = v

!print *, j, v(1)

endif
counter = counter + 1

enddo










if (plot .EQ. 1) then


!write( ID, '(f10.2)' )  dummy_angle
write( ID, '(f10.2)' )  E/(2.0_dp*PI*1d9)



save_unit = 1000*(j+1)




open(unit = save_unit, file = trim(adjustl(path))//'RT_'//trim(adjustl(ID))//'.txt', status = 'replace', form='formatted')



do i=1,counter-1
    xC = sqrt(AllData(i,1)**2 + a**2) * sin(AllData(i,2))*cos(AllData(i,3))
    yC = sqrt(AllData(i,1)**2 + a**2) * sin(AllData(i,2))*sin(AllData(i,3))
    zC = AllData(i,1) * cos(AllData(i,2))
    
 !   print *, xC, yC, zC

    write(save_unit,*) xC,yC,zC

enddo


close(save_unit)



endif




enddo
!$OMP END PARALLEL DO



print *, 'Ray Tracing completed succesfully'




end program main
