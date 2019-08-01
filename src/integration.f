module integrate



use parameters
use constants


implicit none

private geodesic, adaptive_shrink, adaptive_grow

public rk

contains


subroutine rk(v)
!Arguments
real(kind=dp), dimension(6), intent(inout) :: v !Variables (r,theta,phi,t,pr,ptheta)

!Other
real(kind=dp), dimension(6) :: dv
integer(kind=dp) :: i
real(kind=dp), dimension(6) :: ytemp1,ytemp2,ytemp3,ytemp4,ytemp5
real(kind=dp), dimension(6) :: yout, yerr,yscal,ratio
real(kind=dp) :: yi, hdv, errmax


!Step 0 --------------------------------------------
call geodesic(v,dv)
do i=1,6
hdv = dh*dv(i)
yi = v(i)
ytemp1(i) = yi + B21 * hdv
ytemp2(i) = yi + B31 * hdv
ytemp3(i) = yi + B41 * hdv
ytemp4(i) = yi + B51 * hdv
ytemp5(i) = yi + B61 * hdv


yout(i) = yi + c1 * hdv
yerr(i) = (c1 - cbar1 ) * hdv
yscal(i) = abs(yi) + abs(hdv) + 1.0d-3

enddo








!Step 1 --------------------------------------------
call geodesic(ytemp1,dv)
do i = 1,6
hdv = dh*dv(i)
ytemp2(i) = ytemp2(i) + B32*hdv
ytemp3(i) = ytemp3(i) + B42*hdv
ytemp4(i) = ytemp4(i) + B52*hdv
ytemp5(i) = ytemp5(i) + B62*hdv
!C and cbar = 0
enddo







!Step 2 --------------------------------------------

call geodesic(ytemp2,dv)
do i =1,6
hdv = dh*dv(i)
ytemp3(i) = ytemp3(i) + B43*hdv
ytemp4(i) = ytemp4(i) + B53*hdv
ytemp5(i) = ytemp5(i) + B63*hdv
yout(i) = yout(i) + c3*hdv
yerr(i) = yerr(i) + (c3 - cbar3 ) * hdv
enddo

!Step 3 --------------------------------------------
call geodesic(ytemp3,dv)
do i =1,6
hdv = dh*dv(i)
ytemp4(i) = ytemp4(i) + B54*hdv
ytemp5(i) = ytemp5(i) + B64*hdv

yout(i) = yout(i) + c4*hdv
yerr(i) = yerr(i) + (c4 - cbar4 ) * hdv

enddo


!Step 4 --------------------------------------------
call geodesic(ytemp4,dv)

do i =1,6
hdv = dh*dv(i)
ytemp5(i) = ytemp5(i) + B65*hdv
yerr(i) = yerr(i) + (-cbar5 ) * hdv
enddo


!Step 5 --------------------------------------------
call geodesic(ytemp5,dv)
do i = 1,6
hdv = dh*dv(i)
yout(i) = yout(i) + c6*hdv
yerr(i) = yerr(i) + (c6 - cbar6 ) * hdv
enddo




!Got all the info we need. Calculate the error
do i = 1,6
ratio(i) = abs( yerr(i) / yscal(i) )
enddo



errmax = max(ratio(1),ratio(2),ratio(3),ratio(4),ratio(5),ratio(6))
errmax = errmax*escal


if (errmax .GT. 1) then
!The error is too big. Reduce the step size and exit without updating the variable vector


call adaptive_shrink(errmax)


else

!The error is OK. Grow the stepsize a little and set the variables for the next integration step
call adaptive_grow(errmax)
v = yout

endif




end subroutine rk




subroutine geodesic(v,dv)
!Arguments
real(kind=dp), dimension(6) :: v,dv
!Other
real(kind=dp) :: r,theta,phi,t,pr,ptheta
real(kind=dp) ::sigma,delta, SD,csc
!Read in coordinate variables
r= v(1)
theta= v(2)
phi= v(3)
t= v(4)
pr= v(5)
ptheta= v(6)
!Define some useful quantities


sigma = r**2 + a**2 * cos(theta)**2
delta = r**2 -2.0_dp*r +a**2
csc = 1.0_dp /sin(theta)
SD = sigma*delta

dv(1) = pr*delta/sigma
dv(2) = ptheta/sigma
dv(3) = (2.0_dp*a*r + (sigma - 2.0_dp*r)*Lz*csc**2 ) / SD
dv(4) = 1.0_dp + (2.0_dp*r*(r**2+a**2)-2.0_dp*a*r*Lz)/SD

dv(5) = ( &
        -kappa*(r-1.0_dp) &
        +2.0_dp*r*(r**2+a**2) &
        -2.0_dp*a*Lz &
        )/SD &
        - 2.0_dp*pr**2*(r-1.0_dp)/sigma

dv(6) = sin(theta)*cos(theta)*(Lz**2 * csc**4 -a**2) / sigma

end subroutine geodesic


subroutine adaptive_shrink(errmax)
real(kind=dp) :: errmax, htemp


htemp = S*dh*(errmax**PSHRINK)

dh = sign(max(abs(htemp), 0.10_dp*abs(dh)),dh)

end subroutine adaptive_shrink

subroutine adaptive_grow(errmax)
real(kind=dp) errmax,ERRCON,hnext

ERRCON = (5.0_dp/S)**(1.0_dp/PGROW)

!print *, errmax ,ERRCON
if (errmax .GT. ERRCON) then
    hnext = S*dh*errmax**PGROW
 !   print *, 'errcon', S*errmax**PGROW
else
    hnext = 5.0_dp*dh
endif

dh = hnext
END SUBROUTINE adaptive_grow

end module integrate
