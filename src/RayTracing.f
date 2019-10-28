module RayTracing


use parameters
use constants
use metric
use NumericalMethods
use derivatives
use IO

implicit none

public RT

private RT_Forward, RT_Backward, GeneralInitialConditions, shoot, aim

contains

subroutine RT()
integer(kind=dp) :: stat,i

real(kind=dp), dimension(:,:), allocatable :: MPDData
real(kind=dp) :: t,r,theta,phi,tau
real(kind=dp), dimension(4) :: xvector, uvector, svector
real(kind=dp) :: p1,p2,p3,p4
real(kind=dp) :: mag, rstart
real(kind=dp), dimension(6) :: v
real(kind=dp), dimension(4) :: c
real(kind=dp),dimension(:,:),allocatable :: PlotData !Big array to save all data
integer(kind=dp) :: counter , jj
real(kind=dp) :: xC, yC, zC,cpt
real(kind=dp),dimension(4,4) :: metric_covar
real(kind=dp) :: gtt, gpp, gtp, aa,bb,cc, Eobs
character(len=200) :: FileName,IDStr
integer(kind=dp) :: FileID




allocate(MPDData(OrbitNrows,OrbitNcols))

!Load the MPD Data file
open(unit=10, file=MPDBinaryData , form='unformatted',action='read')
read(10,iostat=stat) MPDData
close(10)


!Find the intersecting ray

!For each row in array
do i = 1,OrbitNrows

    !Get the data as a position vector t, r, theta, phi
    xvector(1:4) = MPDData(i,1:4)
 
!xvector(1) = 5.00000000000000000000000000000000000         
!xvector(2) = 1484.55810146794412658580140724669970         
!xvector(3) = 1.57079632679489661923132169163975140        
!xvector(4) = 9.42477796076937971538793014983850839




xvector(1) =  1.00000000000000000000000000000000000         
xvector(2) = 394.416820652839295349451495054830836         
xvector(3)=1.57079632679489661923132169163975140         
xvector(4)=1.88495559215387594307758602996770164      


    print *, 'XVector:', xvector

    !Do either forward or backward Ray Tracing

    if (RayTracingDirection .EQ. +1.0_dp) then
        !Forward
        call RT_Forward()
    else if (RayTracingDirection .EQ. -1.0_dp) then
        !Backward
        call RT_Backward(xvector)
    else
        !Exit
        print *, 'Error: You need to set the ray tracing method'

    endif



enddo





end subroutine RT



subroutine RT_Forward()

        print *, 'Forward Ray Tracing. Doing Nothing'
end subroutine RT_Forward

subroutine RT_Backward(xi)
!Arguments
real(kind=dp), dimension(4), intent(in) :: xi !position of MPD
!Other
real(kind=dp) :: xT, yT, zT, mm
real(kind=dp) :: alpha,beta
real(kind=dp),dimension(6) :: IO
real(kind=dp),dimension(4) :: globals
real(kind=dp),dimension(3) :: OT
integer(kind=dp) :: plot !Do you want to save ray path?
integer(kind=dp) :: RayClass

!The target points 
mm = sqrt(xi(2)**2 + a**2)
xT = mm * sin(xi(3))*cos(xi(4))
yT = mm * sin(xi(3))*sin(xi(4))
zT = mm * cos(xi(3))

IO(1) = xT
IO(2) = yT
IO(3) = zT

RayClass = 1 !Minkowski Initialisation
RayClass = 3
!TO DO()
!Secondary rays
!Speed optimization
!Robustness for edge cases e.g. hits BH

!Initial guess at the coordinates of the image plane
!alpha = -10.0_dp*yT !Falls below horizon for alpha=-yT (=-5.0_dp)

!Intiailly ray class = minkowski


11 call setAB(IO,OT,RayClass)

print *, 'IC:', RayClass, IO(4), IO(5)
!Fire a ray. How it does i.e. ds is written to IO(6)
plot = 0

call shoot(IO,plot)
print *, 'Initial ds =', IO(6)

if (IO(6) .EQ. -1.0_dp) then
!Fell into BH. Poor choice of initial ray. Search for more bent rays instead
RayClass = RayClass + 1
goto 11
endif




!Now run some conjugate gradient descent (non-linear) using the AIM subroutine
OT(1) = 0.010_dp ! Initial stepsize alpha for backtracing
globals = 0.0_dp

do while (IO(6) .GT. 1.0d-6)
call aim(IO,globals,OT)
enddo


!Trace the ray - useful for data+plotting. Or just write to file for use later

!The mininum has been sucessfully found. Now tracing the ray
print *, 'The mininum has been sucessfully found. Now tracing the ray'
plot=1
call shoot(IO,plot)








if (xT .LT. 0.0_dp .and. RayClass .NE. 3) then
RayClass = 3 
goto 11
endif



end subroutine RT_Backward


subroutine aim(IO,globals,OT)
!Arguments
real(kind=dp),dimension(6),intent(inout) :: IO
real(kind=dp),dimension(4),intent(inout) :: globals
real(kind=dp),dimension(:),intent(inout) :: OT !i.e. other
!Other
real(kind=dp),parameter :: gbit = 1.0e-18_dp
real(kind=dp) :: dsA, dsB, gA, gB, dsORIG, zeta,hA, hB
real(kind=dp) :: etaRT, factor, alpha0, beta0,dsBEST,aBEST,bBEST
real(kind=dp) :: tau, c,t, pA, pB,norm,m, dstemp !Armijo
integer(kind=dp), parameter :: plot = 0 !never want to plot when  aiming.





!Load original ds
alpha0 = IO(4)
beta0 = IO(5)
dsORIG = IO(6)




!Alpha gradient gA
IO(4) = IO(4) + gbit
call shoot(IO,plot)
dsA = IO(6)
gA = (dsA - dsORIG)/gbit
gA = - gA 
IO(4) = IO(4) - gbit



!Beta gradient gB
IO(5) = IO(5) + gbit
call shoot(IO,plot)
dsB = IO(6)
gB = (dsB - dsORIG)/gbit
gB = -gB 
IO(5) = IO(5) - gbit







!Now get the alpha/beta adjustment direction

if (globals(1) .EQ. 0.0_dp) then
!It is the first time, just set zeta = 0
zeta = 0.0_dp
else
zeta = (gA*gA +gB*gB)/(globals(1)**2 + globals(2)**2)
endif


hA = gA +zeta*globals(3) 
hB = gB +zeta*globals(4)


!Got the direction. Now perform a line search

!set inital optimization params
etaRT = OT(1)
factor = 2.0_dp

!loop line search

dsBEST = 1d20
aBEST = alpha0
bBEST = beta0






tau = 0.10_dp
c = 0.0010_dp

!Testing backtrack line search

norm = sqrt(hA**2 + hB**2)
pA = hA/norm
pB = hB/norm

m = gA*pA + gB*pB !Is this the correct defenition?
t = -c*m

etaRT = OT(1)
dstemp = IO(6)


!print *, alpha0, etaRT*hA
!stop

02 do

!if success on first go, try a larger step size?

        IO(4) = alpha0 + etaRT*hA
        IO(5) = beta0 + etaRT*hB
 
        !Check sign of alpha doesnt change
        if ( sign(1.0_dp,IO(4)/alpha0) .NE. 1.0_dp) then
        globals = 0
        etaRT = etaRT*0.10_dp
        goto 02
        endif





        call shoot(IO,plot)

        print *, 'i:', IO(4:6)



       ! if ((IO(6) - dstemp) .LT. etaRT*t) then
        
        if (IO(6) .LT. dstemp) then


        !Found something good. EXIT



!Update globals
globals(1) = gA
globals(2) = gB
globals(3) = hA
globals(4) = hB

OT(1) = etaRT*1.10_dp !Try small increase in stepsize for next time?

if (etaRT .LT. 1e-10_dp) then
!etaRT is very small.
!Reset to normal
OT(1) =1.0_Dp
else
OT(1) = etaRT  !etaRT alwas decas Is this good?
endif




print *, '------'
        return



        else


        etaRT = etaRT*tau

        if (etaRT .LT. 1e-10_dp) then
        !Reset
        IO(6) = dsORIG
        IO(5) = beta0
        IO(4) = alpha0
        globals = 0.0_dp
        OT(1) = 1.0_dp
        return
        endif

        endif





enddo






stop


end subroutine aim




subroutine shoot(IO,plot)
!Arguments
real(kind=dp),dimension(6),intent(inout) :: IO
integer(kind=dp), intent(in) :: plot
!Other
real(kind=dp) :: xT, yT, zT, xP, yP, zP,mm
real(kind=dp) :: alpha,beta
real(kind=dp) :: xprime, yprime, zprime
real(kind=dp) :: w, r, theta,phi,t
real(kind=dp) :: rdot, thetadot, phidot, zdot, sig,u,vv
real(kind=dp), dimension(7) :: ray !Ray initial conditions
real(kind=dp), dimension(6) :: v !The variables e.g. r, theta phi etc
real(kind=dp), dimension(4) :: c !The constants L. kappa etc
real(kind=dp), dimension(2) :: b !the BackArray - contains stuff which relates to intersections and backwards ray tracing
integer(kind=dp) :: counter,i
real(kind=dp) :: Rstart, xOUT,dx, ds2
integer(kind=dp), parameter :: nrows = 1e8,ncols = 3
real(kind=dp), dimension(:,:), allocatable :: PlotArray
character(len=300) :: RTFile


if (plot .EQ. 1) then
!Set up to plot
allocate(PlotArray(nrows,ncols))
call create_RT_plotfile(IO(4),IO(5),RTFile)
print *, 'create RT plotfile:', RTFile
endif

!Read in array
xT = IO(1)
yT = IO(2)
zT = IO(3)
alpha = IO(4)
beta = IO(5)




!Convert to the primed Cartesian frame
xprime = sqrt(Robs**2.0_dp +a**2.0_dp) * sin(ThetaObs) - beta*cos(ThetaObs)
yprime = alpha
zprime = Robs*cos(ThetaObs) + beta*sin(ThetaObs)




!Convert it to Boyer Lindquist
w = xprime**2.0_dp +yprime**2.0_dp +zprime **2.0_dp - a**2.0_dp
r = sqrt((w+sqrt(w**2.0_dp+4.0_dp*a**2.0_dp*zprime**2.0_dp))/2.0_dp)
theta = acos(zprime/r)
phi = atan2(yprime,xprime)
t=0.0_dp


!And then get the derivatives
sig = r**2.0_dp +(a*cos(theta))**2.0_dp
u = sqrt(r**2.0_dp+a**2.0_dp)
vv= -sin(ThetaObs)*cos(phi)
zdot = -1.0_dp



rdot = -zdot*(-u**2.0*cos(ThetaObs)*cos(theta)+r*u*vv*sin(theta))/sig
thetadot = -zdot*(cos(ThetaObs)*r*sin(theta)+u*vv*cos(theta))/sig	
phidot = -zdot*sin(ThetaObs)*sin(phi)/(u*sin(theta))



!Write to array
ray(1) = t
ray(2) = r
ray(3) = theta
ray(4) = phi
ray(5) = rdot
ray(6) = thetadot
ray(7) = phidot

!Set up initial conditions
call GeneralInitialConditions(ray,v,c)
Rstart = v(1)
counter = 1

b(1) = xT
b(2) = 0.0_dp







!Iterate
do while (v(1) .GT. Rhor .and. v(1) .LT. rstart*10.0_dp)
!Do a timestep

call rk_geodesic(v,c,b)
!print *, plot, v(1)
if (plot .EQ. 1) then
!Save to array for plotting
PlotArray(counter,1:3) = v(1:3) 
endif
counter = counter + 1


!Check to see if condition is satisfied
if (b(2) .EQ. 1.0_dp) then

!Calculate ds2
mm = sqrt(v(1)**2 + a**2)
xP = mm*sin(v(2))*cos(v(3))
yP = mm*sin(v(2))*sin(v(3))
zP = v(1)*cos(v(2))

ds2 = (xP - xT)**2 + (yP - yT)**2 + (zP-zT)**2
IO(6) = ds2



!Save the plot array


if (plot .EQ. 1) then

!MAKE THIS A SUBROUTINE

open(unit =10, file = RTFile,form='formatted', status='replace')
do i=1,counter-1
mm = sqrt(PlotArray(i,1)**2 + a**2)
xP = mm*sin(PlotArray(i,2))*cos(PlotArray(i,3))
yP = mm*sin(PlotArray(i,2))*sin(PlotArray(i,3))
zP = PlotArray(i,1)*cos(PlotArray(i,2))
write(10,*) xP , yP, zP
enddo
close(10)
print *, 'saved'
deallocate(PlotArray)
endif





return
endif

enddo


!If it gets this far it means that it fell below the horizon
!Return error code
IO(6) = -1.0_dp
print *, 'Fell below horizon'






!open(unit =10, file = RTFile,form='formatted')
!do i=1,counter-1
!mm = sqrt(PlotArray(i,1)**2 + a**2)
!xP = mm*sin(PlotArray(i,2))*cos(PlotArray(i,3))
!yP = mm*sin(PlotArray(i,2))*sin(PlotArray(i,3))
!zP = PlotArray(i,1)*cos(PlotArray(i,2))
!write(10,*) xP , yP, zP
!enddo
!close(10)
!print *, 'saved'
!deallocate(PlotArray)

end subroutine shoot





subroutine GeneralInitialConditions(ray,v,c)
!Argument
real(kind=dp), dimension(7) :: ray !Ray initial conditions
real(kind=dp), dimension(6), intent(out) :: v
real(kind=dp), dimension(4), intent(out) :: c
!Other
real(kind=dp) :: r,theta,rdot,thetadot,phidot, sigma,delta
real(kind=dp) :: E2, E, Enorm, Eprime2, Eprime,Eobs
real(kind=dp) :: B2, fr,ft,omega2
real(kind=dp) :: Lz, pr, ptheta, phi,kappa, En, s1

!Load the data
r = ray(2)
theta = ray(3)
phi = ray(4)
rdot = ray(5)
thetadot = ray(6)
phidot = ray(7)

!Setup some defenitions
sigma = r**2 +a**2*cos(theta)**2
delta = r**2 -2.0_dp*r + a**2



!Define the plasma frequency
B2 =  N*4.0_dp*PI*electron_charge**2 / electron_mass
call plasma_fr(r,fr)
call plasma_ft(r,ft)
omega2 = B2 * (fr+ft)/sigma

!Construct the energy is the conserved energy?
!Eobs = 16.0_dp
!Enorm = (sigma-2.0_dp*r)*(rdot**2/delta + thetadot**2) + delta*(sin(theta)*phidot)**2
!Eprime2 = (Eobs**2)/Enorm
!Eprime = sqrt(Eprime2)

!Correct to ensure correct Energies
!rdot = rdot * Eprime
!thetadot = thetadot*Eprime
!phidot = phidot*Eprime


!If you want you can check that E = Eobs
!E2 = (sigma-2.0_dp*r)*(rdot**2/delta + thetadot**2) + delta*(sin(theta)*phidot)**2
!E = sqrt(E2)


! Compute the energy
s1 = sigma-2.0*r
E2 = s1*(rdot**2.0/delta +thetadot**2.0 +omega2/sigma) + delta*sin(theta)**2.0*phidot**2.0
En  = sqrt(E2)






!Get the angular momentum
Lz = (sigma*delta*phidot - 2.0_dp*a*r*En)*sin(theta)**2 / (s1)


!Get the momenta (non constant)
pr = rdot * sigma/delta
ptheta = sigma*thetadot

!Normalise to E = 1
pr = pr/En
ptheta = ptheta/En
Lz = Lz/En
B2 = B2 / E2




!Define one last constant and export
kappa = ptheta**2 + Lz**2/sin(theta)**2 + a**2*sin(theta)**2

v(1) = r
v(2) = theta
v(3) = phi
v(4) = 0.0_dp !t
v(5) = pr
v(6) = ptheta


!Dont declare these globally as need to be careful when running in parallel
c(1) = Lz
c(2) = kappa
c(3) = B2
c(4) = 1.0d-3 !Initial stepsize for RT


!print *, 'Iliad IC:'
!print *, v(1:3)
!print *, v(4:6)
!print *, c, En



end subroutine GeneralInitialConditions



subroutine setAB(IO,OT, RayClass)
!Arguments
real(kind=dp),dimension(:), intent(inout) :: IO
integer(kind=dp), intent(in) :: RayClass
real(kind=dp),dimension(3),intent(inout) :: OT
!Other
real(kind=dp) :: xT,yT, zT, alpha, beta


print *, 'RayClass = ', RayClass

xT = IO(1)
yT = IO(2)
zT = IO(3)

if (RayClass .EQ. 1) then

alpha = yT
beta = xT*cos(ThetaObs) + zT*sin(ThetaObs)
OT(1) = 1.0_dp
else if (RayClass .EQ. 2) then

alpha = sign(50.0_dp, yT)
beta = xT*cos(ThetaObs) + zT*sin(ThetaObs)
OT(1) = 0.010_dp
else if (RayClass .EQ. 3) then

alpha = -sign(50.0_dp, yT)
beta = xT*cos(ThetaObs) + zT*sin(ThetaObs)
OT(1) = 0.010_dp

endif

IO(4) = alpha
IO(5) = beta

end subroutine setAB


end module RayTracing
