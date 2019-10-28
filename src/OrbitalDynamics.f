module OrbitalDynamics


use parameters
use constants
implicit none


public SimpleKeplerOrbit


private

contains


subroutine SimpleKeplerOrbit()
integer(kind=dp) , parameter :: res=10
integer(kind=dp) , parameter :: ncol=4
real(kind=dp) :: r, phi,dphi,phi0, theta
real(kind=dp), dimension(:,:), allocatable :: out_array
integer(kind=dp) :: i


allocate(out_array(res,ncol))


dphi = Norbit*2.0_dp*PI / real(res,kind=dp)
phi0 = 0.0_dp
theta = PI/2.0_dp

do i=0,res-1

phi = phi0 + real(i,kind=dp)*dphi
r = semi_latus / (1.0_dp + eccentricity*cos(phi))

out_array(i+1,1) = i !t
out_array(i+1,2) = r !r
out_array(i+1,3) = theta !theta
out_array(i+1,4) = phi !phi

enddo


!Write binary data to file
open(unit=10, file=MPDBinaryData,status='replace', form='unformatted')
write(10) out_array
close(10)

OrbitNrows = res
OrbitNcols = 4

!Write formatted data to file
!Used by tools/PlotOrbit.py

open(unit=20, file=MPDFormatData , form='formatted',access='stream')
do i=1,res
write(20,*) out_array(i,2), out_array(i,3), out_array(i,4), a
enddo
close(20)


!deallocate(out_array)
!allocate(out_array(res,ncol))

!open(unit=10, file=MPDBinaryData, form='unformatted',action='read')
!read(10) out_array
!close(10)


!print *, out_array(1,:)

!call ToTextFile()

!Convert binary to formatted for plotting
!call convert_from_binary(MPDBinaryData, res,ncol)

end subroutine SimpleKeplerOrbit


end module OrbitalDynamics
