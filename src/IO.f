module IO

use parameters
use constants



implicit none

private

public FileOpen, ToTextFile, create_RT_plotfile, convert_from_binary

contains






subroutine convert_from_binary(f,nrows,ncols)
!Argument
character(len=200) :: f !This is the binary data file path
integer(kind=dp) :: ncols,nrows
!Other
real(kind=dp), dimension(:,:), allocatable :: DataArray
integer(kind=dp) :: stat


print *, f
!Open the .dat file
open(unit=10, file=f , form='unformatted')

!Allocate the data array
allocate(DataArray(nrows, ncols))



print *, shape(DataArray)
!Read in data
read(10,iostat=stat) DataArray

print *, DataArray(2,:)




end subroutine convert_from_binary






subroutine ToTextFile()

real(kind=dp), dimension(13) :: row
real(kind=dp), dimension(:,:), allocatable :: array
integer(kind=dp) :: i,stat


allocate(array(nrows,4))

!Open the .dat file
open(unit=10, file=MPDBinaryData , form='unformatted',access='stream')

!Open the .txt file
open(unit=20, file=MPDFormatData , form='formatted',access='stream')




stat = 0

do while (stat .EQ. 0)

    read(10,iostat=stat) array
    
    print *, array(1,:)
    stop





    do i=1,nrows
    
    if (array(i,2) .EQ. 0.0_dp) then
    exit
    endif

    
    write(20,*) array(i,2), array(i,3), array(i,4), a
    !print *,  array(i,2), array(i,3), array(i,4), a
    !print *, stat
    
    
    enddo
    array = 0.0_dp



enddo


close(10)
close(20)





end subroutine ToTextFile


subroutine FileOpen(f)
character(len=200) :: f
logical :: res


inquire( file=f, exist=res )
if (res) then

open(unit=10, file=MPDBinaryData,status='old',position='append', form='unformatted',access='stream')
        

else

open(unit=10, file=MPDBinaryData,status='replace', form='unformatted',access='stream')


endif


end subroutine FileOpen


subroutine create_RT_plotfile(alpha,beta,RTfile)
!Arguments
real(kind=dp),intent(in) :: alpha,beta
character(len=300), intent(out) :: RTFile
!Other
character(len=300) :: aSTR, bSTR
character(len = 20), parameter :: FMT1 = "(F10.2)"


write(aSTR, FMT1) alpha
write(bSTR, FMT1) beta


RTFile = trim(adjustl(RTPath))//'RTFile_alpha='//trim(adjustl(aSTR))//'_beta='//trim(adjustl(bSTR))//'.txt'
RTFile = trim(adjustl(RTFile))



end subroutine create_RT_plotfile

end module IO
