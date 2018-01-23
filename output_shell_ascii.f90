!==================================================================================================!
!  ***  *  *  ***   **** ***** !                                                                   !
! *     *  * *   * *       *   !                                                                   !
! *  *  **** *   *  ***    *   !                                                                   !
! *   * *  * *   *     *   *   !                                                                   !
!  ***  *  *  ***  ****    *   !                                                       C. Thieulot !
!==================================================================================================!

subroutine output_shell_ascii

use structures

implicit none

integer i

!==================================================================================================!
!@@ \subsection{\tt output\_shell\_ascii}
!==================================================================================================!

call cpu_time(t3)

if (generate_ascii_output)then

open(unit=123,file='OUT/ASCII/shell_xyz.vtu',status='replace',form='formatted')
write(123,*) shell%ncell,shell%np,shell%nv
do i=1,shell%np
write(123,'(3es15.5)') shell%x(i),shell%y(i),shell%z(i)
end do
close(123)

open(unit=123,file='OUT/ASCII/shell_icon.vtu',status='replace',form='formatted')
write(123,*) shell%ncell,shell%np,shell%nv
do i=1,shell%ncell
write(123,*) shell%icon(:,i)
end do
close(123)

end if

call cpu_time(t4) ; write(*,'(a,f6.3,a)') 'output_shell_ascii:',t4-t3,'s'

end subroutine

!==================================================================================================!
