!==================================================================================================!
!  ***  *  *  ***   **** ***** !                                                                   !
! *     *  * *   * *       *   !                                                                   !
! *  *  **** *   *  ***    *   !                                                                   !
! *   * *  * *   *     *   *   !                                                                   !
!  ***  *  *  ***  ****    *   !                                                       C. Thieulot !
!==================================================================================================!

subroutine output_blocks_ascii

use structures

implicit none

character(len=2) cib
integer i,ib

!==================================================================================================!
!==================================================================================================!

call cpu_time(t3)

if (generate_ascii_output) then

do ib=1,nblock
   call int_to_char(cib,2,ib)

   open(unit=123,file='OUT/ASCII/block_xyz'//cib//'.ascii',status='replace',form='formatted')
   write(123,*) block(ib)%ncell, block(ib)%np, block(ib)%nv
   do i=1,block(ib)%np
   write(123,'(3es20.10)') block(ib)%x(i),block(ib)%y(i),block(ib)%z(i)
   end do
   close(123)

   open(unit=123,file='OUT/ASCII/block_icon'//cib//'.ascii',status='replace',form='formatted')
   write(123,*) block(ib)%ncell, block(ib)%np, block(ib)%nv
   do i=1,block(ib)%ncell
   write(123,*) block(ib)%icon(:,i)
   end do
   close(123)

end do

end if

call cpu_time(t4) ; write(*,'(a,f6.3,a)') 'output_blocks_ascii:',t4-t3,'s'

end subroutine

!==================================================================================================!
!==================================================================================================!
