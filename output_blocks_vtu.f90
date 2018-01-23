!==================================================================================================!
!  ***  *  *  ***   **** ***** !                                                                   !
! *     *  * *   * *       *   !                                                                   !
! *  *  **** *   *  ***    *   !                                                                   !
! *   * *  * *   *     *   *   !                                                                   !
!  ***  *  *  ***  ****    *   !                                                       C. Thieulot !
!==================================================================================================!

subroutine output_blocks_vtu

use structures

implicit none
integer ib
integer ic,i
character(len=2) cib

integer, parameter :: uunit=123

!==================================================================================================!
!@@ \subsection{\tt output\_blocks\_vtu}
!==================================================================================================!

call cpu_time(t3)

if (generate_vtu_output) then

do ib=1,nblock

   call int_to_char(cib,2,ib)

   open(unit=123,file='OUT/block'//cib//'.vtu',status='replace',form='formatted')
   write(123,*) '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="BigEndian">'
   write(123,*) '<UnstructuredGrid>'
   write(123,*) '<Piece NumberOfPoints="',block(ib)%np,'" NumberOfCells="',block(ib)%ncell,'">'
   !.............................
   write(123,*) '<Points>'
   call write_positions(block(ib)%np,block(ib)%x,block(ib)%y,block(ib)%z,uunit)
   write(123,*) '</Points>'
   !.............................
   write(123,*) '<PointData Scalars="scalars">'
   write(123,*) '<DataArray type="Float32"  Name="hull" Format="ascii">'
   do i=1,block(ib)%np
   if (block(ib)%hull(i)) then
   write(123,*) 1.
   else
   write(123,*) 0.
   end if
   end do
   write(123,*) '</DataArray>'
   write(123,*) '</PointData>'
   !.............................
   write(123,*) '<Cells>'
   write(123,*) '<DataArray type="Int32" Name="connectivity" Format="ascii">'
   do ic=1,block(ib)%ncell
   write(123,*) block(ib)%icon(1:block(ib)%nv,ic)-1
   end do
   write(123,*) '</DataArray>'
   write(123,*) '<DataArray type="Int32" Name="offsets" Format="ascii">'
   write(123,*) (ic*block(ib)%nv,ic=1,block(ib)%ncell)
   write(123,*) '</DataArray>'
   write(123,*) '<DataArray type="Int32" Name="types" Format="ascii">'
   if (block(ib)%nv==3) write(123,*) (5,ic=1,block(ib)%ncell)
   if (block(ib)%nv==4) write(123,*) (9,ic=1,block(ib)%ncell)
   if (block(ib)%nv==8) write(123,*) (12,ic=1,block(ib)%ncell)
   write(123,*) '</DataArray>'
   write(123,*) '</Cells>'
   write(123,*) '</Piece>'
   write(123,*) '</UnstructuredGrid>'
   write(123,*) '</VTKFile>'
   close(123)

end do

end if

call cpu_time(t4) ; write(*,'(a,f6.3,a)') 'output_blocks:',t4-t3,'s'

end subroutine

!==================================================================================================!
!==================================================================================================!
