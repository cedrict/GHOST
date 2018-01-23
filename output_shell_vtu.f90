!==================================================================================================!
!  ***  *  *  ***   **** ***** !                                                                   !
! *     *  * *   * *       *   !                                                                   !
! *  *  **** *   *  ***    *   !                                                                   !
! *   * *  * *   *     *   *   !                                                                   !
!  ***  *  *  ***  ****    *   !                                                       C. Thieulot !
!==================================================================================================!

subroutine output_shell_vtu

use structures

implicit none

integer, parameter :: uunit=123

!==================================================================================================!
!@@ \subsection{\tt output\_shell\_vtu}
!@@ If {\sl generate\_vtu\_output} is true 
!@@ a vtu file of the shell is created in {\sl OUT}.
!==================================================================================================!

call cpu_time(t3)

if (generate_vtu_output) then

open(unit=123,file='OUT/shell.vtu',status='replace',form='formatted')
write(123,*) '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="BigEndian">'
write(123,*) '<UnstructuredGrid>'
write(123,*) '<Piece NumberOfPoints="',shell%np,'" NumberOfCells="',shell%ncell,'">'
!.............................
write(123,*) '<Points>'
call write_positions(shell%np,shell%x,shell%y,shell%z,uunit)
write(123,*) '</Points>'
!.............................
write(123,*) '<PointData Scalars="scalars">'
call write_field_dp(shell%np,shell%r,'r',uunit)
call write_field_dp(shell%np,shell%theta,'theta',uunit)
call write_field_dp(shell%np,shell%phi,'phi',uunit)
write(123,*) '</PointData>'
!.............................
write(123,*) '<CellData Scalars="scalars">'
call write_field_dp(shell%ncell,shell%area,'area',uunit)
call write_field_int(shell%ncell,shell%blocknumber,'block number',uunit)
write(123,*) '</CellData>'
!.............................
write(123,*) '<Cells>'
call write_icon(shell%nv,shell%ncell,shell%icon,uunit)
call write_offsets(shell%nv,shell%ncell,uunit)
call write_types(shell%nv,shell%ncell,uunit)
write(123,*) '</Cells>'
write(123,*) '</Piece>'
write(123,*) '</UnstructuredGrid>'
write(123,*) '</VTKFile>'
close(123)

end if

call cpu_time(t4) ; write(*,'(a,f6.3,a)') 'output_shell_vtu:',t4-t3,'s'

end subroutine

!==================================================================================================!
!==================================================================================================!
