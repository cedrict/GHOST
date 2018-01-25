!==================================================================================================!
!  ***  *  *  ***   **** ***** !                                                                   !
! *     *  * *   * *       *   !                                                                   !
! *  *  **** *   *  ***    *   !                                                                   !
! *   * *  * *   *     *   *   !                                                                   !
!  ***  *  *  ***  ****    *   !                                                       C. Thieulot !
!==================================================================================================!

subroutine output_hollow_sphere_vtu

use structures

implicit none

integer i,ic,inner,outer

!==================================================================================================!
!@@ \subsection{\tt output\_hollow\_sphere\_vtu}
!@@ If the {\sl generate\_vtu\_output} flag is true, 
!@@ it generates the {\sl hollow\_sphere.vtu} file in the {\sl OUT} folder. 
!==================================================================================================!

call cpu_time(t3)

if (generate_vtu_output)then

open(unit=123,file='OUT/hollow_sphere.vtu',status='replace',form='formatted')
write(123,*) '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="BigEndian">'
write(123,*) '<UnstructuredGrid>'
write(123,*) '<Piece NumberOfPoints="',hollow_sphere%np,'" NumberOfCells="',hollow_sphere%ncell,'">'
!.............................
write(123,*) '<Points>'
write(123,*) '<DataArray type="Float32" NumberOfComponents="3" Format="ascii">'
do i=1,hollow_sphere%np
write(123,'(3es13.5)') hollow_sphere%x(i),hollow_sphere%y(i),hollow_sphere%z(i)
end do
write(123,*) '</DataArray>'
write(123,*) '</Points>'
!.............................
write(123,*) '<PointData Scalars="scalars">'

if (s40rts) then
write(123,*) '<DataArray type="Float32" Name="dlnvs (S40RTS)" Format="ascii">'
do i=1,hollow_sphere%np
write(123,'(es13.5)') hollow_sphere%dlnvs(i)
end do
write(123,*) '</DataArray>'
write(123,*) '<DataArray type="Float32" Name="drho (S40RTS)" Format="ascii">'
do i=1,hollow_sphere%np
write(123,'(es13.5)') hollow_sphere%drho(i)
end do
write(123,*) '</DataArray>'
end if

write(123,*) '<DataArray type="Float32"  Name="r" Format="ascii">'
do i=1,hollow_sphere%np
write(123,'(es11.3)') hollow_sphere%r(i)
end do
write(123,*) '</DataArray>'
write(123,*) '<DataArray type="Float32"  Name="theta" Format="ascii">'
do i=1,hollow_sphere%np
write(123,'(es11.3)') hollow_sphere%theta(i)
end do
write(123,*) '</DataArray>'
write(123,*) '<DataArray type="Float32"  Name="phi" Format="ascii">'
do i=1,hollow_sphere%np
write(123,'(es11.3)') hollow_sphere%phi(i)
end do
write(123,*) '</DataArray>'
write(123,*) '<DataArray type="Float32"  Name="inner" Format="ascii">'
do i=1,hollow_sphere%np
inner=0 ; if (hollow_sphere%inner_node(i)) inner=1
write(123,'(i2)') inner
end do
write(123,*) '</DataArray>'
write(123,*) '<DataArray type="Float32"  Name="outer" Format="ascii">'
do i=1,hollow_sphere%np
outer=0 ; if (hollow_sphere%outer_node(i)) outer=1
write(123,'(i2)') outer
end do
write(123,*) '</DataArray>'
write(123,*) '<DataArray type="Float32"  Name="normal" NumberOfComponents="3" Format="ascii">'
do i=1,hollow_sphere%np
write(123,'(3es11.3)') hollow_sphere%nx(i),hollow_sphere%ny(i),hollow_sphere%nz(i)
end do
write(123,*) '</DataArray>'
write(123,*) '</PointData>'
!.............................
write(123,*) '<CellData Scalars="scalars">'
write(123,*) '<DataArray type="Float32"  Name="volume" Format="ascii">'
do ic=1,hollow_sphere%ncell
write(123,'(es11.3)') hollow_sphere%volume(ic)
end do
write(123,*) '</DataArray>'
write(123,*) '<DataArray type="Float32"  Name="blocknumber" Format="ascii">'
do ic=1,hollow_sphere%ncell
write(123,*) hollow_sphere%blocknumber(ic)
end do
write(123,*) '</DataArray>'
write(123,*) '<DataArray type="Float32"  Name="inner" Format="ascii">'
do ic=1,hollow_sphere%ncell
inner=0 ; if (hollow_sphere%inner_cell(ic)) inner=1
write(123,'(i2)') inner
end do
write(123,*) '</DataArray>'
write(123,*) '<DataArray type="Float32"  Name="outer" Format="ascii">'
do ic=1,hollow_sphere%ncell
outer=0 ; if (hollow_sphere%outer_cell(ic)) outer=1
write(123,'(i2)') outer
end do
write(123,*) '</DataArray>'
write(123,*) '</CellData>'
!.............................
write(123,*) '<Cells>'
write(123,*) '<DataArray type="Int32" Name="connectivity" Format="ascii">'
do ic=1,hollow_sphere%ncell
write(123,*) hollow_sphere%icon(1:hollow_sphere%nv,ic)-1
end do
write(123,*) '</DataArray>'

write(123,*) '<DataArray type="Int32" Name="offsets" Format="ascii">'
write(123,*) (ic*hollow_sphere%nv,ic=1,hollow_sphere%ncell)
write(123,*) '</DataArray>'

write(123,*) '<DataArray type="Int32" Name="types" Format="ascii">'
if (hollow_sphere%nv==6) write(123,*) (13,ic=1,hollow_sphere%ncell)
if (hollow_sphere%nv==8) write(123,'(i3)') (12,ic=1,hollow_sphere%ncell)
write(123,*) '</DataArray>'

write(123,*) '</Cells>'
write(123,*) '</Piece>'
write(123,*) '</UnstructuredGrid>'
write(123,*) '</VTKFile>'
close(123)

end if

call cpu_time(t4) ; write(*,'(a,f6.3,a)') 'output_hollow_sphere_vtu:',t4-t3,'s'

end subroutine

!==================================================================================================!
!==================================================================================================!
