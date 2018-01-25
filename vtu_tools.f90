!==================================================================================================!
!  ***  *  *  ***   **** ***** !                                                                   !
! *     *  * *   * *       *   !                                                                   !
! *  *  **** *   *  ***    *   !                                                                   !
! *   * *  * *   *     *   *   !                                                                   !
!  ***  *  *  ***  ****    *   !                                                       C. Thieulot !
!==================================================================================================!

!==============================================================================!
!@@ \subsection{\tt write\_positions}
!@@ This routine is part of vtu\_tools.f90
!==============================================================================!

subroutine write_positions(np,x,y,z,iunit)

implicit none

integer, intent(in) :: np,iunit
real(8), dimension(np), intent(in) :: x,y,z

integer ip


write(iunit,*) '<DataArray type="Float32" NumberOfComponents="3" Format="ascii">'
do ip=1,np
write(iunit,'(3es30.15)') x(ip),y(ip),z(ip)
end do
write(iunit,*) '</DataArray>'

end subroutine

!==============================================================================!
!@@ \subsection{\tt write\_field\_dp}
!@@ This routine is part of vtu\_tools.f90
!==============================================================================!

subroutine write_field_dp(np,array,fieldname,iunit)

implicit none

integer, intent(in) :: np,iunit
real(8), dimension(np), intent(in) :: array
character(*) fieldname

integer ip

write(iunit,*) '<DataArray type="Float32" Name="'//trim(fieldname)//'" Format="ascii">'
do ip=1,np
write(iunit,'(es35.12)') array(ip)
end do
write(iunit,*) '</DataArray>'

end subroutine

!==============================================================================!
!@@ \subsection{\tt write\_field\_int}
!@@ This routine is part of vtu\_tools.f90
!==============================================================================!

subroutine write_field_int(np,array,fieldname,iunit)

implicit none

integer, intent(in) :: np,iunit
integer, dimension(np), intent(in) :: array
character(*) fieldname

integer ip

write(iunit,*) '<DataArray type="Float32" Name="'//trim(fieldname)//'" Format="ascii">'
do ip=1,np
write(iunit,*) array(ip)
end do
write(iunit,*) '</DataArray>'

end subroutine

!==============================================================================!
!@@ \subsection{\tt write\_icon}
!@@ This routine is part of vtu\_tools.f90
!==============================================================================!

subroutine write_icon(nv,nel,icon,iunit)

implicit none

integer, intent(in) :: nv,nel,iunit
integer, dimension(nv,nel) :: icon

integer iel

write(iunit,*) '<DataArray type="Int32" Name="connectivity" Format="ascii">'
do iel=1,nel
write(iunit,*) icon(1:nv,iel)-1
end do
write(iunit,*) '</DataArray>'

end subroutine


!==============================================================================!
!@@ \subsection{\tt write\_offsets}
!@@ This routine is part of vtu\_tools.f90
!==============================================================================!

subroutine write_offsets(nv,nel,iunit)

implicit none

integer, intent(in) :: nel, iunit,nv
integer iel

write(iunit,*) '<DataArray type="Int32" Name="offsets" Format="ascii">'
write(iunit,*) (iel*nv,iel=1,nel)
write(iunit,*) '</DataArray>'

end subroutine

!==============================================================================!
!@@ \subsection{\tt write\_types}
!@@ This routine is part of vtu\_tools.f90
!==============================================================================!

subroutine write_types(nv,nel,iunit)

implicit none

integer, intent(in) :: nel, iunit,nv
integer iel,celltype

select case(nv)
case(3)
celltype=5
case(4)
celltype=9
case default
   stop 'write_types: discretisation is unknown'
end select

write(iunit,*) '<DataArray type="Int32" Name="types" Format="ascii">'
write(iunit,*) (celltype,iel=1,nel)
write(iunit,*) '</DataArray>'

end subroutine

!==================================================================================================!
!==================================================================================================!
