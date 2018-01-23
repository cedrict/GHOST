!==================================================================================================!
!  ***  *  *  ***   **** ***** !                                                                   !
! *     *  * *   * *       *   !                                                                   !
! *  *  **** *   *  ***    *   !                                                                   !
! *   * *  * *   *     *   *   !                                                                   !
!  ***  *  *  ***  ****    *   !                                                       C. Thieulot !
!==================================================================================================!

subroutine compute_normals_hollow_sphere

use structures

implicit none

integer i,inode,ic
real(8) :: bx,by,bz,cx,cy,cz 
real(8) :: ax,ay,az
real(8) :: x1,x2,x3,x4,x5,x6,x7,x8
real(8) :: y1,y2,y3,y4,y5,y6,y7,y8
real(8) :: z1,z2,z3,z4,z5,z6,z7,z8

!==================================================================================================!
!@@ \subsection{\tt compute\_normals\_hollow\_sphere}
!==================================================================================================!

call cpu_time(t3)

do i=1,hollow_sphere%np
   hollow_sphere%inner_node(i)=(abs(hollow_sphere%r(i)-R1)<1.d-6*R1) 
   hollow_sphere%outer_node(i)=(abs(hollow_sphere%r(i)-R2)<1.d-6*R2) 
end do

write(iunit,*) 'nb of nodes out inner boundary',count(hollow_sphere%inner_node)
write(iunit,*) 'nb of nodes out outer boundary',count(hollow_sphere%outer_node)

!===================================

hollow_sphere%inner_cell=.false.
hollow_sphere%outer_cell=.false.

do ic=1,hollow_sphere%ncell
   do i=1,hollow_sphere%nv
      inode=hollow_sphere%icon(i,ic)
      if (hollow_sphere%inner_node(inode)) hollow_sphere%inner_cell(ic)=.true.
      if (hollow_sphere%outer_node(inode)) hollow_sphere%outer_cell(ic)=.true.
   end do
end do

write(iunit,*) 'nb of cells out inner boundary',count(hollow_sphere%inner_cell)
write(iunit,*) 'nb of cells out outer boundary',count(hollow_sphere%outer_cell)

!===================================

select case(mtype)

case(06,12) ! quadrilaterals 

   do ic=1,hollow_sphere%ncell

      if (hollow_sphere%inner_cell(ic)) then

         x1=hollow_sphere%x(hollow_sphere%icon(1,ic)) 
         y1=hollow_sphere%y(hollow_sphere%icon(1,ic))  
         z1=hollow_sphere%z(hollow_sphere%icon(1,ic))  
         x2=hollow_sphere%x(hollow_sphere%icon(2,ic)) 
         y2=hollow_sphere%y(hollow_sphere%icon(2,ic))  
         z2=hollow_sphere%z(hollow_sphere%icon(2,ic))  
         x3=hollow_sphere%x(hollow_sphere%icon(3,ic)) 
         y3=hollow_sphere%y(hollow_sphere%icon(3,ic))  
         z3=hollow_sphere%z(hollow_sphere%icon(3,ic))  
         x4=hollow_sphere%x(hollow_sphere%icon(4,ic)) 
         y4=hollow_sphere%y(hollow_sphere%icon(4,ic))  
         z4=hollow_sphere%z(hollow_sphere%icon(4,ic))  

         !node 1 
         bx=x4-x1 ; by=y4-y1 ; bz=z4-z1
         cx=x2-x1 ; cy=y2-y1 ; cz=z2-z1
         call cross_product(bx,by,bz,cx,cy,cz,ax,ay,az)
         call normalise(ax,ay,az) 
         inode=hollow_sphere%icon(1,ic)
         hollow_sphere%nx(inode)=hollow_sphere%nx(inode)+ax
         hollow_sphere%ny(inode)=hollow_sphere%ny(inode)+ay
         hollow_sphere%nz(inode)=hollow_sphere%nz(inode)+az

         !node 2 
         bx=x1-x2 ; by=y1-y2 ; bz=z1-z2
         cx=x3-x2 ; cy=y3-y2 ; cz=z3-z2
         call cross_product(bx,by,bz,cx,cy,cz,ax,ay,az)
         call normalise(ax,ay,az) 
         inode=hollow_sphere%icon(2,ic)
         hollow_sphere%nx(inode)=hollow_sphere%nx(inode)+ax
         hollow_sphere%ny(inode)=hollow_sphere%ny(inode)+ay
         hollow_sphere%nz(inode)=hollow_sphere%nz(inode)+az

         !node 3 
         bx=x2-x3 ; by=y2-y3 ; bz=z2-z3
         cx=x4-x3 ; cy=y4-y3 ; cz=z4-z3
         call cross_product(bx,by,bz,cx,cy,cz,ax,ay,az)
         call normalise(ax,ay,az) 
         inode=hollow_sphere%icon(3,ic)
         hollow_sphere%nx(inode)=hollow_sphere%nx(inode)+ax
         hollow_sphere%ny(inode)=hollow_sphere%ny(inode)+ay
         hollow_sphere%nz(inode)=hollow_sphere%nz(inode)+az

         !node 4
         bx=x3-x4 ; by=y3-y4 ; bz=z3-z4
         cx=x1-x4 ; cy=y1-y4 ; cz=z1-z4
         call cross_product(bx,by,bz,cx,cy,cz,ax,ay,az)
         call normalise(ax,ay,az) 
         inode=hollow_sphere%icon(4,ic)
         hollow_sphere%nx(inode)=hollow_sphere%nx(inode)+ax
         hollow_sphere%ny(inode)=hollow_sphere%ny(inode)+ay
         hollow_sphere%nz(inode)=hollow_sphere%nz(inode)+az

      end if

      if (hollow_sphere%outer_cell(ic)) then

         x5=hollow_sphere%x(hollow_sphere%icon(5,ic)) 
         y5=hollow_sphere%y(hollow_sphere%icon(5,ic))  
         z5=hollow_sphere%z(hollow_sphere%icon(5,ic))  
         x6=hollow_sphere%x(hollow_sphere%icon(6,ic)) 
         y6=hollow_sphere%y(hollow_sphere%icon(6,ic))  
         z6=hollow_sphere%z(hollow_sphere%icon(6,ic))  
         x7=hollow_sphere%x(hollow_sphere%icon(7,ic)) 
         y7=hollow_sphere%y(hollow_sphere%icon(7,ic))  
         z7=hollow_sphere%z(hollow_sphere%icon(7,ic))  
         x8=hollow_sphere%x(hollow_sphere%icon(8,ic)) 
         y8=hollow_sphere%y(hollow_sphere%icon(8,ic))  
         z8=hollow_sphere%z(hollow_sphere%icon(8,ic))  

         !node 5
         bx=x6-x5 ; by=y6-y5 ; bz=z6-z5
         cx=x8-x5 ; cy=y8-y5 ; cz=z8-z5
         call cross_product(bx,by,bz,cx,cy,cz,ax,ay,az)
         call normalise(ax,ay,az) 
         inode=hollow_sphere%icon(5,ic)
         hollow_sphere%nx(inode)=hollow_sphere%nx(inode)+ax
         hollow_sphere%ny(inode)=hollow_sphere%ny(inode)+ay
         hollow_sphere%nz(inode)=hollow_sphere%nz(inode)+az

         !node 6
         bx=x7-x6 ; by=y7-y6 ; bz=z7-z6
         cx=x5-x6 ; cy=y5-y6 ; cz=z5-z6
         call cross_product(bx,by,bz,cx,cy,cz,ax,ay,az)
         call normalise(ax,ay,az) 
         inode=hollow_sphere%icon(6,ic)
         hollow_sphere%nx(inode)=hollow_sphere%nx(inode)+ax
         hollow_sphere%ny(inode)=hollow_sphere%ny(inode)+ay
         hollow_sphere%nz(inode)=hollow_sphere%nz(inode)+az

         !node 7
         bx=x8-x7 ; by=y8-y7 ; bz=z8-z7
         cx=x6-x7 ; cy=y6-y7 ; cz=z6-z7
         call cross_product(bx,by,bz,cx,cy,cz,ax,ay,az)
         call normalise(ax,ay,az) 
         inode=hollow_sphere%icon(7,ic)
         hollow_sphere%nx(inode)=hollow_sphere%nx(inode)+ax
         hollow_sphere%ny(inode)=hollow_sphere%ny(inode)+ay
         hollow_sphere%nz(inode)=hollow_sphere%nz(inode)+az

         !node 8
         bx=x5-x8 ; by=y5-y8 ; bz=z5-z8
         cx=x7-x8 ; cy=y7-y8 ; cz=z7-z8
         call cross_product(bx,by,bz,cx,cy,cz,ax,ay,az)
         call normalise(ax,ay,az) 
         inode=hollow_sphere%icon(8,ic)
         hollow_sphere%nx(inode)=hollow_sphere%nx(inode)+ax
         hollow_sphere%ny(inode)=hollow_sphere%ny(inode)+ay
         hollow_sphere%nz(inode)=hollow_sphere%nz(inode)+az

      end if

   end do

case(20) ! triangles

   do ic=1,hollow_sphere%ncell

      if (hollow_sphere%inner_cell(ic)) then

         x1=hollow_sphere%x(hollow_sphere%icon(1,ic)) 
         y1=hollow_sphere%y(hollow_sphere%icon(1,ic))  
         z1=hollow_sphere%z(hollow_sphere%icon(1,ic))  
         x2=hollow_sphere%x(hollow_sphere%icon(2,ic)) 
         y2=hollow_sphere%y(hollow_sphere%icon(2,ic))  
         z2=hollow_sphere%z(hollow_sphere%icon(2,ic))  
         x3=hollow_sphere%x(hollow_sphere%icon(3,ic)) 
         y3=hollow_sphere%y(hollow_sphere%icon(3,ic))  
         z3=hollow_sphere%z(hollow_sphere%icon(3,ic))  

         !node 1 
         bx=x3-x1 ; by=y3-y1 ; bz=z3-z1
         cx=x2-x1 ; cy=y2-y1 ; cz=z2-z1
         call cross_product(bx,by,bz,cx,cy,cz,ax,ay,az)
         call normalise(ax,ay,az) 
         inode=hollow_sphere%icon(1,ic)
         hollow_sphere%nx(inode)=hollow_sphere%nx(inode)+ax
         hollow_sphere%ny(inode)=hollow_sphere%ny(inode)+ay
         hollow_sphere%nz(inode)=hollow_sphere%nz(inode)+az

         !node 2 
         bx=x1-x2 ; by=y1-y2 ; bz=z1-z2
         cx=x3-x2 ; cy=y3-y2 ; cz=z3-z2
         call cross_product(bx,by,bz,cx,cy,cz,ax,ay,az)
         call normalise(ax,ay,az) 
         inode=hollow_sphere%icon(2,ic)
         hollow_sphere%nx(inode)=hollow_sphere%nx(inode)+ax
         hollow_sphere%ny(inode)=hollow_sphere%ny(inode)+ay
         hollow_sphere%nz(inode)=hollow_sphere%nz(inode)+az

         !node 3 
         bx=x2-x3 ; by=y2-y3 ; bz=z2-z3
         cx=x1-x3 ; cy=y1-y3 ; cz=z1-z3
         call cross_product(bx,by,bz,cx,cy,cz,ax,ay,az)
         call normalise(ax,ay,az) 
         inode=hollow_sphere%icon(3,ic)
         hollow_sphere%nx(inode)=hollow_sphere%nx(inode)+ax
         hollow_sphere%ny(inode)=hollow_sphere%ny(inode)+ay
         hollow_sphere%nz(inode)=hollow_sphere%nz(inode)+az

      end if


      if (hollow_sphere%outer_cell(ic)) then

         x4=hollow_sphere%x(hollow_sphere%icon(4,ic)) 
         y4=hollow_sphere%y(hollow_sphere%icon(4,ic))  
         z4=hollow_sphere%z(hollow_sphere%icon(4,ic))  
         x5=hollow_sphere%x(hollow_sphere%icon(5,ic)) 
         y5=hollow_sphere%y(hollow_sphere%icon(5,ic))  
         z5=hollow_sphere%z(hollow_sphere%icon(5,ic))  
         x6=hollow_sphere%x(hollow_sphere%icon(6,ic)) 
         y6=hollow_sphere%y(hollow_sphere%icon(6,ic))  
         z6=hollow_sphere%z(hollow_sphere%icon(6,ic))  

         !node 4
         bx=x5-x4 ; by=y5-y4 ; bz=z5-z4
         cx=x6-x4 ; cy=y6-y4 ; cz=z6-z4
         call cross_product(bx,by,bz,cx,cy,cz,ax,ay,az)
         call normalise(ax,ay,az) 
         inode=hollow_sphere%icon(4,ic)
         hollow_sphere%nx(inode)=hollow_sphere%nx(inode)+ax
         hollow_sphere%ny(inode)=hollow_sphere%ny(inode)+ay
         hollow_sphere%nz(inode)=hollow_sphere%nz(inode)+az

         !node 5
         bx=x6-x5 ; by=y6-y5 ; bz=z6-z5
         cx=x4-x5 ; cy=y4-y5 ; cz=z4-z5
         call cross_product(bx,by,bz,cx,cy,cz,ax,ay,az)
         call normalise(ax,ay,az) 
         inode=hollow_sphere%icon(5,ic)
         hollow_sphere%nx(inode)=hollow_sphere%nx(inode)+ax
         hollow_sphere%ny(inode)=hollow_sphere%ny(inode)+ay
         hollow_sphere%nz(inode)=hollow_sphere%nz(inode)+az

         !node 6
         bx=x4-x6 ; by=y4-y6 ; bz=z4-z6
         cx=x5-x6 ; cy=y5-y6 ; cz=z5-z6
         call cross_product(bx,by,bz,cx,cy,cz,ax,ay,az)
         call normalise(ax,ay,az) 
         inode=hollow_sphere%icon(6,ic)
         hollow_sphere%nx(inode)=hollow_sphere%nx(inode)+ax
         hollow_sphere%ny(inode)=hollow_sphere%ny(inode)+ay
         hollow_sphere%nz(inode)=hollow_sphere%nz(inode)+az

      end if

   end do

end select


do i=1,hollow_sphere%np
   if (hollow_sphere%inner_node(i)) then
      call normalise(hollow_sphere%nx(i),hollow_sphere%ny(i),hollow_sphere%nz(i))
   end if
   if (hollow_sphere%outer_node(i)) then
      call normalise(hollow_sphere%nx(i),hollow_sphere%ny(i),hollow_sphere%nz(i))
   end if
end do


call cpu_time(t4) ; write(*,'(a,f6.3,a)') 'compute_normals_hollow_sphere:',t4-t3,'s'

end subroutine

!==================================================================================================!
!==============================================================================!

subroutine cross_product(bx,by,bz,cx,cy,cz,ax,ay,az)
implicit none
real(8), intent(in) :: bx,by,bz,cx,cy,cz 
real(8), intent(out) :: ax,ay,az

ax=by*cz-bz*cy
ay=bz*cx-bx*cz
az=bx*cy-by*cx

end subroutine

!==============================================================================!

subroutine normalise(ax,ay,az)
implicit none
real(8), intent(inout) :: ax,ay,az
real(8) norm

norm=sqrt(ax**2+ay**2+az**2)
ax=ax/norm
ay=ay/norm
az=az/norm

end subroutine

!==================================================================================================!
