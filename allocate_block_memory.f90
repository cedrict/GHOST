!==================================================================================================!
!  ***  *  *  ***   **** ***** !                                                                   !
! *     *  * *   * *       *   !                                                                   !
! *  *  **** *   *  ***    *   !                                                                   !
! *   * *  * *   *     *   *   !                                                                   !
!  ***  *  *  ***  ****    *   !                                                       C. Thieulot !
!==================================================================================================!

subroutine allocate_block_memory

use structures

implicit none

integer ib

!==================================================================================================!
!@@ \subsection{\tt allocate\_block\_memory}
!@@ This subroutine assigns to every block {\sl i} a number of node 
!@@ points {\sl block(i)\%np} and a number of cells {\sl block(i)\%ncell}  
!@@ The number of vertices per cell is also stored in {\sl block(i)\%nv}.
!@@ It then loops over all blocks and allocates the required arrays 
!@@ to store the mesh nodes position in Cartesian and spherical coordinates
!@@ and the connectivity.  
!==================================================================================================!

allocate(block(nblock))

select case(mtype)

case(06,12) ! cubed sphere & citcom mesh

   block(1:nblock)%ncell=level**2
   block(1:nblock)%np=(level+1)**2
   block(1:nblock)%nv=4

case(20) ! icosahedron

   block(1:nblock)%ncell=level**2
   block(1:nblock)%np=(level+1)*(level+2)/2
   block(1:nblock)%nv=3

case default

   stop 'allocate_block_memory: wrong value for mtype'

end select

do ib=1,nblock
   allocate(block(ib)%x(block(ib)%np))
   allocate(block(ib)%y(block(ib)%np))
   allocate(block(ib)%z(block(ib)%np))
   allocate(block(ib)%icon(block(ib)%nv,block(ib)%ncell))
   allocate(block(ib)%hull(block(ib)%np))
   allocate(block(ib)%r(block(ib)%np))
   allocate(block(ib)%theta(block(ib)%np))
   allocate(block(ib)%phi(block(ib)%np))
end do

write(iunit,*) 'np per block   =',block(1)%np
write(iunit,*) 'ncell per block=',block(1)%ncell

end subroutine

!==================================================================================================!
!==================================================================================================!
