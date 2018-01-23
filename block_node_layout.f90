!==================================================================================================!
!  ***  *  *  ***   **** ***** !                                                                   !
! *     *  * *   * *       *   !                                                                   !
! *  *  **** *   *  ***    *   !                                                                   !
! *   * *  * *   *     *   *   !                                                                   !
!  ***  *  *  ***  ****    *   !                                                       C. Thieulot !
!==================================================================================================!

subroutine block_node_layout

use structures

implicit none

integer counter,nelx,nely,i,j,ib,nplevel,il,offset

!==================================================================================================!
!@@ \subsection{\tt block\_node\_layout}
!@@ This subroutine generates the position of the mesh node points for each 
!@@ block as well as their connectivity (i.e. the list of nodes making 
!@@ each element. 
!@@ \begin{center}
!@@ \includegraphics[width=8cm]{images/basics}\\
!@@ {\small a) block layout for mtype=3,4; b) block layout for mtype=5}
!@@ \end{center}
!@@ In the case that the block is a quadrilateral, it is assumed that 
!@@ it is made of $l \times l$ quadrilateral cells. In the case that the block is a 
!@@ triangle, it is made of $l^2$ triangular cells.
!==================================================================================================!

call cpu_time(t3)

do ib=1,nblock

select case(mtype)

case(06,12) ! cubed sphere, citcom mesh

   nelx=level
   nely=level
   Lx=1.d0
   Ly=1.d0
   counter=0
   do j=0,nely
   do i=0,nelx
      counter=counter+1
      block(ib)%x(counter)=dble(i)*Lx/dble(nelx)
      block(ib)%y(counter)=dble(j)*Ly/dble(nely)
      block(ib)%z(counter)=0.d0
   end do
   end do


   counter=0
   do j=1,nely
   do i=1,nelx
      counter=counter+1
      block(ib)%icon(1,counter)=i+(j-1)*(nelx+1)
      block(ib)%icon(2,counter)=i+1+(j-1)*(nelx+1)
      block(ib)%icon(3,counter)=i+1+j*(nelx+1)
      block(ib)%icon(4,counter)=i+j*(nelx+1)
   end do
   end do

case(20) ! icosahedron

   nplevel=level+1
   counter=0
   do il=1,level+1
      do i=1,nplevel
         counter=counter+1
         block(ib)%x(counter)=dble(i-1)/level
         block(ib)%y(counter)=dble(il-1)/level
         block(ib)%z(counter)=0.d0
      end do
      nplevel=nplevel-1
   end do

   counter=0
   offset=0
   nplevel=level
   do il=1,level
      do i=1,nplevel-1
         counter=counter+1
         block(ib)%icon(1,counter)= i            +offset
         block(ib)%icon(2,counter)= i+1          +offset
         block(ib)%icon(3,counter)= i+(nplevel+1)+offset
         counter=counter+1
         block(ib)%icon(1,counter)= i+1            +offset
         block(ib)%icon(2,counter)= i+(nplevel+1)+1+offset
         block(ib)%icon(3,counter)= i+(nplevel+1)  +offset
      end do
      counter=counter+1
      block(ib)%icon(1,counter)= nplevel            +offset
      block(ib)%icon(2,counter)= nplevel+1          +offset
      block(ib)%icon(3,counter)= nplevel+(nplevel+1)+offset
      offset=offset+nplevel+1
      nplevel=nplevel-1
   end do


case default

   stop 'block_node_layout: wrong value for mtype'

end select

end do

call cpu_time(t4) ; write(*,'(a,f6.3,a)') 'block_node_layout:',t4-t3,'s'

end subroutine

!==================================================================================================!
!==================================================================================================!
