!==================================================================================================!
!  ***  *  *  ***   **** ***** !                                                                   !
! *     *  * *   * *       *   !                                                                   !
! *  *  **** *   *  ***    *   !                                                                   !
! *   * *  * *   *     *   *   !                                                                   !
!  ***  *  *  ***  ****    *   !                                                       C. Thieulot !
!==================================================================================================!

subroutine merge_blocks

use structures

implicit none

integer ib,ip,nnp,jp,counter,nncell,ic,i
integer(4), dimension(:),allocatable :: pointto
integer, dimension(:), allocatable :: compact
real(8),dimension(:),allocatable :: tempx,tempy,tempz
real(8) distance,gxip,gyip,gzip
logical,dimension(:),allocatable :: doubble
logical(1),dimension(:),allocatable :: sides

!==================================================================================================!
!@@ \subsection{\tt merge\_blocks}
!@@ This subroutine loops over all the (shell) blocks and merges them together to 
!@@ arrive at a full shell. In order to do so it needs to remove point duplicates 
!@@ and also merges the connectivity arrays of all the blocks.
!@@ All shell blocks (x,y,z coordinates and hull information) are first merged 
!@@ into temporary arrays (tempx, tempy,tempz,sides) of size nblock*nnp.
!@@ A logical array doubble of the same size is then set to true when two points 
!@@ of two different blocks share the same position. At the same time 
!@@ another integer array pointto is used: if node m is colocated with point n then 
!@@ pointto(m)=n (for n$<$m).
!@@ The 'real' total number of nodes of the assembled shell can then be computed while
!@@ the total number of cells simply is the number of blocks multiplied by the 
!@@ number of cells each block contains.
!@@ All arrays pertaining to this assembled shell can then be allocated. 
!@@ Finally the doubble and pointto arrays are used to construct a single set of coordinates and 
!@@ a connectivity array for an assembled shell without duplicates. 
!==================================================================================================!

call cpu_time(t3)

if (nblock/=1) then

nnp=block(1)%np
nncell=block(1)%ncell

allocate(tempx(nblock*nnp))
allocate(tempy(nblock*nnp))
allocate(tempz(nblock*nnp))
allocate(sides(nblock*nnp))

do ib=1,nblock
  tempx((ib-1)*nnp+1:(ib-1)*nnp+nnp)=block(ib)%x
  tempy((ib-1)*nnp+1:(ib-1)*nnp+nnp)=block(ib)%y
  tempz((ib-1)*nnp+1:(ib-1)*nnp+nnp)=block(ib)%z
  sides((ib-1)*nnp+1:(ib-1)*nnp+nnp)=block(ib)%hull
end do

allocate(doubble(nblock*nnp)) ; doubble=.false.
allocate(pointto(nblock*nnp))

do ip=1,nblock*nnp
   pointto(ip)=ip
end do

distance=1.d-6*R2

counter=0
do ip=2,nblock*nnp
   if (sides(ip)) then 
   gxip=tempx(ip)
   gyip=tempy(ip)
   gzip=tempz(ip)
   do jp=1,ip-1
      if (sides(jp)) then 
      if (abs(gxip-tempx(jp))<distance .and. &
          abs(gyip-tempy(jp))<distance .and. &
          abs(gzip-tempz(jp))<distance       ) then 
          doubble(ip)=.true.
          pointto(ip)=jp
          exit 
      end if
      end if
   end do
   end if
end do

shell%np=nblock*nnp-count(doubble)

shell%ncell=nblock*block(1)%ncell

write(iunit,*) 'count(doubble)=',count(doubble)
write(iunit,*) 'shell%np',shell%np
write(iunit,*) 'shell%ncell',shell%ncell

shell%nv=block(1)%nv

allocate(shell%x(shell%np))
allocate(shell%y(shell%np))
allocate(shell%z(shell%np))
allocate(shell%icon(shell%nv,shell%ncell))
allocate(shell%r(shell%np))
allocate(shell%theta(shell%np))
allocate(shell%phi(shell%np))
allocate(shell%area(shell%ncell))
allocate(shell%blocknumber(shell%ncell))

counter=0
do ip=1,nblock*nnp
   if (.not.doubble(ip)) then 
      counter=counter+1
      shell%x(counter)=tempx(ip)
      shell%y(counter)=tempy(ip)
      shell%z(counter)=tempz(ip)
   end if
end do

do ib=1,nblock
shell%icon(1:block(ib)%nv,(ib-1)*nncell+1:(ib-1)*nncell+nncell)=block(ib)%icon+(ib-1)*nnp
shell%blocknumber((ib-1)*nncell+1:(ib-1)*nncell+nncell)=ib
end do


do ic=1,shell%ncell
   do i=1,shell%nv
      shell%icon(i,ic)=pointto(shell%icon(i,ic))
   end do
end do

allocate(compact(nblock*nnp))

counter=0
do ip=1,nblock*nnp
   if (.not.doubble(ip)) then
      counter=counter+1
      compact(ip)=counter
   end if
end do

do ic=1,shell%ncell
   do i=1,shell%nv
      shell%icon(i,ic)=compact(shell%icon(i,ic))
   end do
end do

end if

call cpu_time(t4) ; write(*,'(a,f8.3,a)') 'merge_blocks:',t4-t3,'s'

end subroutine

!==================================================================================================!
!==================================================================================================!

