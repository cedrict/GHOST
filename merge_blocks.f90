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
logical(1),dimension(:),allocatable :: sides
real(8),dimension(:),allocatable :: tempx,tempy,tempz
integer(4), dimension(:),allocatable :: pointto
logical,dimension(:),allocatable :: doubble
real(8) distance,gxip,gyip,gzip
integer, dimension(:), allocatable :: compact

!==================================================================================================!
!@@ \subsection{\tt merge\_blocks}
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

distance=1.d-6*1!outer_radius

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

