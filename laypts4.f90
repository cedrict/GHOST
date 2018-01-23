!==================================================================================================!
!  ***  *  *  ***   **** ***** !                                                                   !
! *     *  * *   * *       *   !                                                                   !
! *  *  **** *   *  ***    *   !                                                                   !
! *   * *  * *   *     *   *   !                                                                   !
!  ***  *  *  ***  ****    *   !                                                       C. Thieulot !
!==================================================================================================!

subroutine laypts4 (x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,np,x,y,z,hull)

use structures

implicit none
integer, intent(in) :: np
real(8), intent(in) :: x1,x2,x3,x4,y1,y2,y3,y4,z1,z2,z3,z4
real(8), intent(out) :: x(np)
real(8), intent(out) :: y(np)
real(8), intent(out) :: z(np)
logical, intent(out) :: hull(np)

integer counter,i,j
real(8) r,s,N1,N2,N3,N4,x0,y0

!==================================================================================================!

counter=0
do j=0,level
do i=0,level
   counter=counter+1

   !equidistant
   r=-1.d0+2.d0/level*i
   s=-1.d0+2.d0/level*j
   !write(567,*) r,s

   !equiangular
   if (equiangular) then
   x0=-pi4+(i)*2*pi4/level
   y0=-pi4+(j)*2*pi4/level
   r=tan(x0)
   s=tan(y0)
   end if
   !write(678,*) r,s

   N1=0.25d0*(1.d0-r)*(1.d0-s)
   N2=0.25d0*(1.d0+r)*(1.d0-s)
   N3=0.25d0*(1.d0+r)*(1.d0+s)
   N4=0.25d0*(1.d0-r)*(1.d0+s)
   x(counter)=x1*N1+x2*N2+x3*N3+x4*N4
   y(counter)=y1*N1+y2*N2+y3*N3+y4*N4
   z(counter)=z1*N1+z2*N2+z3*N3+z4*N4

   if (i==0) hull(counter)=.true.
   if (j==0) hull(counter)=.true.
   if (i==level) hull(counter)=.true.
   if (j==level) hull(counter)=.true.

end do
end do

end subroutine

!==================================================================================================!
!==================================================================================================!
