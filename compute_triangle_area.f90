!==================================================================================================!
!  ***  *  *  ***   **** ***** !                                                                   !
! *     *  * *   * *       *   !                                                                   !
! *  *  **** *   *  ***    *   !                                                                   !
! *   * *  * *   *     *   *   !                                                                   !
!  ***  *  *  ***  ****    *   !                                                       C. Thieulot !
!==================================================================================================!

subroutine compute_triangle_area(x1,y1,z1,x2,y2,z2,x3,y3,z3,area)

implicit none

real(8), intent(in) :: x1,y1,z1,x2,y2,z2,x3,y3,z3
real(8), intent(out) :: area
real(8) a,b,c,s
real(8), external :: distance

!==================================================================================================!
!@@ \subsection{\tt compute\_triangle\_area}
!@@ {\tt https://en.wikipedia.org/wiki/Heron\%27s\_formula }
!@@ Heron's formula states that the area of a triangle whose sides have lengths a, b, and c is
!@@ \[
!@@ A=\sqrt{s(s-a)(s-b)(s-c)}
!@@ \]
!@@ where s is the semiperimeter of the triangle; that is,
!@@ \[
!@@ s=\frac{a+b+c}{2}
!@@ \]
!==================================================================================================!

a=distance(x1,y1,z1,x2,y2,z2)
b=distance(x1,y1,z1,x3,y3,z3)
c=distance(x2,y2,z2,x3,y3,z3)

s=(a+b+c)/2.d0

area=sqrt(s*(s-a)*(s-b)*(s-c))

end subroutine

!==============================================================================!

function distance(x1,y1,z1,x2,y2,z2)
implicit none
real(8) x1,y1,z1,x2,y2,z2,distance

distance=sqrt( (x1-x2)**2+(y1-y2)**2+(z1-z2)**2  )

end function

!==================================================================================================!
!==================================================================================================!
