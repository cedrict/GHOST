!==================================================================================================!
!  ***  *  *  ***   **** ***** !                                                                   !
! *     *  * *   * *       *   !                                                                   !
! *  *  **** *   *  ***    *   !                                                                   !
! *   * *  * *   *     *   *   !                                                                   !
!  ***  *  *  ***  ****    *   !                                                       C. Thieulot !
!==================================================================================================!

subroutine project_on_sphere(radius,x,y,z,r,theta,phi)

implicit none

real(8), intent(in) :: radius
real(8), intent(inout) :: x,y,z
real(8), intent(out) :: r,theta,phi

!==================================================================================================!
!@@ \subsection{\tt project\_on\_sphere}
!@@ This subroutine gets a radius $r$ value as argument and the current coordinates
!@@ of a point. It returns the new coordinates of the point on a sphere 
!@@ of radius $r$ with the same $\theta,\phi$.
!==================================================================================================!

r=sqrt(x**2+y**2+z**2)

theta=atan2(y,x)

phi=acos(z/r)

r=radius

x=r*cos(theta)*sin(phi)
y=r*sin(theta)*sin(phi)
z=r*cos(phi)

end subroutine

!==================================================================================================!
