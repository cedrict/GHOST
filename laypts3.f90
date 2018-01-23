!==================================================================================================!
!  ***  *  *  ***   **** ***** !                                                                   !
! *     *  * *   * *       *   !                                                                   !
! *  *  **** *   *  ***    *   !                                                                   !
! *   * *  * *   *     *   *   !                                                                   !
!  ***  *  *  ***  ****    *   !                                                       C. Thieulot !
!==================================================================================================!

subroutine laypts3(x1,y1,z1,x2,y2,z2,x3,y3,z3,np,x,y,z,hull)
implicit none
real(8), intent(in) :: x1,x2,x3,y1,y2,y3,z1,z2,z3
integer, intent(in) :: np
real(8), intent(out) :: x(np)
real(8), intent(out) :: y(np)
real(8), intent(out) :: z(np)
logical, intent(out) :: hull(np)

real(8), parameter :: eps = 1.d-8 
integer ip
real(8) N1,N2,N3

hull=.false.

do ip=1,np
   if (abs(x(ip))<eps) hull(ip)=.true. 
   if (abs(y(ip))<eps) hull(ip)=.true. 
   if (abs(1.d0-x(ip)-y(ip))<eps) hull(ip)=.true. 
   N1=1-x(ip)-y(ip)
   N2=x(ip)
   N3=y(ip)
   x(ip)=x1*N1+x2*N2+x3*N3
   y(ip)=y1*N1+y2*N2+y3*N3
   z(ip)=z1*N1+z2*N2+z3*N3
end do

end subroutine

!==============================================================================!

subroutine laypts3b(theta1,phi1,theta2,phi2,theta3,phi3,np,x,y,z,hull)
implicit none
real(8), intent(in) :: theta1,phi1 
real(8), intent(in) :: theta2,phi2
real(8), intent(in) :: theta3,phi3 
integer, intent(in) :: np
real(8), intent(out) :: x(np)
real(8), intent(out) :: y(np)
real(8), intent(out) :: z(np)
logical, intent(out) :: hull(np)

real(8), parameter :: eps = 1.d-8 
integer ip
real(8) N1,N2,N3,theta,phi

hull=.false.

do ip=1,np
   if (abs(x(ip))<eps) hull(ip)=.true. 
   if (abs(y(ip))<eps) hull(ip)=.true. 
   if (abs(1.d0-x(ip)-y(ip))<eps) hull(ip)=.true. 
   N1=1-x(ip)-y(ip)
   N2=x(ip)
   N3=y(ip)
   theta=theta1*N1+theta2*N2+theta3*N3
   phi=phi1*N1+phi2*N2+phi3*N3
   ! radius is 1
   x(ip)=sin(theta)*cos(phi)
   y(ip)=sin(theta)*sin(phi)
   z(ip)=cos(theta)
end do

end subroutine

!==================================================================================================!
!==================================================================================================!
