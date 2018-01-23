!==================================================================================================!
!  ***  *  *  ***   **** ***** !                                                                   !
! *     *  * *   * *       *   !                                                                   !
! *  *  **** *   *  ***    *   !                                                                   !
! *   * *  * *   *     *   *   !                                                                   !
!  ***  *  *  ***  ****    *   !                                                       C. Thieulot !
!==================================================================================================!

subroutine compute_gravity_on_line 

use structures

implicit none

integer i
integer, parameter :: npts=256
real(8) r
real(8) theta_g,phi_g
real(8), dimension(:), allocatable :: x,y,z
real(8), dimension(:), allocatable :: gx,gy,gz,U,g
real(8), dimension(:), allocatable :: g_theory,U_theory

!==================================================================================================!
!@@ \subsection{\tt compute\_gravity\_on\_line}
!@@ This subroutine computes on a line parametrised by $0\leq r \leq 2R_2$,
!@@ $\theta=13^\circ$ and $\phi=17^\circ$. This line is discretised 
!@@ over npts points. At each point is the subroutine compute\_gravity\_at\_point called.
!==================================================================================================!

if (compute_gravity) then

allocate(x(npts))
allocate(y(npts))
allocate(z(npts))
allocate(gx(npts))
allocate(gy(npts))
allocate(gz(npts))
allocate(g(npts))
allocate(g_theory(npts))
allocate(U_theory(npts))
allocate(U(npts))

call cpu_time(t3)

open(unit=123,file='gravity_along_a_line.dat',status='replace',action='write')

do i=1,npts

   r=2*(i-1)*R2/(npts-1)
   theta_g=13./180.*pi
   phi_g=17./180.*pi

   x(i)=r*sin(theta_g)*cos(phi_g)
   y(i)=r*sin(theta_g)*sin(phi_g)
   z(i)=r*cos(theta_g)

   call compute_gravity_at_point(x(i),y(i),z(i),gx(i),gy(i),gz(i),U(i))

   g(i)=sqrt(gx(i)**2+gy(i)**2+gz(i)**2)

   if (r<=R1) then
      g_theory(i)=0
      U_theory(i)=2.d0*pi*Ggrav*rho*(R1**2-R2**2)
   elseif (r<=R2) then
      g_theory(i)=4.d0/3.d0*pi*Ggrav*(r-R1**3/r**2)*rho
      U_theory(i)=4.d0/3.d0*pi*Ggrav*(r**2/2+R1**3/r)*rho-2*pi*rho*Ggrav*R2**2
   else
      g_theory(i)=4.d0/3.d0*pi*Ggrav*(R2**3-R1**3)/r**2 *rho
      U_theory(i)=-4.d0/3.d0*pi*Ggrav*(R2**3-R1**3)/r *rho
   end if

   write(123,*) r,gx(i),gy(i),gz(i),g(i),U(i),&
                g_theory(i),g(i)-g_theory(i),&
                U_theory(i),U(i)-U_theory(i)

end do

close(123)

write(444,*) hollow_sphere%ncell,&
             maxval(abs(g-g_theory)),&
             sum(abs(g-g_theory)),&
             sqrt(sum((g-g_theory)**2)),&
             maxval(abs(U-U_theory)),&
             sum(abs(U-U_theory)),&
             sqrt(sum((U-U_theory)**2))

end if

call cpu_time(t4) ; write(*,'(a,f10.3,a)') 'compute_gravity_on_line:',t4-t3,'s'

end subroutine

!==================================================================================================!
!==================================================================================================!
