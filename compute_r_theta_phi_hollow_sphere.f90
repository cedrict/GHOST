!==================================================================================================!
!  ***  *  *  ***   **** ***** !                                                                   !
! *     *  * *   * *       *   !                                                                   !
! *  *  **** *   *  ***    *   !                                                                   !
! *   * *  * *   *     *   *   !                                                                   !
!  ***  *  *  ***  ****    *   !                                                       C. Thieulot !
!==================================================================================================!

subroutine compute_r_theta_phi_hollow_sphere

use structures

implicit none

integer ip

!==================================================================================================!
!@@ \subsection{\tt compute\_r\_theta\_phi\_hollow\_sphere}
!@@ This subroutine computes $r$, $\theta$ and $\phi$ from $x,y,z$ for all 
!@@ nodes of the hollow\_sphere mesh with
!@@ \[
!@@ r=\sqrt{x^2+y^2+z^2}
!@@ \]
!@@ \[
!@@ \theta=\cos^{-1}(z/r)
!@@ \]
!@@ \[
!@@ \phi=\tan^{-1}_2 (y/x)
!@@ \]
!==================================================================================================!

call cpu_time(t3)

do ip=1,hollow_sphere%np
   hollow_sphere%r(ip)=sqrt(hollow_sphere%x(ip)**2+hollow_sphere%y(ip)**2+hollow_sphere%z(ip)**2)
   hollow_sphere%phi(ip)=atan2(hollow_sphere%y(ip),hollow_sphere%x(ip))
   hollow_sphere%theta(ip)=acos(hollow_sphere%z(ip)/hollow_sphere%r(ip))
end do

call cpu_time(t4) ; write(*,'(a,f6.3,a)') 'compute_r_theta_phi_hollow_sphere:',t4-t3,'s'

end subroutine 

!==================================================================================================!
!==================================================================================================!
