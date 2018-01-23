!==================================================================================================!
!  ***  *  *  ***   **** ***** !                                                                   !
! *     *  * *   * *       *   !                                                                   !
! *  *  **** *   *  ***    *   !                                                                   !
! *   * *  * *   *     *   *   !                                                                   !
!  ***  *  *  ***  ****    *   !                                                       C. Thieulot !
!==================================================================================================!

subroutine compute_r_theta_phi_shell

use structures

implicit none

integer ip

!==================================================================================================!
!@@ \subsection{\tt compute\_r\_theta\_phi\_shell}
!@@ This subroutine computes $r$, $\theta$ and $\phi$ from $x,y,z$ for all 
!@@ nodes of the shell mesh with
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

do ip=1,shell%np
   shell%r(ip)=sqrt(shell%x(ip)**2+shell%y(ip)**2+shell%z(ip)**2)
   shell%phi(ip)=atan2(shell%y(ip),shell%x(ip))
   shell%theta(ip)=acos(shell%z(ip)/shell%r(ip))
end do

call cpu_time(t4) ; write(*,'(a,f6.3,a)') 'compute_r_theta_phi_shell:',t4-t3,'s'

end subroutine 

!==================================================================================================!
!==================================================================================================!
