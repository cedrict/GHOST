!==================================================================================================!
!  ***  *  *  ***   **** ***** !                                                                   !
! *     *  * *   * *       *   !                                                                   !
! *  *  **** *   *  ***    *   !                                                                   !
! *   * *  * *   *     *   *   !                                                                   !
!  ***  *  *  ***  ****    *   !                                                       C. Thieulot !
!==================================================================================================!

subroutine output_shell_lat_lon

use structures 

implicit none

integer ip

!==================================================================================================!
!@@ \subsection{\tt output\_shell\_lat\_lon}
!@@ If the {\sl generate\_ascii\_output} flag is set to true
!@@ this subroutine exports $\theta$ and $\phi$ values of all points on the 
!@@ shell in the file {\sl shell\_theta\_phi.dat} in {\sl OUT}.
!==================================================================================================!

if (generate_ascii_output) then
open(unit=777,file='OUT/shell_theta_phi.dat',action='write',status='replace')
do ip=1,shell%np
   write(777,*) shell%theta(ip),shell%phi(ip)
end do
close(777)
end if

end subroutine

!==================================================================================================!
!==================================================================================================!
