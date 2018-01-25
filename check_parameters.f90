!==================================================================================================!
!  ***  *  *  ***   **** ***** !                                                                   !
! *     *  * *   * *       *   !                                                                   !
! *  *  **** *   *  ***    *   !                                                                   !
! *   * *  * *   *     *   *   !                                                                   !
!  ***  *  *  ***  ****    *   !                                                       C. Thieulot !
!==================================================================================================!

subroutine check_parameters

use structures
use s40rts_tomography

implicit none

!==================================================================================================!
!@@ \subsection{\tt check\_parameters}
!@@ This routine performs a basic check on the nblock, level and nlayer parameter values.
!==================================================================================================!

write(iunit,*) 'nblock=',nblock
if (level>1) then
   write(iunit,*) 'level=',level
else
   stop 'level too small'
end if

if (nlayer>1) then
   write(iunit,*) 'nlayer=',nlayer
else
   stop 'nlayer too small'
end if

if (s40rts) then
   R1=rcmb
   R2=rmoho
   write(iunit,*) 'R1=',R1
   write(iunit,*) 'R2=',R2
end if

end subroutine

!==================================================================================================!
!==================================================================================================!
