!==================================================================================================!
!  ***  *  *  ***   **** ***** !                                                                   !
! *     *  * *   * *       *   !                                                                   !
! *  *  **** *   *  ***    *   !                                                                   !
! *   * *  * *   *     *   *   !                                                                   !
!  ***  *  *  ***  ****    *   !                                                       C. Thieulot !
!==================================================================================================!

subroutine compute_PREM_density(radius,densprem)

implicit none
    
real(8), intent(in) :: radius 
real(8), intent(out):: densprem
real(8) x

!==================================================================================================!
!@@ \subsection{\tt compute\_PREM\_density}
!@@ This routine gets a radius $r<6371$km and returns the density 
!@@ as given by the PREM model \cite{dzan81}. 
!==================================================================================================!


x=radius/6371.d3

if (radius>6371d3) then

   densprem=0

elseif (radius<=3630.d3) then

   densprem=7.9565-6.4761*x+5.5283*x**2-3.0807*x**3

elseif (radius<=5600.d3 ) then

   densprem=7.9565-6.4761*x+5.5283*x**2-3.0807*x**3

elseif (radius<=5701.d3 ) then

   densprem=7.9565-6.4761*x+5.5283*x**2-3.0807*x**3

elseif (radius<=5771.d3 ) then

   densprem=5.3197-1.4836*x

elseif (radius<=5971.d3 ) then

   densprem=11.2494-8.0298*x

elseif (radius<=6151.d3 ) then

   densprem=7.1089-3.8045*x

elseif (radius<=6291.d3 ) then

   densprem=2.6910+0.6924*x

elseif (radius<=6346.d3 ) then

   densprem=2.6910+0.6924*x

elseif (radius<=6356.d3 ) then

   densprem=2.9

elseif (radius<=6368.d3 ) then

   densprem=2.6

else

   densprem=1.020

end if

densprem=densprem*1000

end subroutine

!==================================================================================================!
