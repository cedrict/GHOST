!==================================================================================================!
!  ***  *  *  ***   **** ***** !                                                                   !
! *     *  * *   * *       *   !                                                                   !
! *  *  **** *   *  ***    *   !                                                                   !
! *   * *  * *   *     *   *   !                                                                   !
!  ***  *  *  ***  ****    *   !                                                       C. Thieulot !
!==================================================================================================!

subroutine int_to_char(charac,ncharac,integ)
implicit none
character(len=*) :: charac
integer integ,ncharac

select case (ncharac)
   case (1) 
      write (charac,'(i1)') integ
   case (2) 
      write (charac,'(i2)') integ
      if (integ.lt.10) charac(1:1)='0'
   case (3) 
      write (charac,'(i3)') integ
      if (integ.lt.100) charac(1:1)='0'
      if (integ.lt.10) charac(1:2)='00'
   case default
      stop 'value ncharac too big'
end select

end subroutine int_to_char

!==================================================================================================!
!==================================================================================================!
