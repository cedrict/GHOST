!==================================================================================================!
!  ***  *  *  ***   **** ***** !                                                                   !
! *     *  * *   * *       *   !                                                                   !
! *  *  **** *   *  ***    *   !                                                                   !
! *   * *  * *   *     *   *   !                                                                   !
!  ***  *  *  ***  ****    *   !                                                       C. Thieulot !
!==================================================================================================!

subroutine read_command_line_options

use structures

implicit none

character(len=255) :: arg
integer :: option_ID=1
integer :: argc,numarg,i
logical :: err_detected

!==================================================================================================!
!@@ \subsection{\tt read\_command\_line\_options}
!@@ This subroutine reads in the command line arguments
!==================================================================================================!

numarg = command_argument_count()

if (numarg>0) then

   argc=command_argument_count() 

   do while (option_ID <= argc)
   call getarg(option_ID,arg)

   if (arg == '-level' ) then
      option_ID=option_ID+1
      call getarg(option_ID,arg)
      read(arg,*) level
      write(*,'(a,i3)') 'read: level = ',level

   elseif (arg == '-mtype' ) then
      option_ID=option_ID+1
      call getarg(option_ID,arg)
      read(arg,*) mtype
      write(*,'(a,i3)') 'read: mtype = ',mtype

   elseif (arg == '-nlayer' ) then
      option_ID=option_ID+1
      call getarg(option_ID,arg)
      read(arg,*) nlayer
      write(*,'(a,i3)') 'read: nlayer = ',nlayer

   elseif (arg == '-equiangular' ) then
      option_ID=option_ID+1
      call getarg(option_ID,arg)
      read(arg,*) i
      equiangular=(i==1)
      write(*,'(a,i3)') 'read: equiangular = ',i
   elseif (arg == '--help' ) then
      print *,'**********'
      print *,'* SPIDER *'
      print *,'**********'
      print *,'available options:'
      print *,'-level [integer number]'
      print *,'-mtype [integer number]'
      print *,'-nlayer [integer number]'
      print *,'-equiangular [0 or 1]'
      stop ' '
   else
      err_detected=.true.
      exit
   end if

   option_ID=option_ID+1

   end do

end if

end subroutine

!==================================================================================================!
!==================================================================================================!
