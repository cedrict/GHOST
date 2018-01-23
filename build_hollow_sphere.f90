!==================================================================================================!
!  ***  *  *  ***   **** ***** !                                                                   !
! *     *  * *   * *       *   !                                                                   !
! *  *  **** *   *  ***    *   !                                                                   !
! *   * *  * *   *     *   *   !                                                                   !
!  ***  *  *  ***  ****    *   !                                                       C. Thieulot !
!==================================================================================================!

subroutine build_hollow_sphere

use structures

implicit none

integer ilayer,ip,ibeg,iend
real(8) radius,r,theta,phi

!==================================================================================================!
!@@ \subsection{\tt build\_hollow\_sphere}
!@@ This subroutine first makes a copy of the current shell into {\sl shell\_temp}
!@@ It then computes the number of mesh nodes and cells for the hollow sphere
!@@ and allocate all arrays accordingly.
!@@ It then loops over layers, places the temporary shell at the required 
!@@ radius, and appends it to the current hollow sphere mesh from the inside out.
!==================================================================================================!

call cpu_time(t3)

shell_temp%nv=shell%nv
shell_temp%np=shell%np
shell_temp%ncell=shell%ncell
allocate(shell_temp%x(shell_temp%np))
allocate(shell_temp%y(shell_temp%np))
allocate(shell_temp%z(shell_temp%np))
allocate(shell_temp%icon(shell_temp%nv,shell_temp%ncell))
allocate(shell_temp%r(shell_temp%np))
allocate(shell_temp%theta(shell_temp%np))
allocate(shell_temp%phi(shell_temp%np))

hollow_sphere%nv=shell%nv*2 
hollow_sphere%np=(nlayer+1)*shell%np 
hollow_sphere%ncell=nlayer*shell%ncell
allocate(hollow_sphere%x(hollow_sphere%np))
allocate(hollow_sphere%y(hollow_sphere%np))
allocate(hollow_sphere%z(hollow_sphere%np))
allocate(hollow_sphere%icon(hollow_sphere%nv,hollow_sphere%ncell))
allocate(hollow_sphere%r(hollow_sphere%np))
allocate(hollow_sphere%theta(hollow_sphere%np))
allocate(hollow_sphere%phi(hollow_sphere%np))
allocate(hollow_sphere%volume(hollow_sphere%ncell))
allocate(hollow_sphere%blocknumber(hollow_sphere%ncell))
allocate(hollow_sphere%inner_node(hollow_sphere%np))
allocate(hollow_sphere%outer_node(hollow_sphere%np))
allocate(hollow_sphere%nx(hollow_sphere%np))
allocate(hollow_sphere%ny(hollow_sphere%np))
allocate(hollow_sphere%nz(hollow_sphere%np))
allocate(hollow_sphere%inner_cell(hollow_sphere%ncell))
allocate(hollow_sphere%outer_cell(hollow_sphere%ncell))

write(iunit,*) 'hollow_sphere%nv=',hollow_sphere%nv
write(iunit,*) 'hollow_sphere%np=',hollow_sphere%np
write(iunit,*) 'hollow_sphere%ncell=',hollow_sphere%ncell

do ilayer=1,nlayer+1

   shell_temp=shell

   radius=R1+(R2-R1)/nlayer*(ilayer-1) 
   do ip=1,shell%np
   call project_on_sphere(radius,shell_temp%x(ip),shell_temp%y(ip),shell_temp%z(ip),r,theta,phi)
   end do

   ibeg=1+(ilayer-1)*shell%np
   iend=ilayer*shell%np
   hollow_sphere%x(ibeg:iend)=shell_temp%x
   hollow_sphere%y(ibeg:iend)=shell_temp%y
   hollow_sphere%z(ibeg:iend)=shell_temp%z

end do
  
do ilayer=1,nlayer
   ibeg=1+(ilayer-1)*shell%ncell
   iend=ilayer*shell%ncell
   hollow_sphere%icon(1:shell%nv,ibeg:iend)=shell_temp%icon(1:shell%nv,1:shell%ncell) + (ilayer-1)*shell%np
   hollow_sphere%icon(shell%nv+1:2*shell%nv,ibeg:iend)=shell_temp%icon(1:shell%nv,1:shell%ncell) +ilayer*shell%np
   hollow_sphere%blocknumber(ibeg:iend)=shell%blocknumber(1:shell%ncell)
end do

call cpu_time(t4) ; write(*,'(a,f6.3,a)') 'build_hollow_sphere:',t4-t3,'s'

end subroutine

!==================================================================================================!
!==================================================================================================!
