!==================================================================================================!
!  ***  *  *  ***   **** ***** !                                                                   !
! *     *  * *   * *       *   !                                                                   !
! *  *  **** *   *  ***    *   !                                                                   !
! *   * *  * *   *     *   *   !                                                                   !
!  ***  *  *  ***  ****    *   !                                                       C. Thieulot !
!==================================================================================================!

subroutine compute_shell_area

use structures

implicit none

integer ic
real(8) x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,xC,yC,zC
real(8) area1,area2,area3,area4,area_stddev,area_avrg

!==================================================================================================!
!@@ \subsection{\tt compute\_shell\_area}
!==================================================================================================!

call cpu_time(t3)

select case(mtype)
case(06,12) ! squares on sphere
   do ic=1,shell%ncell
      x1=shell%x(shell%icon(1,ic)) ; y1=shell%y(shell%icon(1,ic)) ; z1=shell%z(shell%icon(1,ic)) 
      x2=shell%x(shell%icon(2,ic)) ; y2=shell%y(shell%icon(2,ic)) ; z2=shell%z(shell%icon(2,ic)) 
      x3=shell%x(shell%icon(3,ic)) ; y3=shell%y(shell%icon(3,ic)) ; z3=shell%z(shell%icon(3,ic)) 
      x4=shell%x(shell%icon(4,ic)) ; y4=shell%y(shell%icon(4,ic)) ; z4=shell%z(shell%icon(4,ic)) 
      xC=(x1+x2+x3+x4)/4.d0
      yC=(y1+y2+y3+y4)/4.d0
      zC=(z1+z2+z3+z4)/4.d0
      call compute_triangle_area(xC,yC,zC,x1,y1,z1,x2,y2,z2,area1) ! triangle C12
      call compute_triangle_area(xC,yC,zC,x2,y2,z2,x3,y3,z3,area2) ! triangle C23
      call compute_triangle_area(xC,yC,zC,x3,y3,z3,x4,y4,z4,area3) ! triangle C34
      call compute_triangle_area(xC,yC,zC,x4,y4,z4,x1,y1,z1,area4) ! triangle C41
      shell%area(ic)=area1+area2+area3+area4
   end do
case(20) ! triangles on sphere
   do ic=1,shell%ncell
      x1=shell%x(shell%icon(1,ic)) ; y1=shell%y(shell%icon(1,ic)) ; z1=shell%z(shell%icon(1,ic)) 
      x2=shell%x(shell%icon(2,ic)) ; y2=shell%y(shell%icon(2,ic)) ; z2=shell%z(shell%icon(2,ic)) 
      x3=shell%x(shell%icon(3,ic)) ; y3=shell%y(shell%icon(3,ic)) ; z3=shell%z(shell%icon(3,ic)) 
      call compute_triangle_area(x1,y1,z1,x2,y2,z2,x3,y3,z3,shell%area(ic))
   end do
case default
   stop 'compute_shell_area: mtype not supported'
end select

write(iunit,*) 'shell area (m/M)',minval(shell%area),maxval(shell%area)
write(iunit,*) 'shell area (ratio)',maxval(shell%area)/minval(shell%area)
write(iunit,*) 'shell area (total)',sum(shell%area),4*pi,abs(sum(shell%area)-4*pi)/4/pi

area_avrg=sum(shell%area)/shell%ncell
area_stddev=sqrt( sum((shell%area-area_avrg)**2)/shell%ncell   )

write(iunit,*) 'shell area (stddev)',area_stddev

write(888,*) shell%np,&
             sum(shell%area),&
             4*pi,&
             abs(sum(shell%area)-4*pi)/4/pi,&
             area_avrg,&
             area_stddev ; call flush(888)

call cpu_time(t4) ; write(*,'(a,f6.3,a)') 'compute_shell_area:',t4-t3,'s'

end subroutine

!==================================================================================================!
!==================================================================================================!
