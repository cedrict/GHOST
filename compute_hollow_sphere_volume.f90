!==================================================================================================!
!  ***  *  *  ***   **** ***** !                                                                   !
! *     *  * *   * *       *   !                                                                   !
! *  *  **** *   *  ***    *   !                                                                   !
! *   * *  * *   *     *   *   !                                                                   !
!  ***  *  *  ***  ****    *   !                                                       C. Thieulot !
!==================================================================================================!

subroutine compute_hollow_sphere_volume
use structures
implicit none
integer ic,k,iq,inode,inoode(8),nv
real(8) xe(8),ye(8),ze(8),Vol
real(8), external :: hexahedron_volume 
real(8) jcb3D(3,3),jcob,JxW,Vel

!==================================================================================================!
!@@ \subsection{\tt compute\_hollow\_sphere\_volume}
!@@ This routine loops over all cells and adds up their volumes. 
!@@ Each cell volume is obtaines by means of a Gauss quadrature.
!@@ Measurements are written in {\sl fort.777}. 
!==================================================================================================!

call cpu_time(t3)

   nv=hollow_sphere%nv

   do ic=1,hollow_sphere%ncell
      inoode(1:nv)=hollow_sphere%icon(1:nv,ic)
      Vel=0.d0
      do iq=1,hollow_sphere%nqel   
         jcb3D=0.d0    
         do k=1,hollow_sphere%nv
         jcb3D(1,1)=jcb3D(1,1)+dNNNdr(k,iq)*hollow_sphere%x(inoode(k)) 
         jcb3D(1,2)=jcb3D(1,2)+dNNNdr(k,iq)*hollow_sphere%y(inoode(k)) 
         jcb3D(1,3)=jcb3D(1,3)+dNNNdr(k,iq)*hollow_sphere%z(inoode(k)) 
         jcb3D(2,1)=jcb3D(2,1)+dNNNds(k,iq)*hollow_sphere%x(inoode(k)) 
         jcb3D(2,2)=jcb3D(2,2)+dNNNds(k,iq)*hollow_sphere%y(inoode(k)) 
         jcb3D(2,3)=jcb3D(2,3)+dNNNds(k,iq)*hollow_sphere%z(inoode(k))  
         jcb3D(3,1)=jcb3D(3,1)+dNNNdt(k,iq)*hollow_sphere%x(inoode(k))   
         jcb3D(3,2)=jcb3D(3,2)+dNNNdt(k,iq)*hollow_sphere%y(inoode(k))    
         jcb3D(3,3)=jcb3D(3,3)+dNNNdt(k,iq)*hollow_sphere%z(inoode(k))    
         enddo    
         jcob=jcb3D(1,1)*jcb3D(2,2)*jcb3D(3,3) &    
             +jcb3D(1,2)*jcb3D(2,3)*jcb3D(3,1) &    
             +jcb3D(2,1)*jcb3D(3,2)*jcb3D(1,3) &    
             -jcb3D(1,3)*jcb3D(2,2)*jcb3D(3,1) &    
             -jcb3D(1,2)*jcb3D(2,1)*jcb3D(3,3) &    
             -jcb3D(2,3)*jcb3D(3,2)*jcb3D(1,1)    
         JxW=jcob*weight
         Vel=Vel+ JxW * sum(NNN(1:nv,iq))
      end do
      hollow_sphere%volume(ic)=Vel
   end do

write(iunit,*) 'hollow_sphere volume (m/M)',minval(hollow_sphere%volume),maxval(hollow_sphere%volume)
write(iunit,*) 'hollow_sphere (ratio)',maxval(hollow_sphere%volume)/minval(hollow_sphere%volume)
write(iunit,*) 'hollow_sphere (total)',sum(hollow_sphere%volume),4.d0/3.d0*pi*(2**3-1**3)

Vol=4.d0/3.d0*pi*(R2**3-R1**3)

write(777,*) hollow_sphere%np,sum(hollow_sphere%volume),Vol,&
             abs(sum(hollow_sphere%volume)-Vol)/Vol

call cpu_time(t4) ; write(*,'(a,f6.3,a)') 'compute_hollow_sphere_volume:',t4-t3,'s'

end subroutine

!==================================================================================================!
!==================================================================================================!
