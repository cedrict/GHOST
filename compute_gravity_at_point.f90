!==================================================================================================!
!  ***  *  *  ***   **** ***** !                                                                   !
! *     *  * *   * *       *   !                                                                   !
! *  *  **** *   *  ***    *   !                                                                   !
! *   * *  * *   *     *   *   !                                                                   !
!  ***  *  *  ***  ****    *   !                                                       C. Thieulot !
!==================================================================================================!

subroutine compute_gravity_at_point(x,y,z,gx,gy,gz,U)

use structures

implicit none

real(8), intent(in) :: x,y,z
real(8), intent(out) :: gx,gy,gz,U

integer ic,iiq,iq,k,inoode(8),nqel,nv
real(8) KK,dist,jcob,JxW
real(8) jcb3D(3,3)

!==================================================================================================!
!@@ \subsection{\tt compute\_gravity\_at\_point}
!@@ This subroutine receives as argument the coordinates {\tt x,y,z} of a point and 
!@@ returns the gravity vector components {\tt gx,gy,gz} at this location 
!@@ by means of 
!@@ \[
!@@ {\bm g}({\bm r}') = \int_\Omega {\cal G} \frac{\rho({\bm r})-\rho_0}{|{\bm r}'-{\bm r}|^3} ({\bm r}'-{\bm r}) d\Omega
!@@ \] 
!@@ where $\rho_0=${\tt refdensgrav}.
!@@ The gravitational potential is computed as follows:
!@@ \[
!@@ U({\bm r}') = -\int_\Omega {\cal G} \frac{\rho({\bm r})-\rho_0}{|{\bm r}'-{\bm r}|}  d\Omega
!@@ \]
!==================================================================================================!

nv=hollow_sphere%nv
nqel=hollow_sphere%nqel

gx=0.d0               
gy=0.d0              
gz=0.d0 
U=0.d0            
do ic=1,hollow_sphere%ncell
   inoode(1:nv)=hollow_sphere%icon(1:nv,ic)
   do iq=1,nqel                                              
      iiq=(ic-1)*nqel+iq                             
      jcb3D=0.d0    
      do k=1,nv
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
      dist=sqrt( (x-hollow_sphere%xq(iiq))**2 &
               + (y-hollow_sphere%yq(iiq))**2 &
               + (z-hollow_sphere%zq(iiq))**2 )  
      KK=rho*JxW/dist**3
      gx=gx+KK*(x-hollow_sphere%xq(iiq))
      gy=gy+KK*(y-hollow_sphere%yq(iiq))
      gz=gz+KK*(z-hollow_sphere%zq(iiq))
      U=U-rho*JxW/dist
   end do
end do                                       

gx=gx*Ggrav
gy=gy*Ggrav
gz=gz*Ggrav
U=U*Ggrav

end subroutine

!==================================================================================================!
!==================================================================================================!
