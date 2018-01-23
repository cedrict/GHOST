!==================================================================================================!
!  ***  *  *  ***   **** ***** !                                                                   !
! *     *  * *   * *       *   !                                                                   !
! *  *  **** *   *  ***    *   !                                                                   !
! *   * *  * *   *     *   *   !                                                                   !
!  ***  *  *  ***  ****    *   !                                                       C. Thieulot !
!==================================================================================================!

subroutine generate_quadrature_points

use structures

implicit none

integer inoode(8),ic,nv
integer iq1,iq2,iq3,iq,iiq,counter
real(8) rq,sq,tq

call cpu_time(t3)

!==================================================================================================!
!@@ \subsection{\tt generate\_quadrature\_points}
!@@ \begin{center}
!@@ \includegraphics[width=6cm]{images/elements}
!@@ \end{center}
!@@ For hexahedral elements the linear shape functions are:
!@@ \begin{eqnarray}
!@@ N_1(r,s,t) &=& (1-r)(1-s)(1-t)/8 \nonumber\\ 
!@@ N_2(r,s,t) &=& (1+r)(1-s)(1-t)/8 \nonumber\\ 
!@@ N_3(r,s,t) &=& (1+r)(1+s)(1-t)/8 \nonumber\\ 
!@@ N_4(r,s,t) &=& (1-r)(1+s)(1-t)/8 \nonumber\\ 
!@@ N_5(r,s,t) &=& (1-r)(1-s)(1+t)/8 \nonumber\\ 
!@@ N_6(r,s,t) &=& (1+r)(1-s)(1+t)/8 \nonumber\\ 
!@@ N_7(r,s,t) &=& (1+r)(1+s)(1+t)/8 \nonumber\\ 
!@@ N_8(r,s,t) &=& (1-r)(1+s)(1+t)/8 \nonumber
!@@ \end{eqnarray}
!@@ For the triangular prism elements the linear shape functions are:
!@@ \begin{eqnarray}
!@@ N_1(r,s,t) &=& (1-r-s)(1-t)/2 \nonumber\\ 
!@@ N_2(r,s,t) &=& r(1-t)/2 \nonumber\\ 
!@@ N_3(r,s,t) &=& s(1-t)/2 \nonumber\\ 
!@@ N_4(r,s,t) &=& (1-r-s)(1+t)/2 \nonumber\\ 
!@@ N_5(r,s,t) &=& r(1+t)/2 \nonumber\\ 
!@@ N_6(r,s,t) &=& s(1+t)/2 \nonumber
!@@ \end{eqnarray}
!==================================================================================================!

select case(mtype)

case(06,12)

   weight=1

   nv=8
   hollow_sphere%nqel=8

   hollow_sphere%nq=hollow_sphere%nqel*hollow_sphere%ncell

   allocate(hollow_sphere%xq(hollow_sphere%nq))
   allocate(hollow_sphere%yq(hollow_sphere%nq))
   allocate(hollow_sphere%zq(hollow_sphere%nq))

   counter=0    
   do iq3=1,2    
   do iq2=1,2    
   do iq1=1,2    
      counter=counter+1   
      rq=qcoords2(iq1)    
      sq=qcoords2(iq2)    
      tq=qcoords2(iq3)    
      call compute_NNN_Q1(rq,sq,tq,nv,NNN(1:nv,counter))
      call compute_dNNNdr_Q1(rq,sq,tq,nv,dNNNdr(1:nv,counter))
      call compute_dNNNds_Q1(rq,sq,tq,nv,dNNNds(1:nv,counter))
      call compute_dNNNdt_Q1(rq,sq,tq,nv,dNNNdt(1:nv,counter))
   end do    
   end do    
   end do    

   do ic=1,hollow_sphere%ncell
      inoode(1:nv)=hollow_sphere%icon(1:nv,ic)
      do iq=1,hollow_sphere%nqel
      iiq=(ic-1)*hollow_sphere%nqel+iq    
         hollow_sphere%xq(iiq)=sum(NNN(1:nv,iq)*hollow_sphere%x(inoode(1:nv)))
         hollow_sphere%yq(iiq)=sum(NNN(1:nv,iq)*hollow_sphere%y(inoode(1:nv)))
         hollow_sphere%zq(iiq)=sum(NNN(1:nv,iq)*hollow_sphere%z(inoode(1:nv)))
      end do
   end do

case(20)

   weight=1./6.

   ! quadrature points in plane rs are (1/6, 2/3),(1/6, 1/6),(2/3, 1/6)

   nv=6
   hollow_sphere%nqel=6

   hollow_sphere%nq=hollow_sphere%nqel*hollow_sphere%ncell

   allocate(hollow_sphere%xq(hollow_sphere%nq))
   allocate(hollow_sphere%yq(hollow_sphere%nq))
   allocate(hollow_sphere%zq(hollow_sphere%nq))

   do iq=1,hollow_sphere%nqel
      select case(iq)
      case(1) ; rq= 1.d0/6.d0  ; sq= 2.d0/3.d0 ; tq=qcoords2(1)
      case(2) ; rq= 2.d0/3.d0  ; sq= 1.d0/6.d0 ; tq=qcoords2(1)
      case(3) ; rq= 1.d0/6.d0  ; sq= 1.d0/6.d0 ; tq=qcoords2(1)
      case(4) ; rq= 1.d0/6.d0  ; sq= 2.d0/3.d0 ; tq=qcoords2(2)
      case(5) ; rq= 2.d0/3.d0  ; sq= 1.d0/6.d0 ; tq=qcoords2(2)
      case(6) ; rq= 1.d0/6.d0  ; sq= 1.d0/6.d0 ; tq=qcoords2(2)
      end select
      call compute_NNN_T1(rq,sq,tq,nv,NNN(1:nv,iq))
      call compute_dNNNdr_T1(rq,sq,tq,nv,dNNNdr(1:nv,iq))
      call compute_dNNNds_T1(rq,sq,tq,nv,dNNNds(1:nv,iq))
      call compute_dNNNdt_T1(rq,sq,tq,nv,dNNNdt(1:nv,iq))
   end do

   do ic=1,hollow_sphere%ncell
      inoode(1:nv)=hollow_sphere%icon(1:nv,ic)
      do iq=1,hollow_sphere%nqel
      iiq=(ic-1)*hollow_sphere%nqel+iq    
         hollow_sphere%xq(iiq)=sum(NNN(1:nv,iq)*hollow_sphere%x(inoode(1:nv)))
         hollow_sphere%yq(iiq)=sum(NNN(1:nv,iq)*hollow_sphere%y(inoode(1:nv)))
         hollow_sphere%zq(iiq)=sum(NNN(1:nv,iq)*hollow_sphere%z(inoode(1:nv)))
      end do
   end do

case default

   stop 'generate_quadrature_points: mtype not supported'

end select

!end if

call cpu_time(t4) ; write(*,'(a,f6.3,a)') 'generate_qpts:',t4-t3,'s'

end subroutine

!==================================================================================================!
!==================================================================================================!

subroutine compute_NNN_Q1(r,s,t,nv,NNN)

implicit none

real(8), intent(in) :: r,s,t
integer, intent(in) :: nv
real(8), intent(out) :: NNN(nv)

NNN(1)=0.125d0*(1.d0-r)*(1.d0-s)*(1.d0-t)                     
NNN(2)=0.125d0*(1.d0+r)*(1.d0-s)*(1.d0-t)                   
NNN(3)=0.125d0*(1.d0+r)*(1.d0+s)*(1.d0-t)                 
NNN(4)=0.125d0*(1.d0-r)*(1.d0+s)*(1.d0-t)               
NNN(5)=0.125d0*(1.d0-r)*(1.d0-s)*(1.d0+t)             
NNN(6)=0.125d0*(1.d0+r)*(1.d0-s)*(1.d0+t)           
NNN(7)=0.125d0*(1.d0+r)*(1.d0+s)*(1.d0+t)         
NNN(8)=0.125d0*(1.d0-r)*(1.d0+s)*(1.d0+t)       

end subroutine

!==============================================================================!

subroutine compute_dNNNdr_Q1(r,s,t,mpe,dNNNdr)

implicit none

real(8), intent(in) :: r,s,t
integer, intent(in) :: mpe
real(8), intent(out) :: dNNNdr(mpe)

dNNNdr(1)=-0.125d0*(1.d0-s)*(1.d0-t)                     
dNNNdr(2)=+0.125d0*(1.d0-s)*(1.d0-t)                   
dNNNdr(3)=+0.125d0*(1.d0+s)*(1.d0-t)                 
dNNNdr(4)=-0.125d0*(1.d0+s)*(1.d0-t)               
dNNNdr(5)=-0.125d0*(1.d0-s)*(1.d0+t)             
dNNNdr(6)=+0.125d0*(1.d0-s)*(1.d0+t)           
dNNNdr(7)=+0.125d0*(1.d0+s)*(1.d0+t)         
dNNNdr(8)=-0.125d0*(1.d0+s)*(1.d0+t)       

end subroutine

!==============================================================================!

subroutine compute_dNNNds_Q1(r,s,t,mpe,dNNNds)

implicit none

real(8), intent(in) :: r,s,t
integer, intent(in) :: mpe
real(8), intent(out) :: dNNNds(mpe)

dNNNds(1)=-0.125d0*(1.d0-r)*(1.d0-t)                     
dNNNds(2)=-0.125d0*(1.d0+r)*(1.d0-t)                   
dNNNds(3)=+0.125d0*(1.d0+r)*(1.d0-t)                 
dNNNds(4)=+0.125d0*(1.d0-r)*(1.d0-t)               
dNNNds(5)=-0.125d0*(1.d0-r)*(1.d0+t)             
dNNNds(6)=-0.125d0*(1.d0+r)*(1.d0+t)           
dNNNds(7)=+0.125d0*(1.d0+r)*(1.d0+t)         
dNNNds(8)=+0.125d0*(1.d0-r)*(1.d0+t)       

end subroutine

!==============================================================================!

subroutine compute_dNNNdt_Q1(r,s,t,mpe,dNNNdt)

implicit none

real(8), intent(in) :: r,s,t
integer, intent(in) :: mpe
real(8), intent(out) :: dNNNdt(mpe)
   
dNNNdt(1)=-0.125d0*(1.d0-r)*(1.d0-s)
dNNNdt(2)=-0.125d0*(1.d0+r)*(1.d0-s)
dNNNdt(3)=-0.125d0*(1.d0+r)*(1.d0+s)
dNNNdt(4)=-0.125d0*(1.d0-r)*(1.d0+s)
dNNNdt(5)=+0.125d0*(1.d0-r)*(1.d0-s)
dNNNdt(6)=+0.125d0*(1.d0+r)*(1.d0-s)
dNNNdt(7)=+0.125d0*(1.d0+r)*(1.d0+s)
dNNNdt(8)=+0.125d0*(1.d0-r)*(1.d0+s)

end subroutine

!==============================================================================!

subroutine compute_NNN_T1(r,s,t,nv,NNN)

implicit none

real(8), intent(in) :: r,s,t
integer, intent(in) :: nv
real(8), intent(out) :: NNN(nv)

NNN(1)=(1-r-s)*0.5d0*(1.d0-t)
NNN(2)=(r    )*0.5d0*(1.d0-t)
NNN(3)=(s    )*0.5d0*(1.d0-t)
NNN(4)=(1-r-s)*0.5d0*(1.d0+t)
NNN(5)=(r    )*0.5d0*(1.d0+t)
NNN(6)=(s    )*0.5d0*(1.d0+t)

end subroutine

!==============================================================================!

subroutine compute_dNNNdr_T1(r,s,t,nv,dNNNdr)

implicit none

real(8), intent(in) :: r,s,t
integer, intent(in) :: nv
real(8), intent(out) :: dNNNdr(nv)

dNNNdr(1)=-0.5d0*(1.d0-t)
dNNNdr(2)=+0.5d0*(1.d0-t)
dNNNdr(3)=0
dNNNdr(4)=-0.5d0*(1.d0+t)
dNNNdr(5)=+0.5d0*(1.d0+t)
dNNNdr(6)=0

end subroutine

!==============================================================================!

subroutine compute_dNNNds_T1(r,s,t,nv,dNNNds)

implicit none

real(8), intent(in) :: r,s,t
integer, intent(in) :: nv
real(8), intent(out) :: dNNNds(nv)

dNNNds(1)=-0.5d0*(1.d0-t)
dNNNds(2)=0
dNNNds(3)=+0.5d0*(1.d0-t)
dNNNds(4)=-0.5d0*(1.d0+t)
dNNNds(5)=0
dNNNds(6)=+0.5d0*(1.d0+t)

end subroutine

!==============================================================================!

subroutine compute_dNNNdt_T1(r,s,t,nv,dNNNdt)

implicit none

real(8), intent(in) :: r,s,t
integer, intent(in) :: nv
real(8), intent(out) :: dNNNdt(nv)

dNNNdt(1)=-(1-r-s)*0.5d0
dNNNdt(2)=-(r    )*0.5d0
dNNNdt(3)=-(s    )*0.5d0
dNNNdt(4)=+(1-r-s)*0.5d0
dNNNdt(5)=+(r    )*0.5d0
dNNNdt(6)=+(s    )*0.5d0

end subroutine

!==================================================================================================!





