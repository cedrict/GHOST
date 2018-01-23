!==================================================================================================!
!  ***  *  *  ***   **** ***** !                                                                   !
! *     *  * *   * *       *   !                                                                   !
! *  *  **** *   *  ***    *   !                                                                   !
! *   * *  * *   *     *   *   !                                                                   !
!  ***  *  *  ***  ****    *   !                                                       C. Thieulot !
!==================================================================================================!

subroutine read_s40rts

use structures
use S40RTS_tomography

implicit none

integer i,j,k,ip
real(8) radius,latitude,longitude,densprem
real, external :: sdot,splh


!==================================================================================================!
!@@ \subsection{\tt read\_s40rts}
!==================================================================================================!

call cpu_time(t3)

if (s40rts) then

allocate(hollow_sphere%dlnvs(hollow_sphere%np)) ; hollow_sphere%dlnvs=0
allocate(hollow_sphere%drho(hollow_sphere%np))  ; hollow_sphere%drho=0

!------------------------------------------------------------
! reading the data file 

open(unit=10,file='DATA/S40RTS.sph',action='read',status='old')
read(10,*) lmax,dum,nsmx
ndmx = nsmx
ndmn = 4
natd=(lmax+1)**2
ind=(ndmn-1)*natd+1
do i=ndmn,ndmx
   do j=0,lmax
      ind1=ind+2*j
      read(10,'(11e12.4)',end=100)(x(k),k=ind,ind1)
      ind=ind1+1
   enddo
enddo
goto 200
100  stop 'incompatible sph header'
200 continue
close(10)
if(lmax.gt.MXLH) stop 'lmax.gt.MXLH'
if(ndmx.gt.MXSPL) stop 'ndmx.gt.MXSPL'

natd=(lmax+1)**2

!--------------------------------------------

do ip=1,hollow_sphere%np

   latitude =pi/2.d0-hollow_sphere%theta(ip)
   longitude=        hollow_sphere%phi(ip)

   latitude =latitude /pi*180.
   longitude=longitude/pi*180.

   radius=hollow_sphere%r(ip)

   d0=0 ; wk1=0 ; wk2=0 ; wk3=0

   ! Calculate the spline basis functions at a regular hollow_sphere

   call splhsetup()

   ! calculate Y_k(latitude,longitude) at random depth level

   call ylm(real(latitude),real(longitude),lmax,d0,wk1,wk2,wk3)

   ! calculate spline coefficients

   nbeg=max(4,ndmn)
   nend=ndmx-3
   do i=nbeg,ndmx
   ind=(i-1)*natd+1
   spl(i-3)=sdot(natd,d0,1,x(ind),1)
   enddo

   ! xd takes value -1 at the CMB and +1 at the moho
   ! if xd is outside this range, the splh returns 0

   xd=-1.+2.*(radius-rcmb)/(rmoho-rcmb)
   do ipp=nbeg-3,nend
   hollow_sphere%dlnvs(ip)=hollow_sphere%dlnvs(ip)+spl(ipp)*splh(ipp-1,xd)
   enddo

   call compute_PREM_density(hollow_sphere%r(ip),densprem)

   hollow_sphere%drho(ip)=densprem*xi_s40rts*hollow_sphere%dlnvs(ip)

end do

end if

call cpu_time(t4) ; write(*,'(a,f10.3,a)') 'read_s40rts:',t4-t3,'s'

end subroutine

!==================================================================================================!
