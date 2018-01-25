!==================================================================================================!
!  ***  *  *  ***   **** ***** !                                                                   !
! *     *  * *   * *       *   !                                                                   !
! *  *  **** *   *  ***    *   !                                                                   !
! *   * *  * *   *     *   *   !                                                                   !
!  ***  *  *  ***  ****    *   !                                                       C. Thieulot !
!==================================================================================================!

subroutine map_blocks

use structures

implicit none

!==============================================================================!
!@@ \subsection{\tt map\_blocks}
!@@ \subsubsection{HS06 - cubed sphere}
!@@ See section 3.1 of \cite{puli07}, or \cite{sado72,natl05}
!@@ 
!@@ \includegraphics[width=8cm]{images/sado72}

!@@ \subsubsection{HS12 - CitcomS mesh}

!@@ See \cite{zhzm00,zhmt08}.

!@@ \includegraphics[width=8cm]{images/citcoms}

!@@ \subsubsection{HS20 - Icosahedron}
!@@ Let $t=(1+\sqrt{5})/2$, then the vertices of the icosahedron are 
!@@ \[
!@@ {\bm x}_0=\frac{1}{\sqrt{1+t^2}}(t,1,0) \quad\quad {\bm x}_6=\frac{1}{\sqrt{1+t^2}}(-1,0,t)
!@@ \]
!@@ \[
!@@ {\bm x}_1=\frac{1}{\sqrt{1+t^2}}(-t,1,0) \quad\quad {\bm x}_7=\frac{1}{\sqrt{1+t^2}}(-1,0,-t)
!@@ \]
!@@ \[
!@@ {\bm x}_2=\frac{1}{\sqrt{1+t^2}}(t,-1,0) \quad\quad {\bm x}_8=\frac{1}{\sqrt{1+t^2}}(0,t,1)
!@@ \]
!@@ \[
!@@ {\bm x}_3=\frac{1}{\sqrt{1+t^2}}(-t,-1,0) \quad\quad {\bm x}_9=\frac{1}{\sqrt{1+t^2}}(0,-t,1)
!@@ \]
!@@ \[
!@@ {\bm x}_4=\frac{1}{\sqrt{1+t^2}}(1,0,t) \quad\quad {\bm x}_{10}=\frac{1}{\sqrt{1+t^2}}(0,t,-1)
!@@ \]
!@@ \[
!@@ {\bm x}_5=\frac{1}{\sqrt{1+t^2}}(1,0,-t) \quad\quad {\bm x}_{11}=\frac{1}{\sqrt{1+t^2}}(0,-t,-1)
!@@ \]
!@@ It is obvious that all points lie on a sphere of unit radius.
!@@
!@@ The triangles are as follows:
!@@ \[  T_0=[0,8,4]  \quad\quad T_{10}=[2,9,11]   \]
!@@ \[  T_1=[0,5,10] \quad\quad T_{11}=[3,11,9]   \]
!@@ \[  T_2=[2,4,9]  \quad\quad T_{12}=[4,2,0]   \]
!@@ \[  T_3=[2,11,5] \quad\quad T_{13}=[5,0,2]   \]
!@@ \[  T_4=[1,6,8]  \quad\quad T_{14}=[6,1,3]   \]
!@@ \[  T_5=[1,10,7] \quad\quad T_{15}=[7,3,1]   \]
!@@ \[  T_6=[3,9,6]  \quad\quad T_{16}=[8,6,4]   \]
!@@ \[  T_7=[3,7,11] \quad\quad T_{17}=[9,4,6]   \]
!@@ \[  T_8=[0,10,8] \quad\quad T_{18}=[10,5,7]   \]
!@@ \[  T_9=[1,8,10] \quad\quad T_{19}=[11,7,5]   \]
!@@ Graphically, the mesh has the following structure:
!@@ \begin{center}
!@@ \includegraphics[width=12cm]{images/flat_icosahedron}
!@@ \end{center}
!==============================================================================!


integer counter,i,j,ib,ip
real(8) alpha_min,alpha_extent
real(8) beta_min,beta_extent
real(8) d_alpha,d_beta,radius
real(8) alphaa,betaa,X,Y,delta
real(8) xA,yA,zA,rA,thetaA,phiA
real(8) xB,yB,zB,rB,thetaB,phiB
real(8) xC,yC,zC,rC,thetaC,phiC
real(8) xD,yD,zD,rD,thetaD,phiD
real(8) xE,yE,zE,rE,thetaE,phiE
real(8) xF,yF,zF,rF,thetaF,phiF
real(8) xG,yG,zG,rG,thetaG,phiG
real(8) xH,yH,zH,rH,thetaH,phiH
real(8) xI,yI,zI,thetaI,phiI
real(8) xJ,yJ,zJ,rJ,thetaJ,phiJ
real(8) xK,yK,zK,rK,thetaK,phiK
real(8) xL,yL,zL,thetaL,phiL
real(8) xM,yM,zM,rM,thetaM,phiM
real(8) xN,yN,zN,rN,thetaN,phiN
real(8) xP,yP,zP,rP,thetaP,phiP
real(8) xQ,yQ,zQ,rQ,thetaQ,phiQ
real(8) t,tt

call cpu_time(t3)

select case(mtype)

case(06) ! cubed sphere

   alpha_min=-pi/4.d0
   alpha_extent=pi/2.d0
   beta_min=-pi/4.d0
   beta_extent=pi/2.d0
   d_alpha=(alpha_extent)/level
   d_beta =(beta_extent)/level
   radius=1

   counter=0    
   do j=0,level
   do i=0,level
   counter=counter+1   

   !equidistant
   X=-1.d0+2.d0/level*i
   Y=-1.d0+2.d0/level*j

   !equiangular
   if (equiangular) then
   alphaa=alpha_min+d_alpha*i
   betaa =beta_min +d_beta *j
   X=tan(alphaa)
   Y=tan(betaa)
   end if

   delta=1+X**2+Y**2   
   block(1)%x(counter)= radius/sqrt(delta)  
   block(1)%y(counter)= radius/sqrt(delta)*X
   block(1)%z(counter)= radius/sqrt(delta)*Y
   block(2)%x(counter)=-radius/sqrt(delta)*X
   block(2)%y(counter)= radius/sqrt(delta)
   block(2)%z(counter)= radius/sqrt(delta)*Y
   !block(3)%x(counter)=-radius/sqrt(delta)
   !block(3)%y(counter)= radius/sqrt(delta)*X
   !block(3)%z(counter)= radius/sqrt(delta)*Y
   !block(4)%x(counter)=-radius/sqrt(delta)*X
   !block(4)%y(counter)=-radius/sqrt(delta)
   !block(4)%z(counter)= radius/sqrt(delta)*Y
   block(5)%x(counter)=-radius/sqrt(delta)*Y
   block(5)%y(counter)= radius/sqrt(delta)*X
   block(5)%z(counter)= radius/sqrt(delta)  
   !block(6)%x(counter)=-radius/sqrt(delta)*Y
   !block(6)%y(counter)= radius/sqrt(delta)*X
   !block(6)%z(counter)=-radius/sqrt(delta)  
   if (i==0 .or. i==level .or. j==0 .or. j==level) then
   block(1)%hull(counter)=.true.
   block(2)%hull(counter)=.true.
   block(3)%hull(counter)=.true.
   block(4)%hull(counter)=.true.
   block(5)%hull(counter)=.true.
   block(6)%hull(counter)=.true.
   end if
   end do    
   end do    

   ! block3 is rotation of block 1
   ! by 180 degrees around axis Z
   block(3)=block(1)
   block(3)%x=-block(1)%x
   block(3)%y=-block(1)%y

   ! block4 is rotation of block 2
   ! by 180 degrees around axis Z
   block(4)=block(2)
   block(4)%x=-block(2)%x
   block(4)%y=-block(2)%y

   ! block6 is rotation of block 5
   ! by 180 degrees around axis Y
   block(6)=block(5)
   block(6)%x=-block(5)%x
   block(6)%z=-block(5)%z

case(12) ! citcom mesh

   radius=1

   !----------------------
   ! four corners

   xA=-1.d0
   yA= 0.d0
   zA=-1.d0/sqrt(2.d0)

   xB=+1.d0
   yB= 0.d0 
   zB=-1.d0/sqrt(2.d0)

   xC= 0.d0 
   yC=-1.d0
   zC= 1.d0/sqrt(2.d0)

   xD= 0.d0 
   yD=+1.d0
   zD= 1.d0/sqrt(2.d0)

   !----------------------
   ! middles of faces

   xM=(xA+xB+xC)/3.d0
   yM=(yA+yB+yC)/3.d0
   zM=(zA+zB+zC)/3.d0

   xN=(xA+xD+xC)/3.d0
   yN=(yA+yD+yC)/3.d0
   zN=(zA+zD+zC)/3.d0

   xP=(xA+xD+xB)/3.d0
   yP=(yA+yD+yB)/3.d0
   zP=(zA+zD+zB)/3.d0

   xQ=(xC+xD+xB)/3.d0
   yQ=(yC+yD+yB)/3.d0
   zQ=(zC+zD+zB)/3.d0

   !----------------------
   ! middle of edges

   xF=(xB+xC)/2.d0
   yF=(yB+yC)/2.d0
   zF=(zB+zC)/2.d0

   xG=(xA+xC)/2.d0
   yG=(yA+yC)/2.d0
   zG=(zA+zC)/2.d0

   xE=(xB+xA)/2.d0
   yE=(yB+yA)/2.d0
   zE=(zB+zA)/2.d0

   xH=(xD+xC)/2.d0
   yH=(yD+yC)/2.d0
   zH=(zD+zC)/2.d0

   xJ=(xD+xA)/2.d0
   yJ=(yD+yA)/2.d0
   zJ=(zD+zA)/2.d0

   xK=(xD+xB)/2.d0
   yK=(yD+yB)/2.d0
   zK=(zD+zB)/2.d0

   ! making sure points A...Q are on the sphere

   call project_on_sphere(radius,xA,yA,zA,rA,thetaA,phiA)
   call project_on_sphere(radius,xB,yB,zB,rB,thetaB,phiB)
   call project_on_sphere(radius,xC,yC,zC,rC,thetaC,phiC)
   call project_on_sphere(radius,xD,yD,zD,rD,thetaD,phiD)
   call project_on_sphere(radius,xE,yE,zE,rE,thetaE,phiE)
   call project_on_sphere(radius,xF,yF,zF,rF,thetaF,phiF)
   call project_on_sphere(radius,xG,yG,zG,rG,thetaG,phiG)
   call project_on_sphere(radius,xH,yH,zH,rH,thetaH,phiH)
   call project_on_sphere(radius,xJ,yJ,zJ,rJ,thetaJ,phiJ)
   call project_on_sphere(radius,xK,yK,zK,rK,thetaK,phiK)
   call project_on_sphere(radius,xM,yM,zM,rM,thetaM,phiM)
   call project_on_sphere(radius,xN,yN,zN,rN,thetaN,phiN)
   call project_on_sphere(radius,xP,yP,zP,rP,thetaP,phiP)
   call project_on_sphere(radius,xQ,yQ,zQ,rQ,thetaQ,phiQ)

   call laypts4(xM,yM,zM,xG,yG,zG,xA,yA,zA,xE,yE,zE,block(01)%np,block(01)%x,block(01)%y,block(01)%z,block(01)%hull) ! cell #01 MGAE
   call laypts4(xF,yF,zF,xM,yM,zM,xE,yE,zE,xB,yB,zB,block(02)%np,block(02)%x,block(02)%y,block(02)%z,block(02)%hull) ! cell #02 MEBF 
   call laypts4(xC,yC,zC,xG,yG,zG,xM,yM,zM,xF,yF,zF,block(03)%np,block(03)%x,block(03)%y,block(03)%z,block(03)%hull) ! cell #03 MFCG
   call laypts4(xG,yG,zG,xN,yN,zN,xJ,yJ,zJ,xA,yA,zA,block(04)%np,block(04)%x,block(04)%y,block(04)%z,block(04)%hull) ! cell #04 NJAG
   call laypts4(xC,yC,zC,xH,yH,zH,xN,yN,zN,xG,yG,zG,block(05)%np,block(05)%x,block(05)%y,block(05)%z,block(05)%hull) ! cell #05 NGCH
   call laypts4(xH,yH,zH,xD,yD,zD,xJ,yJ,zJ,xN,yN,zN,block(06)%np,block(06)%x,block(06)%y,block(06)%z,block(06)%hull) ! cell #06 NHDJ
   call laypts4(xA,yA,zA,xJ,yJ,zJ,xP,yP,zP,xE,yE,zE,block(07)%np,block(07)%x,block(07)%y,block(07)%z,block(07)%hull) ! cell #07 PEAJ
   call laypts4(xJ,yJ,zJ,xD,yD,zD,xK,yK,zK,xP,yP,zP,block(08)%np,block(08)%x,block(08)%y,block(08)%z,block(08)%hull) ! cell #08 PJDK
   call laypts4(xP,yP,zP,xK,yK,zK,xB,yB,zB,xE,yE,zE,block(09)%np,block(09)%x,block(09)%y,block(09)%z,block(09)%hull) ! cell #09 PKBE
   call laypts4(xQ,yQ,zQ,xK,yK,zK,xD,yD,zD,xH,yH,zH,block(10)%np,block(10)%x,block(10)%y,block(10)%z,block(10)%hull) ! cell #10 QKDH
   call laypts4(xQ,yQ,zQ,xH,yH,zH,xC,yC,zC,xF,yF,zF,block(11)%np,block(11)%x,block(11)%y,block(11)%z,block(11)%hull) ! cell #11 QHCF
   call laypts4(xQ,yQ,zQ,xF,yF,zF,xB,yB,zB,xK,yK,zK,block(12)%np,block(12)%x,block(12)%y,block(12)%z,block(12)%hull) ! cell #12 QFBK

   ! make sure all points end up on sphere

   do ib=1,nblock
      do ip=1,block(ib)%np
         call project_on_sphere(radius,block(ib)%x(ip),block(ib)%y(ip),block(ib)%z(ip),&
                                block(ib)%r(ip),block(ib)%theta(ip),block(ib)%phi(ip))
      end do
   end do

case(20) ! icosahedron

   t=(1+sqrt(5.d0))/2.d0
   tt=1./sqrt(1+t**2)

   xA=+t*tt ; yA=+1*tt ; zA=0     ; thetaA=acos(zA) ; phiA=atan2(yA,xA)
   xB=-t*tt ; yB=+1*tt ; zB=0     ; thetaB=acos(zB) ; phiB=atan2(yB,xB)
   xC=+t*tt ; yC=-1*tt ; zC=0     ; thetaC=acos(zC) ; phiC=atan2(yC,xC)
   xD=-t*tt ; yD=-1*tt ; zD=0     ; thetaD=acos(zD) ; phiD=atan2(yD,xD)
   xE=+1*tt ; yE=0     ; zE=+t*tt ; thetaE=acos(zE) ; phiE=atan2(yE,xE)
   xF=+1*tt ; yF=0     ; zF=-t*tt ; thetaF=acos(zF) ; phiF=atan2(yF,xF)
   xG=-1*tt ; yG=0     ; zG=+t*tt ; thetaG=acos(zG) ; phiG=atan2(yG,xG)
   xH=-1*tt ; yH=0     ; zH=-t*tt ; thetaH=acos(zH) ; phiH=atan2(yH,xH)
   xI=0     ; yI=+t*tt ; zI=+1*tt ; thetaI=acos(zI) ; phiI=atan2(yI,xI)
   xJ=0     ; yJ=-t*tt ; zJ=+1*tt ; thetaJ=acos(zJ) ; phiJ=atan2(yJ,xJ)
   xK=0     ; yK=+t*tt ; zK=-1*tt ; thetaK=acos(zK) ; phiK=atan2(yK,xK)
   xL=0     ; yL=-t*tt ; zL=-1*tt ; thetaL=acos(zL) ; phiL=atan2(yL,xL)

   !A B C D E F G H I J K  L
   !0 1 2 3 4 5 6 7 8 9 10 11

   ! T00 0-8-4  : AIE
   ! T01 0-5-10 : AFK
   ! T02 2-4-9  : CEJ
   ! T03 2-11-5 : CLF
   ! T04 1-6-8  : BGI
   ! T05 1-10-7 : BKH
   ! T06 3-9-6  : DJG
   ! T07 3-7-11 : DHL
   ! T08 0-10-8 : AKI
   ! T09 1-8-10 : BIK
   ! T10 2-9-11 : CJL
   ! T11 3-11-9 : DLJ
   ! T12 4-2-0  : ECA
   ! T13 5-0-2  : FAC
   ! T14 6-1-3  : GBD
   ! T15 7-3-1  : HDB
   ! T16 8-6-4  : IGE
   ! T17 9-4-6  : JEG
   ! T18 10-5-7 : KFH
   ! T19 11-7-5 : LHF

   call laypts3(xA,yA,zA,xI,yI,zI,xE,yE,zE,block(1)%np,block(1)%x,block(1)%y,block(1)%z,block(1)%hull)      ! T00
   call laypts3(xA,yA,zA,xF,yF,zF,xK,yK,zK,block(2)%np,block(2)%x,block(2)%y,block(2)%z,block(2)%hull)      ! T01
   call laypts3(xC,yC,zC,xE,yE,zE,xJ,yJ,zJ,block(3)%np,block(3)%x,block(3)%y,block(3)%z,block(3)%hull)      ! T02
   call laypts3(xC,yC,zC,xL,yL,zL,xF,yF,zF,block(4)%np,block(4)%x,block(4)%y,block(4)%z,block(4)%hull)      ! T03
   call laypts3(xB,yB,zB,xG,yG,zG,xI,yI,zI,block(5)%np,block(5)%x,block(5)%y,block(5)%z,block(5)%hull)      ! T04
   call laypts3(xB,yB,zB,xK,yK,zK,xH,yH,zH,block(6)%np,block(6)%x,block(6)%y,block(6)%z,block(6)%hull)      ! T05
   call laypts3(xD,yD,zD,xJ,yJ,zJ,xG,yG,zG,block(7)%np,block(7)%x,block(7)%y,block(7)%z,block(7)%hull)      ! T06
   call laypts3(xD,yD,zD,xH,yH,zH,xL,yL,zL,block(8)%np,block(8)%x,block(8)%y,block(8)%z,block(8)%hull)      ! T07
   call laypts3(xA,yA,zA,xK,yK,zK,xI,yI,zI,block(9)%np,block(9)%x,block(9)%y,block(9)%z,block(9)%hull)      ! T08
   call laypts3(xB,yB,zB,xI,yI,zI,xK,yK,zK,block(10)%np,block(10)%x,block(10)%y,block(10)%z,block(10)%hull) ! T09
   call laypts3(xC,yC,zC,xJ,yJ,zJ,xL,yL,zL,block(11)%np,block(11)%x,block(11)%y,block(11)%z,block(11)%hull) ! T10
   call laypts3(xD,yD,zD,xL,yL,zL,xJ,yJ,zJ,block(12)%np,block(12)%x,block(12)%y,block(12)%z,block(12)%hull) ! T11
   call laypts3(xE,yE,zE,xC,yC,zC,xA,yA,zA,block(13)%np,block(13)%x,block(13)%y,block(13)%z,block(13)%hull) ! T12
   call laypts3(xF,yF,zF,xA,yA,zA,xC,yC,zC,block(14)%np,block(14)%x,block(14)%y,block(14)%z,block(14)%hull) ! T13
   call laypts3(xG,yG,zG,xB,yB,zB,xD,yD,zD,block(15)%np,block(15)%x,block(15)%y,block(15)%z,block(15)%hull) ! T14
   call laypts3(xH,yH,zH,xD,yD,zD,xB,yB,zB,block(16)%np,block(16)%x,block(16)%y,block(16)%z,block(16)%hull) ! T15
   call laypts3(xI,yI,zI,xG,yG,zG,xE,yE,zE,block(17)%np,block(17)%x,block(17)%y,block(17)%z,block(17)%hull) ! T16
   call laypts3(xJ,yJ,zJ,xE,yE,zE,xG,yG,zG,block(18)%np,block(18)%x,block(18)%y,block(18)%z,block(18)%hull) ! T17
   call laypts3(xK,yK,zK,xF,yF,zF,xH,yH,zH,block(19)%np,block(19)%x,block(19)%y,block(19)%z,block(19)%hull) ! T18
   call laypts3(xL,yL,zL,xH,yH,zH,xF,yF,zF,block(20)%np,block(20)%x,block(20)%y,block(20)%z,block(20)%hull) ! T19

   radius=1
   do ib=1,nblock
      do ip=1,block(ib)%np
         call project_on_sphere(radius,block(ib)%x(ip),block(ib)%y(ip),block(ib)%z(ip),&
                                block(ib)%r(ip),block(ib)%theta(ip),block(ib)%phi(ip))
      end do
   end do


   !call laypts3b(thetaA,phiA,thetaI,phiI,thetaE,phiE,block(1)%np,block(1)%x,block(1)%y,block(1)%z,block(1)%hull)      ! T00
   !call laypts3b(thetaA,phiA,thetaF,phiF,thetaK,phiK,block(2)%np,block(2)%x,block(2)%y,block(2)%z,block(2)%hull)      ! T01
   !call laypts3b(thetaC,phiC,thetaE,phiE,thetaJ,phiJ,block(3)%np,block(3)%x,block(3)%y,block(3)%z,block(3)%hull)      ! T02
   !call laypts3b(thetaC,phiC,thetaL,phiL,thetaF,phiF,block(4)%np,block(4)%x,block(4)%y,block(4)%z,block(4)%hull)      ! T03
   !call laypts3b(thetaB,phiB,thetaG,phiG,thetaI,phiI,block(5)%np,block(5)%x,block(5)%y,block(5)%z,block(5)%hull)      ! T04
   !call laypts3b(thetaB,phiB,thetaK,phiK,thetaH,phiH,block(6)%np,block(6)%x,block(6)%y,block(6)%z,block(6)%hull)      ! T05
   !call laypts3b(thetaD,phiD,thetaJ,phiJ,thetaG,phiG,block(7)%np,block(7)%x,block(7)%y,block(7)%z,block(7)%hull)      ! T06
   !call laypts3b(thetaD,phiD,thetaH,phiH,thetaL,phiL,block(8)%np,block(8)%x,block(8)%y,block(8)%z,block(8)%hull)      ! T07
   !call laypts3b(thetaA,phiA,thetaK,phiK,thetaI,phiI,block(9)%np,block(9)%x,block(9)%y,block(9)%z,block(9)%hull)      ! T08
   !call laypts3b(thetaB,phiB,thetaI,phiI,thetaK,phiK,block(10)%np,block(10)%x,block(10)%y,block(10)%z,block(10)%hull) ! T09
   !call laypts3b(thetaC,phiC,thetaJ,phiJ,thetaL,phiL,block(11)%np,block(11)%x,block(11)%y,block(11)%z,block(11)%hull) ! T10
   !call laypts3b(thetaD,phiD,thetaL,phiL,thetaJ,phiJ,block(12)%np,block(12)%x,block(12)%y,block(12)%z,block(12)%hull) ! T11
   !call laypts3b(thetaE,phiE,thetaC,phiC,thetaA,phiA,block(13)%np,block(13)%x,block(13)%y,block(13)%z,block(13)%hull) ! T12
   !call laypts3b(thetaF,phiF,thetaA,phiA,thetaC,phiC,block(14)%np,block(14)%x,block(14)%y,block(14)%z,block(14)%hull) ! T13
   !call laypts3b(thetaG,phiG,thetaB,phiB,thetaD,phiD,block(15)%np,block(15)%x,block(15)%y,block(15)%z,block(15)%hull) ! T14
   !call laypts3b(thetaH,phiH,thetaD,phiD,thetaB,phiB,block(16)%np,block(16)%x,block(16)%y,block(16)%z,block(16)%hull) ! T15
   !call laypts3b(thetaI,phiI,thetaG,phiG,thetaE,phiE,block(17)%np,block(17)%x,block(17)%y,block(17)%z,block(17)%hull) ! T16
   !call laypts3b(thetaJ,phiJ,thetaE,phiE,thetaG,phiG,block(18)%np,block(18)%x,block(18)%y,block(18)%z,block(18)%hull) ! T17
   !call laypts3b(thetaK,phiK,thetaF,phiF,thetaH,phiH,block(19)%np,block(19)%x,block(19)%y,block(19)%z,block(19)%hull) ! T18
   !call laypts3b(thetaL,phiL,thetaH,phiH,thetaF,phiF,block(20)%np,block(20)%x,block(20)%y,block(20)%z,block(20)%hull) ! T19

case default 

   stop 'map_blocks: mtype value not supported'

end select

call cpu_time(t4) ; write(*,'(a,f6.3,a)') 'map_blocks:',t4-t3,'s'

end subroutine

!==================================================================================================!
!==================================================================================================!
