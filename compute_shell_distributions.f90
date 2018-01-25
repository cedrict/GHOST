!==================================================================================================!
!  ***  *  *  ***   **** ***** !                                                                   !
! *     *  * *   * *       *   !                                                                   !
! *  *  **** *   *  ***    *   !                                                                   !
! *   * *  * *   *     *   *   !                                                                   !
!  ***  *  *  ***  ****    *   !                                                       C. Thieulot !
!==================================================================================================!

subroutine compute_shell_distributions

use structures

implicit none

integer ic,i
integer, parameter :: nbin=10
integer area_distr(0:nbin)
integer nedge,iedge
real(8) area_min,area_max,area_avrg
real(8) x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4
real(8) dist1,dist2,dist3,dist4,edge_avrg,edge_stddev
real(8), external :: distance
real(8), dimension(:), allocatable :: dist_edges

!==================================================================================================!
!@@ \subsection{\tt compute\_shell\_distributions}
!@@ This routine computes the average edge length of a shell and its standard deviation.
!@@ Results are generated in the {\sl fort.889} file.
!==================================================================================================!

area_min=minval(shell%area)*0.99
area_max=maxval(shell%area)*1.01

write(iunit,*) 'area_min',area_min
write(iunit,*) 'area_max',area_max

area_avrg=sum(shell%area)/shell%ncell

write(iunit,*) 'area_avrg',area_avrg
      
area_distr=0
do ic=1,shell%ncell
   i=int((shell%area(ic)-area_min)/(area_max-area_min)*nbin)
   area_distr(i)=area_distr(i)+1
end do

open(unit=123,file='OUT/shell_area_distribution.dat')
do i=1,nbin
   if (area_distr(i)/=0) write(123,*) ((i-1.)/(nbin-1.)*(area_max-area_min)+area_min)/area_avrg,area_distr(i)
end do
close(123)

!-------------------------------
! all edges are counted twice since they all 
! share an element. 
 
nedge=shell%nv*shell%ncell 
allocate(dist_edges(nedge))

iedge=0
do ic=1,shell%ncell
   select case(mtype)
   case(06,12)
      x1=shell%x(shell%icon(1,ic)) ; y1=shell%y(shell%icon(1,ic)) ; z1=shell%z(shell%icon(1,ic)) 
      x2=shell%x(shell%icon(2,ic)) ; y2=shell%y(shell%icon(2,ic)) ; z2=shell%z(shell%icon(2,ic)) 
      x3=shell%x(shell%icon(3,ic)) ; y3=shell%y(shell%icon(3,ic)) ; z3=shell%z(shell%icon(3,ic)) 
      x4=shell%x(shell%icon(4,ic)) ; y4=shell%y(shell%icon(4,ic)) ; z4=shell%z(shell%icon(4,ic)) 
      dist1=distance(x1,y1,z1,x2,y2,z2) ; iedge=iedge+1 ; dist_edges(iedge)=dist1
      dist2=distance(x2,y2,z2,x3,y3,z3) ; iedge=iedge+1 ; dist_edges(iedge)=dist2
      dist3=distance(x3,y3,z3,x4,y4,z4) ; iedge=iedge+1 ; dist_edges(iedge)=dist3
      dist4=distance(x4,y4,z4,x1,y1,z1) ; iedge=iedge+1 ; dist_edges(iedge)=dist4
   case(20)
      x1=shell%x(shell%icon(1,ic)) ; y1=shell%y(shell%icon(1,ic)) ; z1=shell%z(shell%icon(1,ic)) 
      x2=shell%x(shell%icon(2,ic)) ; y2=shell%y(shell%icon(2,ic)) ; z2=shell%z(shell%icon(2,ic)) 
      x3=shell%x(shell%icon(3,ic)) ; y3=shell%y(shell%icon(3,ic)) ; z3=shell%z(shell%icon(3,ic)) 
      dist1=distance(x1,y1,z1,x2,y2,z2) ; iedge=iedge+1 ; dist_edges(iedge)=dist1
      dist2=distance(x2,y2,z2,x3,y3,z3) ; iedge=iedge+1 ; dist_edges(iedge)=dist2
      dist3=distance(x3,y3,z3,x1,y1,z1) ; iedge=iedge+1 ; dist_edges(iedge)=dist3
      nedge=nedge+3
   case default
      stop 'compute_shell_distributions: mtype not supported'
   end select

end do

edge_avrg=sum(dist_edges)/nedge
edge_stddev=sqrt( sum((dist_edges-edge_avrg)**2)/nedge  )

write(iunit,*) 'shell edge (avrg)',edge_avrg
write(iunit,*) 'shell edge (stddev)',edge_stddev

write(889,*) shell%np,&
             edge_avrg,&
             edge_stddev ; call flush(889)

end subroutine

!==================================================================================================!
!==================================================================================================!
