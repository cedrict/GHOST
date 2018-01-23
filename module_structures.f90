!==================================================================================================!
!  ***  *  *  ***   **** ***** !                                                                   !
! *     *  * *   * *       *   !                                                                   !
! *  *  **** *   *  ***    *   !                                                                   !
! *   * *  * *   *     *   *   !                                                                   !
!  ***  *  *  ***  ****    *   !                                                       C. Thieulot !
!==================================================================================================!

module structures

real(8), parameter :: Ggrav =6.6738480d-11
real(8), parameter :: pi  = 3.14159265358979323846264338327950288d0
real(8), parameter :: pi4  = 3.14159265358979323846264338327950288d0*0.25d0
real(8), parameter :: sqrt2 = 1.414213562373095048801688724209d0
real(8), parameter :: sqrt3 = 1.732050807568877293527446341505d0
real(8), dimension(2), parameter :: qcoords2=(/-1.d0/sqrt3,+1.d0/sqrt3/)

real(8) t1,t2,t3,t4
real(8) R1,R2
real(8) Lx,Ly,Lz
real(8) NNN(8,8),dNNNdr(8,8),dNNNds(8,8),dNNNdt(8,8)
real(8) rho
real(8) weight

integer :: ncellx,ncelly,ncellz
integer :: nnx,nny,nnz
integer :: mtype ! mesh type
integer :: nblock
integer :: level
integer :: nlayer! nb of shell layers in radial direction
integer :: ndim
integer, parameter :: iunit=66

logical generate_vtu_output
logical generate_ascii_output
logical compute_gravity 
logical equiangular
logical s40rts

type mesh 
   integer ncell ! nb of cells
   integer np    ! nb of points 
   integer nv    ! nb of vertices per cell
   integer nq    ! nb of Gauss integration points in the mesh
   integer nqel  ! nb of Gauss integration points per element
   integer, dimension(:,:), allocatable :: icon      ! connectivity array
   integer, dimension(:), allocatable :: blocknumber ! block number
   real(8), dimension(:), allocatable :: x           ! x coordinates
   real(8), dimension(:), allocatable :: y           ! y coordinates
   real(8), dimension(:), allocatable :: z           ! z coordinates
   real(8), dimension(:), allocatable :: r           ! r coordinates
   real(8), dimension(:), allocatable :: theta       ! theta coordinates
   real(8), dimension(:), allocatable :: phi         ! phi coordinates
   real(8), dimension(:), allocatable :: area        ! areas of shell cell/blocks
   real(8), dimension(:), allocatable :: volume      ! volumes of cells 
   real(8), dimension(:), allocatable :: xq          ! x coordinates of quad. points
   real(8), dimension(:), allocatable :: yq          ! y coordinates of quad. points
   real(8), dimension(:), allocatable :: zq          ! z coordinates of quad. points
   real(8), dimension(:), allocatable :: nx          ! x component of normal vectors
   real(8), dimension(:), allocatable :: ny          ! y component of normal vectors
   real(8), dimension(:), allocatable :: nz          ! z component of normal vectors
   real(8), dimension(:), allocatable :: dlnvs       ! velocity of s40rts 
   real(8), dimension(:), allocatable :: drho        ! density anomaly  
   logical, dimension(:), allocatable :: inner_node  !       
   logical, dimension(:), allocatable :: outer_node  !
   logical, dimension(:), allocatable :: inner_cell  !
   logical, dimension(:), allocatable :: outer_cell  !
   logical, dimension(:), allocatable :: hull        ! 
end type   

type(mesh), dimension(:), allocatable :: block
type(mesh) :: shell
type(mesh) :: shell_temp
type(mesh) :: hollow_sphere

end module

!==================================================================================================!
