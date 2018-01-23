!==================================================================================================!
!  ***  *  *  ***   **** ***** !                                                                   !
! *     *  * *   * *       *   !                                                                   !
! *  *  **** *   *  ***    *   !                                                                   !
! *   * *  * *   *     *   *   !                                                                   !
!  ***  *  *  ***  ****    *   !                                                       C. Thieulot !
!==================================================================================================!

program GHOST

use structures

implicit none

!======================================!
! default values
!======================================!

mtype=06
nlayer=8
level=8
rho=1.d3
equiangular=.true.
s40rts=.true.

!======================================!

call read_command_line_options   

!======================================!
! MAIN PARAMETERS
!======================================!
! mtype = 06 3: cubed sphere 
! mtype = 12 4: citcom mesh
! mtype = 20 5: icosahedron
!======================================!

select case(mtype)
case(06)
   nblock=6
   ndim=3
   R1=1
   R2=2
case(12)
   nblock=12
   ndim=3
   R1=1
   R2=2
case(20)
   nblock=20
   ndim=3
   R1=1
   R2=2
case default
   stop 'unknown mtype'
end select

open(unit=iunit,file='Log.out')

generate_vtu_output=.true.
generate_ascii_output=.false.
compute_gravity=.false.

!======================================!

call cpu_time(t1)

call check_parameters
call allocate_block_memory
call block_node_layout
call map_blocks
call output_blocks_vtu
call output_blocks_ascii
call merge_blocks
call compute_r_theta_phi_shell
call compute_shell_area
call compute_shell_distributions
call output_shell_vtu
call output_shell_ascii
call output_shell_lat_lon
call build_hollow_sphere
call compute_r_theta_phi_hollow_sphere
call compute_normals_hollow_sphere
call read_s40rts
call generate_quadrature_points
call compute_gravity_on_line 
call compute_hollow_sphere_volume
call output_hollow_sphere_vtu
!call output_hollow_sphere_qpts_vtu

call cpu_time(t2) 

write(999,*) hollow_sphere%np,t2-t1,hollow_sphere%ncell,level ; call flush(999) 
write(555,*) level,hollow_sphere%np,hollow_sphere%ncell,nblock ; call flush(555)

close(iunit)

end program

!==================================================================================================!
!==================================================================================================!
