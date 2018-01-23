.SUFFIXES:.out .o .s .c .F .f .f90 .e .r .y .yr .ye .l .p .sh .csh .h

F90=gfortran 
#FLAGS= -c -O3 -ffree-line-length-none
FLAGS= -c -ffree-line-length-none -Wall -Wextra -fPIC -fmax-errors=1 -g -fcheck=all -fbacktrace -fbounds-check -ffpe-trap= -Wno-error=implicit-interface -Wno-error=conversion

OBJECTS2D =\
module_structures.o\
module_S40RTS.o\
allocate_block_memory.o\
block_node_layout.o\
build_hollow_sphere.o\
check_parameters.o\
compute_element_volume.o\
compute_r_theta_phi_shell.o\
compute_r_theta_phi_hollow_sphere.o\
compute_hollow_sphere_volume.o\
compute_triangle_area.o\
compute_shell_area.o\
compute_shell_distributions.o\
compute_gravity_at_point.o\
compute_gravity_on_line.o\
compute_normals.o\
compute_PREM_density.o\
generate_quadrature_points.o\
int_to_char.o\
laypts3.o\
laypts4.o\
map_blocks.o\
merge_blocks.o\
output_blocks_vtu.o\
output_blocks_ascii.o\
output_shell_ascii.o\
output_shell_vtu.o\
output_shell_lat_lon.o\
output_hollow_sphere_vtu.o\
project_on_sphere.o\
read_command_line_options.o\
read_s40rts.o S40RTS.o\
vtu_tools.o\
ghost.o 

.f.o:
	$(F90) $(FLAGS) $(INCLUDE) $*.f

.f90.o:
	$(F90) $(FLAGS) $(INCLUDE) $*.f90

code:	$(OBJECTS2D)
	$(F90) $(OPTIONS) $(OBJECTS2D) $(LIBS) -o ghost

clean: 
	rm -f *.o
	rm -f *.dat
	rm -f *.out
	rm -f fort.*
	rm -f OUT/*.dat
	rm -f OUT/*.vtu
	rm -f OUT/ASCII/*.ascii
	rm ghost


