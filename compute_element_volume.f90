!==================================================================================================!
!  ***  *  *  ***   **** ***** !                                                                   !
! *     *  * *   * *       *   !                                                                   !
! *  *  **** *   *  ***    *   !                                                                   !
! *   * *  * *   *     *   *   !                                                                   !
!  ***  *  *  ***  ****    *   !                                                       C. Thieulot !
!==================================================================================================!

function hexahedron_volume (x,y,z)

implicit none

real(8) x(8),y(8),z(8) 
real(8) hexahedron_volume
real(8), external :: triple_product 

!==================================================================================================!
!@@ \subsection{\tt hexahedron\_volume}
!@@ This function computes the volume of any hexahedron, 
!@@ following Eq.(12) of \cite{gran97}.
!==================================================================================================!

hexahedron_volume=( triple_product ( x(7)-x(2)+x(8)-x(1), y(7)-y(2)+y(8)-y(1), z(7)-z(2)+z(8)-z(1), &
                                     x(7)-x(4),           y(7)-y(4),           z(7)-z(4),           &
                                     x(3)-x(1),           y(3)-y(1),           z(3)-z(1)           )&           
                   +triple_product ( x(8)-x(1),           y(8)-y(1),           z(8)-z(1),           &
                                     x(7)-x(4)+x(6)-x(1), y(7)-y(4)+y(6)-y(1), z(7)-z(4)+z(6)-z(1), &
                                     x(7)-x(5),           y(7)-y(5),           z(7)-z(5)           )&           
                   +triple_product ( x(7)-x(2),           y(7)-y(2),           z(7)-z(2),           &
                                     x(6)-x(1),           y(6)-y(1),           z(6)-z(1),           &
                                     x(7)-x(5)+x(3)-x(1), y(7)-y(5)+y(3)-y(1), z(7)-z(5)+z(3)-z(1)  ) )/12.d0
end function

!==================================================================================================!
! In vector calculus, there are two ways of multiplying three vectors together, 
! to make a triple product of vectors. The scalar triple product is defined as 
! the dot product of one of the vectors with the cross product of the other two.
! The scalar triple product can also be understood as the determinant of the 
! 3-by-3 matrix having the three vectors as rows (or columns, since the determinant 
! for a transposed matrix, is the same as the original)
      
function triple_product (Ax,Ay,Az,Bx,By,Bz,Cx,Cy,Cz)
implicit none
real(8) triple_product
real(8) Ax,Ay,Az,Bx,By,Bz,Cx,Cy,Cz
      
triple_product = Ax * ( By * Cz - Bz * Cy ) &
               - Ay * ( Bx * Cz - Bz * Cx ) &
               + Az * ( Bx * Cy - By * Cx )
      
end function

!==================================================================================================!
!==================================================================================================!
