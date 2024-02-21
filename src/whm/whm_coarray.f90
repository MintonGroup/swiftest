! Copyight 2023 - David Minton
! This file is part of Swiftest.
! Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
! Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
! You should have received a copy of the GNU General Public License along with Swiftest. 
! If not, see: https://www.gnu.org/licenses. 

submodule (whm) s_whm_coarray
use coarray
use swiftest
contains

    module subroutine whm_coarray_coclone_pl(self)
        !! author: David A. Minton
        !!
        !! Broadcasts the image 1 object to all other images in a coarray 
        implicit none
        ! Arguments
        class(whm_pl),intent(inout),codimension[*]  :: self  !! WHM pl object

        call coclone(self%eta)
        call coclone(self%xj)
        call coclone(self%vj)
        call coclone(self%muj)
        call coclone(self%ir3j)

        call swiftest_coarray_coclone_pl(self)

        return
    end subroutine whm_coarray_coclone_pl

end submodule s_whm_coarray