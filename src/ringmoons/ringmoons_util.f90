! Copyright 2025 - The Minton Group at Purdue University
! This file is part of Swiftest.
! Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
! Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
! You should have received a copy of the GNU General Public License along with Swiftest. 
! If not, see: https://www.gnu.org/licenses. 

submodule(ringmoons) s_ringmoons_util
   use swiftest
contains
   
    module subroutine ringmoons_util_dealloc_ring(self)
        !! author: David A. Minton
        !!
        !! Deallocates all allocatabale arrays
        implicit none
        ! Arguments
        class(ringmoons_ring),  intent(inout) :: self 
        !! Ringmoons ring object

        return
    end subroutine ringmoons_util_dealloc_ring
end submodule s_ringmoons_util