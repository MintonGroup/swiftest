!! Copyright 2024 - The Minton Group at Purdue University
!! This file is part of Swiftest.
!! Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
!! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!! Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
!! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!! You should have received a copy of the GNU General Public License along with Swiftest. 
!! If not, see: https://www.gnu.org/licenses. 


!! Swiftest submodule to calculate radiation effects on massive bodies

module radiation
    !! author: Kaustub P Anand and David A. Minton
    !!
    !! This module defines functions used for the computation of radiation effects on massive bodies.
    !! Equations taken from Burns, Lamy & Soter (1979) Icarus 40, 1-48.

    use swiftest
    implicit none
    public

    interface radiation
        module subroutine radiation_getacch_pl(self, nbody_system, param)
            implicit none
            ! Arguments
        class(swiftest_body),         intent(inout) :: self
            !! Swiftest body object
        class(swiftest_nbody_system), intent(inout) :: nbody_system
            !! Swiftest nbody system object
        class(swiftest_parameters),   intent(in)    :: param
            !! Current run configuration parameters

        end subroutine radiation_getacch_pl

    end interface



end module radiation