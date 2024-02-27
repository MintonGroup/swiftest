! Copyight 2024 - The Minton Group at Purdue University
! This file is part of Swiftest.
! Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
! Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
! You should have received a copy of the GNU General Public License along with Swiftest. 
! If not, see: https://www.gnu.org/licenses. 

module shgrav
    !! author: David A. Minton and Kaustub Anand
    !!
    !! This module defines functions used for the computation of accelerations based on spherical harmonics representation of the
    !! gravitational potential of the central body. It uses the SHTOOLS library https://shtools.github.io/SHTOOLS/
    use swiftest
    implicit none
    public

    interface
        module subroutine shgrav_acc(body, nbody_system)
            implicit none
            class(swiftest_body), intent(inout) :: body
                !! Swiftest body object
            class(swiftest_nbody_system), intent(inout) :: nbody_system 
                !! Swiftest nbody system object
        end subroutine shgrav_acc
        
    end interface

end module shgrav
 