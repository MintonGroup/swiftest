! Copyright 2026 - The Minton Group at Purdue University
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
        module subroutine shgrav_g_acc_one(GMcb, r_0, phi_cb, rh, c_lm, g_sph, GMpl, aoblcb)
        !! author: Kaustub P. Anand
        !!
        !! Calculate the acceleration terms for one pair of bodies given c_lm, theta, phi, r
        implicit none
        ! Arguments
        real(DP), intent(in) :: GMcb 
            !! GMass of the central body
        real(DP), intent(in) :: r_0 
            !! radius of the central body
        real(DP), intent(in) :: phi_cb 
            !! rotation phase angle of the central body
        real(DP), intent(in), dimension(:) :: rh 
            !! distance vector of body
        real(DP), intent(in), dimension(:, :, :) :: c_lm 
            !! Spherical Harmonic coefficients
        real(DP), intent(out), dimension(NDIM) :: g_sph 
            !! acceleration vector
        real(DP), intent(in),  optional :: GMpl 
            !! Mass of input body if it is not a test particle
        real(DP), dimension(:), intent(inout), optional :: aoblcb
            !! Barycentric acceleration of central body (only for massive input b

        end subroutine shgrav_g_acc_one
        
        module subroutine shgrav_acc(body, nbody_system)
            implicit none
            class(swiftest_body), intent(inout) :: body
                !! Swiftest body object
            class(swiftest_nbody_system), intent(inout) :: nbody_system 
                !! Swiftest nbody system object
        end subroutine shgrav_acc

        module subroutine shgrav_pot_system(self)
            implicit none
            class(swiftest_nbody_system), intent(inout) :: self
                !! Swiftest nbody system object
        end subroutine shgrav_pot_system
        
    end interface

end module shgrav
 
