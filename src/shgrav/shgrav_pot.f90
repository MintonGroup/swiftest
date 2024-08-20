!! Copyright 2024 - The Minton Group at Purdue University
!! This file is part of Swiftest.
!! Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
!! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!! Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
!! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!! You should have received a copy of the GNU General Public License along with Swiftest. 
!! If not, see: https://www.gnu.org/licenses. 


!! Swiftest submodule to calculate higher order terms for gravitational potential given spherical harmonic coefficients (c_lm)

submodule (shgrav) s_shgrav_pot
use swiftest
use SHTOOLS

contains

    module subroutine shgrav_pot_system(self)
        !! author: Kaustub P. Anand
        !!
        !! Compute the gravitational potential energy for a given system of bodies solely due to the gravitational harmonics terms
        !!
        implicit none
        ! Arguments
        class(swiftest_nbody_system), intent(inout) :: self
            !! Swiftest nbody system object
        ! Internals
        integer(I4B) :: i, npl
        real(DP), dimension(self%pl%nbody) :: oblpot_arr

        associate(nbody_system => self, pl => self%pl, cb => self%cb)
            npl = self%pl%nbody

            do i = 1, npl
                if (pl%lmask(i)) then
                    oblpot_arr(i) = shgrav_pot_one(cb%Gmass,  pl%mass(i), cb%radius, cb%rotphase*DEG2RAD, pl%rh(:,i), cb%c_lm) 
                endif
            end do

            nbody_system%oblpot = sum(oblpot_arr, pl%lmask(1:npl)) 
        end associate

        return
    end subroutine shgrav_pot_system

    function shgrav_pot_one(GMcb, Mpl, r_0, phi_cb, rh, c_lm) result(oblpot)
        !! author: Kaustub P. Anand
        !!
        !! Calculate the gravitational potential energy for one pair of bodies given c_lm, theta, phi, r
        !! Similar to shgrav_g_acc_one and swiftest_obl_pot_one
        !!
        implicit none
        ! Arguments
        real(DP), intent(in) :: GMcb 
            !! GMass of the central body
        real(DP), intent(in), optional :: Mpl
            !! Mass of input body
        real(DP), intent(in) :: r_0
            !! radius of the central body
        real(DP), intent(in) :: phi_cb
            !! rotation phase angle of the central body
        real(DP), intent(in), dimension(:) :: rh
            !! distance vector of body
        real(DP), intent(in), dimension(:, :, :) :: c_lm
            !! Spherical Harmonic coefficients

        ! Result
        real(DP) :: oblpot
            !! Gravitational potential energy

        ! Internals
        integer :: l, m 
            !! SPH coefficients
        integer :: l_max 
            !! max Spherical Harmonic l order value
        integer(I4B) :: N, lmindex 
            !! Length of Legendre polynomials and index at a given l, m
        real(DP) :: r_mag  
            !! magnitude of rh
        real(DP) :: phi, phi_bar 
            !! Azimuthal/Phase angle (radians) wrt coordinate axes, and central body rotation phase
        real(DP) :: theta 
            !! Inclination/Zenith angle (radians)
        real(DP) :: plm, plm1 
            !! Associated Legendre polynomials at a given l, m
        real(DP) :: ccss
            !! See definition in source code
        real(DP) :: cos_theta, sin_theta 
            !! cos(theta) and sin(theta)
        real(DP), dimension(:), allocatable  :: p 
            !! Associated Lengendre Polynomials at a given cos(theta)
        real(DP) :: r_fac 
            !! calculation factors

        theta = atan2(sqrt(rh(1)**2 + rh(2)**2), rh(3))
        phi = atan2(rh(2), rh(1)) 
        phi_bar = MOD(phi - phi_cb, 2 * PI) ! represents the phase difference between the central body's and the particle's phase
        r_mag = sqrt(dot_product(rh(:), rh(:)))
        l_max = size(c_lm, 2) - 1
        N = (l_max + 1) * (l_max + 2) / 2
        allocate(p(N))

        cos_theta = cos(theta)
        sin_theta = sin(theta)

        ! check if cos_theta is too small to avoid floating underflow error
        if (abs(cos_theta) < EPSILON(0.0_DP)) then
            call PlmBar(p, l_max, 0.0_DP)
        else
            call PlmBar(p, l_max, cos_theta)
        end if

        oblpot = 0.0_DP
        do l = 1, l_max ! skipping the l = 0 term; It is the spherical body term
            do m = 0, l

                ! If c_lm is too small, skip the iteration to improve performance
                if (abs(c_lm(m+1, l+1, 1)) < epsilon(0.0_DP) .and. abs(c_lm(m+1, l+1, 2)) < epsilon(0.0_DP)) then
                    cycle  
                endif

                ! Associated Legendre Polynomials 
                lmindex = PlmIndex(l, m)  
                plm = p(lmindex)                  ! p_l,m

                r_fac = GMcb * Mpl * r_0**l / r_mag**(l + 1)
                ccss = cos(m * phi_bar) * c_lm(m+1, l+1, 1) &
                        + sin(m * phi_bar) * c_lm(m+1, l+1, 2) ! C_lm * cos(m * phi_bar) + S_lm * sin(m * phi_bar)
                
                oblpot = oblpot + r_fac * plm * ccss
            end do
        end do

        ! deallocate(p)

        return
    end function shgrav_pot_one

end submodule s_shgrav_pot