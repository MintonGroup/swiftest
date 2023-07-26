!! Copyright 2022 - David Minton, Carlisle Wishard, Jennifer Pouplin, Jake Elliott, & Dana Singh
!! This file is part of Swiftest.
!! Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
!! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!! Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
!! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!! You should have received a copy of the GNU General Public License along with Swiftest. 
!! If not, see: https://www.gnu.org/licenses. 


!! Swiftest submodule to calculate higher order terms for gravitational acceleration given spherical harmonic coefficients (c_lm)

submodule (swiftest) s_swiftest_sph
use swiftest
use operators

contains

    module subroutine swiftest_sph_g_acc_one(gm, r_0, phi, theta, rh, c_lm, g_sph)
        !! author: Kaustub P. Anand
        !!
        !! Calculate the acceleration terms for one pair of bodies given c_lm, theta, phi, r
        !!

        implicit none
        ! Arguments
        real(DP), intent(in)        :: gm !! GMass of the central body
        real(DP), intent(in)        :: r_0 !! radius of the central body
        real(DP), intent(in)        :: phi !! Azimuthal/Phase angle (radians)
        real(DP), intent(in)        :: theta !! Inclination/Zenith angle (radians)
        real(DP), intent(in), dimension(NDIM)       :: rh !! distance vector of body
        real(DP), intent(in), dimension(2, :, :)    :: c_lm !! Spherical Harmonic coefficients
        real(DP), intent(out), dimension(NDIM)      :: g_sph !! acceleration vector

        ! Internals
        integer        :: l, m !! SPH coefficients
        real(DP)       :: r_mag !! magnitude of rh

        r_mag = sqrt(dot_product(rh(:), rh(:)))

        ! WHAT DO I DO WITH THE COMPLEX TERM????

        do l = 1, l_max
            do m = 1, l

                swiftest_sph_dylm_dtheta(ylm, ylm1, theta, phi, l, m, dylm_dtheta, dyl_m_dtheta)
                ! m > 0
                g_sph(1) += gm * r_0**l / r_mag**(l + 1) * c_lm[0, l, m] * (dylm_dtheta / (rh(3) * cos(phi)) - (l + 1) * rh(1) * ylm / r_mag**2 - i * m * ylm / rh(2)) ! i = sqrt(-1)
                g_sph(2) += gm * r_0**l / r_mag**(l + 1) * c_lm[0, l, m] * (dylm_dtheta / (rh(3) * sin(phi)) - (l + 1) * rh(2) * ylm / r_mag**2 + i * m * ylm / rh(1)) 
                g_sph(3) += gm * r_0**l / r_mag**(l + 1) * c_lm[0, l, m] * (dylm_dtheta / sqrt(rh(1)**2 + rh(2)**2) - (l + 1) * rh(3) * ylm / r_mag**2) 

                ! m < 0
                g_sph(1) += gm * r_0**l / r_mag**(l + 1) * c_lm[1, l, m] * (dyl_m_dtheta / (rh(3) * cos(phi)) - (l + 1) * rh(1) * ylm * (-1)**m / r_mag**2 - i * m * ylm / rh(2)) ! i = sqrt(-1) 
                g_sph(2) += gm * r_0**l / r_mag**(l + 1) * c_lm[1, l, m] * (dyl_m_dtheta / (rh(3) * sin(phi)) - (l + 1) * rh(2) * ylm * (-1)**m / r_mag**2 + i * m * ylm / rh(1)) 
                g_sph(3) += gm * r_0**l / r_mag**(l + 1) * c_lm[1, l, m] * (dyl_m_dtheta / sqrt(rh(1)**2 + rh(2)**2) - (l + 1) * rh(3) * ylm * (-1)**m / r_mag**2) 

            end do
        end do

        return
        end subroutine swiftest_sph_g_acc_one

    module subroutine swiftest_sph_g_acc_all()
        !! author: Kaustub P. Anand
        !!
        !! Calculate the acceleration terms for all bodies given c_lm
        !!

        end subroutine swiftest_sph_g_acc_all

    module subroutine swiftest_sph_dylm_dtheta(ylm, ylm1, theta, phi, l, m, dylm_dtheta, dyl_m_dtheta)
        !! author: Kaustub P. Anand
        !!
        !! Calculate the derivative of Y_lm with respect to theta dY_lm / dtheta term
        !!

        ! Arguments
        real(DP), intent(in)        :: ylm ! Y_l,m
        real(DP), intent(in)        :: ylm1 ! Y_l,m+1
        real(DP), intent(in)        :: theta ! Zenith angle
        real(DP), intent(in)        :: phi ! Azimuthal angle
        integer, intent(in)         :: l, m ! spherical harmonics numbers
        real(DP), intent(out)       :: dylm_dtheta ! derivative for +m
        real(DP), intent(out)       :: dyl_m_dtheta ! derivative for -m


        ! Internals

        dylm_dtheta = m * cotan(theta) * ylm + sqrt((l - m) * (l + m + 1)) * ylm1 * exp(i * phi)
        dyl_m_dtheta = -m * cotan(theta) * ylm * (-1)**m + sqrt((l + m) * (l - m + 1)) * ylm1 * exp(i * phi)

        return
        end subroutine swiftest_sph_dylm_dtheta

    module subroutine swiftest_sph_ylm(l, m, phi, theta)
        !! author: Kaustub P. Anand
        !!
        !! Calculate the Y_lm for a given l, m, phi, and theta
        !!

        integer, intent(in)         :: l, m ! spherical harmonics numbers
        real(DP), intent(in)        :: theta ! Zenith angle
        real(DP), intent(in)        :: phi ! Azimuthal angle

        return
        end subroutine swiftest_sph_ylm