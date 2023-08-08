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
use SHTOOLS

contains

    module subroutine swiftest_sph_g_acc_one(gm, r_0, phi, theta, rh, c_lm, g_sph)
        !! author: Kaustub P. Anand
        !!
        !! Calculate the acceleration terms for one pair of bodies given c_lm, theta, phi, r
        !!

        implicit none
        ! Arguments
        real(DP), intent(in)        :: gm                       !! GMass of the central body
        real(DP), intent(in)        :: r_0                      !! radius of the central body
        real(DP), intent(in)        :: phi                      !! Azimuthal/Phase angle (radians)
        real(DP), intent(in)        :: theta                    !! Inclination/Zenith angle (radians)
        real(DP), intent(in), dimension(NDIM)       :: rh       !! distance vector of body
        real(DP), intent(in), dimension(2, :, :)    :: c_lm     !! Spherical Harmonic coefficients
        real(DP), intent(out), dimension(NDIM)      :: g_sph    !! acceleration vector

        ! Internals
        integer        :: l, m              !! SPH coefficients
        real(DP)       :: r_mag             !! magnitude of rh
        real(DP)       :: p                 !! Associated Lengendre Polynomials at a given cos(theta)
        integer        :: l_max             !! max Spherical Harmonic l order value

        r_mag = sqrt(dot_product(rh(:), rh(:)))
        l_max = size(c_lm, 2)
        PlBar(p, lmax, cos(theta))      ! Associated Legendre Polynomials

        do l = 1, l_max
            do m = 1, l

                ! Associated Legendre Polynomials   
                plm = p(PlmIndex(l, m))         ! p_l,m
                plm1 = p(PlmIndex(l, m+1))      ! p_l,m+1

                ! Normalization = 4*pi (geodesy) normalized
                N = sqrt((2 * l + 1) * gamma(l - m + 1) / gamma(l + m + 1))

                ! C_lm and S_lm with Cos and Sin of m * phi
                ccss = c_lm(1, l, m) * cos(m * phi) + c_lm(2, l, m) * sin(m * phi)          ! C_lm * cos(m * phi) + S_lm * sin(m * phi)
                cssc = -1.0 * c_lm(1, l, m) * sin(m * phi) + c_lm(2, l, m) * cos(m * phi)   ! - C_lm * sin(m * phi) + S_lm * cos(m * phi) 
                                                                                            ! cssc * m = first derivative of ccss with respect to phi

                ! m > 0
                g_sph(1) -= gm * r_0**l / r_mag**(l + 1) * N * (-1.0 * m * plm * cssc / rh(2) + ccss * (plm * (m * cotan(theta) / (rh(3) * cos(phi)) - (l + 1) * rh(1) / r_mag**2) + plm1 / (rh(3) * cos(phi)))) ! g_x
                g_sph(2) -= gm * r_0**l / r_mag**(l + 1) * N * (m * plm * cssc / rh(1) + ccss * (plm * (m * cotan(theta) / (rh(3) * sin(phi)) - (l + 1) * rh(1) / r_mag**2) + plm1 / (rh(3) * sin(phi)))) ! g_y
                g_sph(3) -= gm * r_0**l / r_mag**(l + 1) * N * (ccss * (plm * (m * cotan(theta) / sqrt(r_mag**2 - rh(3)**3) - (l + 1) * rh(1) / r_mag**2) + plm1 / sqrt(r_mag**2 - rh(3)**2))) ! g_z

            end do
        end do

        return
        end subroutine swiftest_sph_g_acc_one

    module subroutine swiftest_sph_g_acc_all()
        !! author: Kaustub P. Anand
        !!
        !! Calculate the acceleration terms for all bodies given c_lm
        !!

        associate(pl => self, npl => self%nbody, cb => nbody_system%cb)
            do concurrent(i = 1:npl, pl%lmask(i))
                gm = pl%Gmass(i)
                rh = pl%rh(:, i) ! CHECK pl%rh shape
                r = sqrt(rh(1)**2 + rh(2)**2 + rh(3)**2) ! mag of vector function???
                theta = acos(rh(3) / r)
                phi = acos(rh(1) / sqrt(rh(1)**2 + rh(2)**2)) * sign(1, rh(2)) - cb%phase ! CALCULATE CB PHASE VALUE FOR PHI

                call swiftest_sph_g_acc_one(gm, r, phi, theta, rh, cb%c_lm, g_sph)
                pl%ah(:, i) = pl%ah(:, i) + g_sph

        end subroutine swiftest_sph_g_acc_all
