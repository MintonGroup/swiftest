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
use operators
use SHTOOLS

contains

    module subroutine swiftest_sph_g_acc_one(GMcb, r_0, phi, theta, rh, c_lm, g_sph, GMpl, aoblcb)
        !! author: Kaustub P. Anand
        !!
        !! Calculate the acceleration terms for one pair of bodies given c_lm, theta, phi, r
        !!

        implicit none
        ! Arguments
        real(DP), intent(in)        :: GMcb                        !! GMass of the central body
        real(DP), intent(in)        :: r_0                         !! radius of the central body
        real(DP), intent(in)        :: phi                         !! Azimuthal/Phase angle (radians)
        real(DP), intent(in)        :: theta                       !! Inclination/Zenith angle (radians)
        real(DP), intent(in), dimension(:)          :: rh          !! distance vector of body
        real(DP), intent(in), dimension(:, :, :)    :: c_lm        !! Spherical Harmonic coefficients
        real(DP), intent(out), dimension(NDIM)         :: g_sph       !! acceleration vector
        real(DP), intent(in),  optional :: GMpl                    !! Mass of input body if it is not a test particle
        real(DP), dimension(:),   intent(inout), optional :: aoblcb  !! Barycentric acceleration of central body (only for massive input bodies)
     
        ! Internals
        integer        :: l, m              !! SPH coefficients
        integer        :: l_max             !! max Spherical Harmonic l order value
        integer(I4B)   :: N, lmindex        !! Length of Legendre polynomials and index at a given l, m
        real(DP)       :: r_mag             !! magnitude of rh
        real(DP)       :: plm, dplm         !! Associated Legendre polynomials at a given l, m
        real(DP)       :: ccss, cssc        !! See definition in source code
        real(DP), dimension(:), allocatable  :: p, p_deriv   !! Associated Lengendre Polynomials at a given cos(theta)

        g_sph(:) = 0.0_DP
        r_mag = sqrt(dot_product(rh(:), rh(:)))
        l_max = size(c_lm, 2) - 1
        N = (l_max + 1) * (l_max + 2) / 2
        allocate(p(N),p_deriv(N))
        call PlmON_d1(p, p_deriv, l_max, cos(theta))      ! Orthonormalized Associated Legendre Polynomials and the 1st Derivative

        do l = 1, l_max ! skipping the l = 0 term; It is the spherical body term
            do m = 0, l

                ! Associated Legendre Polynomials 
                lmindex = PlmIndex(l, m)  
                plm = p(lmindex)              ! p_l,m
                dplm = p_deriv(lmindex)       ! d(p_l,m)

                ! C_lm and S_lm with Cos and Sin of m * phi
                ccss = c_lm(m+1, l+1, 1) * cos(m * phi) & 
                        + c_lm(m+1, l+1, 2) * sin(m * phi)      ! C_lm * cos(m * phi) + S_lm * sin(m * phi)
                cssc = -1.0_DP * c_lm(m+1, l+1, 1) * sin(m * phi) & 
                        + c_lm(m+1, l+1, 2) * cos(m * phi)      ! - C_lm * sin(m * phi) + S_lm * cos(m * phi) 
                                                                ! cssc * m = first derivative of ccss with respect to phi

                ! m > 0
                g_sph(1) = g_sph(1) - GMcb * r_0**l / r_mag**(l + 1) * (-1.0_DP * cssc * m * plm / rh(2) &
                                                              - ccss * (dplm * sin(theta) / (rh(3) * cos(phi)) &    
                                                                        + plm * (l + 1) * rh(1) / r_mag**2)) ! g_x
                g_sph(2) = g_sph(2) - GMcb * r_0**l / r_mag**(l + 1) * (cssc * m * plm / rh(1) &
                                                              - ccss * (dplm * sin(theta) / (rh(3) * sin(phi)) &
                                                                        + plm * (l + 1) * rh(2) / r_mag**2)) ! g_y
                g_sph(3) = g_sph(3) + GMcb * r_0**l / r_mag**(l + 1) * ccss * (dplm * sin(theta) / sqrt(r_mag**2 - rh(3)**2) &
                                                                     + plm * (l + 1) * rh(3) / r_mag**2) ! g_z
            end do
        end do

        if (present(GMpl) .and. present(aoblcb)) then
            aoblcb(:) = aoblcb(:) - GMpl * g_sph(:) / GMcb
        end if

        return
    end subroutine swiftest_sph_g_acc_one

    module subroutine swiftest_sph_g_acc_pl_all(self, nbody_system)
        !! author: Kaustub P. Anand
        !!
        !! Calculate the acceleration terms for all massive bodies given c_lm
        !!
        implicit none
        ! Arguments
        class(swiftest_pl),           intent(inout) :: self   !! Swiftest massive body object
        class(swiftest_nbody_system), intent(inout) :: nbody_system !! Swiftest nbody system object
        ! Internals
        integer(I4B)    :: i = 1
        real(DP)        :: r_mag, theta, phi    !! magnitude of the position vector, zenith angle, and azimuthal angle
        real(DP), dimension(NDIM)  :: g_sph        !! Gravitational terms from Spherical Harmonics

        associate(pl => self, npl => self%nbody, cb => nbody_system%cb, rh => self%rh)
            cb%aobl(:) = 0.0_DP

            do i = 1, npl
                if (pl%lmask(i)) then
                    r_mag = .mag. rh(:,i)
                    theta = atan2(sqrt(rh(1,i)**2 + rh(2,i)**2), rh(3,i))
                    phi = atan2(rh(2,i), rh(1,i)) - cb%rotphase ! CALCULATE CB PHASE VALUE FOR PHI
    
                    call swiftest_sph_g_acc_one(cb%Gmass, r_mag, phi, theta, rh(:,i), cb%c_lm, g_sph, pl%Gmass(i), cb%aobl)
                    pl%ah(:, i) = pl%ah(:, i) + g_sph(:) - cb%aobl(:)
                    pl%aobl(:, i) = g_sph(:)
                end if
            end do
        end associate
        return 
        end subroutine swiftest_sph_g_acc_pl_all
    
    module subroutine swiftest_sph_g_acc_tp_all(self, nbody_system)
        !! author: Kaustub P. Anand
        !!
        !! Calculate the acceleration terms for all test particles given c_lm
        !!
        implicit none
        ! Arguments
        class(swiftest_tp),           intent(inout) :: self   !! Swiftest test particle object
        class(swiftest_nbody_system), intent(inout) :: nbody_system !! Swiftest nbody system object
        ! Internals
        integer(I4B)    :: i = 1
        real(DP)        :: r_mag, theta, phi       !! magnitude of the position vector, zenith angle, and azimuthal angle
        real(DP), dimension(NDIM)  :: rh           !! Position vector of the test particle
        real(DP), dimension(NDIM)  :: g_sph        !! Gravitational terms from Spherical Harmonics
        real(DP), dimension(NDIM)  :: aoblcb       !! Temporary variable for central body oblateness acceleration

        associate(tp => self, ntp => self%nbody, cb => nbody_system%cb, rh => self%rh)

            if (nbody_system%lbeg) then
                aoblcb = cb%aoblbeg
             else
                aoblcb = cb%aoblend
             end if

            do i = 1, ntp
                if (tp%lmask(i)) then
                    r_mag = .mag. rh(:,i)
                    theta = atan2(sqrt(rh(1,i)**2 + rh(2,i)**2), rh(3,i))
                    phi = atan2(rh(2,i), rh(1,i)) - cb%rotphase

                    call swiftest_sph_g_acc_one(cb%Gmass, r_mag, phi, theta, rh(:,i), cb%c_lm, g_sph)
                    tp%ah(:, i) = tp%ah(:, i) + g_sph(:) - aoblcb(:)
                    tp%aobl(:, i) = g_sph(:)
                end if
            end do
        end associate
        return
        end subroutine swiftest_sph_g_acc_tp_all
    
end submodule s_swiftest_sph