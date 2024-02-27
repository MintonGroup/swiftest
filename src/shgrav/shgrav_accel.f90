!! Copyright 2022 - David Minton, Carlisle Wishard, Jennifer Pouplin, Jake Elliott, & Dana Singh
!! This file is part of Swiftest.
!! Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
!! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!! Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
!! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!! You should have received a copy of the GNU General Public License along with Swiftest. 
!! If not, see: https://www.gnu.org/licenses. 


!! Swiftest submodule to calculate higher order terms for gravitational acceleration given spherical harmonic coefficients (c_lm)

submodule (shgrav) s_shgrav_accel
use swiftest
use SHTOOLS

contains

    subroutine shgrav_g_acc_one(GMcb, r_0, phi_cb, rh, c_lm, g_sph, GMpl, aoblcb)
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
            !! Barycentric acceleration of central body (only for massive input bodies)
     
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
        real(DP) :: ccss, cssc 
            !! See definition in source code
        real(DP) :: cos_theta, sin_theta 
            !! cos(theta) and sin(theta)
        real(DP), dimension(:), allocatable  :: p 
            !! Associated Lengendre Polynomials at a given cos(theta)
        real(DP) :: fac1, fac2, r_fac 
            !! calculation factors 

        g_sph(:) = 0.0_DP
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

        do l = 1, l_max ! skipping the l = 0 term; It is the spherical body term
            do m = 0, l

                ! If c_lm is too small, skip the iteration to improve performance
                if (abs(c_lm(m+1, l+1, 1)) < epsilon(0.0_DP) .and. abs(c_lm(m+1, l+1, 2)) < epsilon(0.0_DP)) then
                    cycle  
                endif

                ! Associated Legendre Polynomials 
                lmindex = PlmIndex(l, m)  
                plm = p(lmindex)                  ! p_l,m

                ! C_lm and S_lm with Cos and Sin of m * phi
                ccss = c_lm(m+1, l+1, 1) * cos(m * phi_bar) & 
                        + c_lm(m+1, l+1, 2) * sin(m * phi_bar)      ! C_lm * cos(m * phi_bar) + S_lm * sin(m * phi_bar)
                cssc = -1 * c_lm(m+1, l+1, 1) * sin(m * phi_bar) & 
                        + c_lm(m+1, l+1, 2) * cos(m * phi_bar)      ! - C_lm * sin(m * phi_bar) + S_lm * cos(m * phi_bar) 
                                                                    ! cssc * m = first derivative of ccss with respect to phi

                if ((m+1) <= l) then
                    lmindex = PlmIndex(l, m+1) 
                    plm1 = p(lmindex) 
                    if(m == 0) then
                        plm1 = plm1 * sqrt(((l + m + 1) * (l - m)) / 2.0) ! renormalize plm1 to the norm of plm
                    else 
                        plm1 = plm1 * sqrt((l + m + 1) * (l - m) * 1.0)       ! renormalize plm1 to the norm of plm
                    end if
                else
                    plm1 = 0.0_DP  
                end if 
                                                                      
                if(abs(sin_theta) < epsilon(1.0_DP)) then
                    fac1 = 0.0_DP
                else
                    fac1 = m * plm / sin_theta
                end if

                fac2 = plm * (l + m + 1) * sin_theta + plm1 * cos_theta
                r_fac = -GMcb * r_0**l / r_mag**(l + 2)

                g_sph(1) = g_sph(1) + r_fac * (cssc * fac1 * sin(phi) + ccss * (fac2 - fac1) * cos(phi))
                g_sph(2) = g_sph(2) + r_fac * (-cssc * fac1 * cos(phi) + ccss * (fac2 - fac1) * sin(phi))
                g_sph(3) = g_sph(3) + r_fac * ccss * (plm * (l + m + 1) * cos_theta - plm1 * sin_theta)
          
            end do
        end do

        if (present(GMpl) .and. present(aoblcb)) then
            aoblcb(:) = aoblcb(:) - GMpl * g_sph(:) / GMcb
        end if

        return
    end subroutine shgrav_g_acc_one

    module subroutine shgrav_acc(body, nbody_system)
        !! author: Kaustub P. Anand
        !!
        !! Calculate the acceleration terms for bodies given c_lm values for the central body
        !!
        implicit none
        ! Arguments
        class(swiftest_body), intent(inout) :: body
            !! Swiftest body object
        class(swiftest_nbody_system), intent(inout) :: nbody_system 
            !! Swiftest nbody system object
        ! Internals
        integer(I4B)    :: i 
        real(DP), dimension(NDIM)  :: g_sph 
            !! Gravitational terms from Spherical Harmonics

        associate(cb => nbody_system%cb)
            cb%aobl(:) = 0.0_DP
            select type(body)
            class is (swiftest_pl)
                do i = 1, body%nbody
                    if (body%lmask(i)) then
                        call shgrav_g_acc_one(cb%Gmass, cb%radius, cb%rotphase, body%rh(:,i), cb%c_lm, body%aobl(:,i), &
                            GMpl=body%Gmass(i), aoblcb=cb%aobl)
                    end if
                end do
            class is (swiftest_tp)
                do i = 1, body%nbody
                    if (body%lmask(i)) then
                        call shgrav_g_acc_one(cb%Gmass, cb%radius, cb%rotphase, body%rh(:,i), cb%c_lm, body%aobl(:,i))
                    end if
                end do
            end select
        end associate
        return 
    end subroutine shgrav_acc
    
end submodule s_shgrav_accel