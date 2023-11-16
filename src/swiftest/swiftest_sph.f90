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

    module subroutine swiftest_sph_g_acc_one(GMcb, r_0, phi_cb, rh, c_lm, g_sph, GMpl, aoblcb)
        !! author: Kaustub P. Anand
        !!
        !! Calculate the acceleration terms for one pair of bodies given c_lm, theta, phi, r
        !!

        implicit none
        ! Arguments
        real(DP), intent(in)        :: GMcb                        !! GMass of the central body
        real(DP), intent(in)        :: r_0                         !! radius of the central body
        real(DP), intent(in)        :: phi_cb                      !! rotation phase of the central body
        real(DP), intent(in), dimension(:)          :: rh          !! distance vector of body
        real(DP), intent(in), dimension(:, :, :)    :: c_lm        !! Spherical Harmonic coefficients
        real(DP), intent(out), dimension(NDIM)         :: g_sph    !! acceleration vector
        real(DP), intent(in),  optional :: GMpl                    !! Mass of input body if it is not a test particle
        real(DP), dimension(:),   intent(inout), optional :: aoblcb!! Barycentric acceleration of central body (only for massive input bodies)
     
        ! Internals
        integer        :: l, m              !! SPH coefficients
        integer        :: l_max             !! max Spherical Harmonic l order value
        integer(I4B)   :: N, lmindex        !! Length of Legendre polynomials and index at a given l, m
        real(DP)       :: r_mag             !! magnitude of rh
        real(DP)       :: phi, phi_bar      !! Azimuthal/Phase angle (radians) wrt coordinate axes, and central body rotation phase
        real(DP)       :: theta             !! Inclination/Zenith angle (radians)
        real(DP)       :: plm, dplm         !! Associated Legendre polynomials at a given l, m
        real(DP)       :: ccss, cssc        !! See definition in source code
        real(DP)       :: cos_theta, sin_theta !! cos(theta) and sin(theta)
        real(DP), dimension(:), allocatable  :: p, p_deriv   !! Associated Lengendre Polynomials at a given cos(theta)
        real(DP)     :: r2, irh, rinv2, t0, t1, t2, t3, fac0, fac1, fac2, fac3, fac4, j2rp2, j4rp4, r_fac, cos_tmp, sin_tmp, sin2_tmp, plm1, sin_phi, cos_phi

        g_sph(:) = 0.0_DP
        theta = atan2(sqrt(rh(1)**2 + rh(2)**2), rh(3))
        phi = atan2(rh(2), rh(1)) 
        phi_bar = MOD(phi - phi_cb, 2 * PI)
        r_mag = sqrt(dot_product(rh(:), rh(:)))
        l_max = size(c_lm, 2) - 1
        N = (l_max + 1) * (l_max + 2) / 2
        allocate(p(N),p_deriv(N))

        ! to compare w s_obl.f90
        j2rp2 = -c_lm(1, 3, 1) * r_0**2
        j4rp4 = -c_lm(1, 4, 1) * r_0**4
        r2 = dot_product(rh(:), rh(:))
        irh = 1.0_DP / sqrt(r2)
        rinv2 = irh**2
        t0 = -GMcb * rinv2 * rinv2 * irh
        t1 = 1.5_DP * j2rp2
        t2 = rh(3) * rh(3) * rinv2
        t3 = 1.875_DP * j4rp4 * rinv2
        fac1 = t0 * (t1 - t3 - (5 * t1 - (14.0_DP - 21.0_DP * t2) * t3) * t2)
        fac2 = 2 * t0 * (t1 - (2.0_DP - (14.0_DP * t2 / 3.0_DP)) * t3)
        fac0 = 4 * PI

        cos_theta = cos(theta)
        sin_theta = sin(theta)



        if(abs(cos_theta) < epsilon(0.0_DP)) then
            cos_theta = 0.0_DP
        end if
        if(abs(sin_theta) < epsilon(0.0_DP)) then
            sin_theta = 0.0_DP
        end if

        ! call PlmBar_d1(p, p_deriv, l_max, cos_theta)      ! Associated Legendre Polynomials and the 1st Derivative
        call PlmBar(p, l_max, cos_theta)


        do l = 1, l_max ! skipping the l = 0 term; It is the spherical body term
            do m = 0, l

                ! Associated Legendre Polynomials 
                lmindex = PlmIndex(l, m)  
                plm = p(lmindex)                ! p_l,m
                ! dplm = p_deriv(lmindex)         ! d(p_l,m)

                ! C_lm and S_lm with Cos and Sin of m * phi
                ccss = c_lm(m+1, l+1, 1) * cos(m * phi_bar) & 
                        + c_lm(m+1, l+1, 2) * sin(m * phi_bar)      ! C_lm * cos(m * phi) + S_lm * sin(m * phi)
                cssc = -1 * c_lm(m+1, l+1, 1) * sin(m * phi_bar) & 
                        + c_lm(m+1, l+1, 2) * cos(m * phi_bar)      ! - C_lm * sin(m * phi) + S_lm * cos(m * phi) 
                                                                ! cssc * m = first derivative of ccss with respect to phi

                ! m > 0
                ! g_sph(1) = g_sph(1) - GMcb * r_0**l / r_mag**(l + 2) * (cssc * m * plm * sin(phi) / sin_theta &
                !                                                         + ccss * sin_theta * cos(phi) &    
                !                                                         * (dplm * cos_theta + plm * (l + 1))) ! g_x
                ! g_sph(2) = g_sph(2) - GMcb * r_0**l / r_mag**(l + 2) * (-1 * cssc * m * plm * cos(phi) / sin_theta &
                !                                                         + ccss * sin_theta * sin(phi) &
                !                                                         * (dplm * cos_theta + plm * (l + 1))) ! g_y
                ! g_sph(3) = g_sph(3) - GMcb * r_0**l / r_mag**(l + 2) * ccss * (-1 * dplm * sin_theta**2  &
                !                                                         + plm * (l + 1) * cos_theta)          ! g_z

                cos_tmp = cos_theta
                sin_tmp = sin_theta
                sin2_tmp = sin(2 * theta)
                sin_phi = sin(phi)
                cos_phi = cos(phi)
                ! g_sph(:) = 0.0_DP 

                ! Alternative form for g_sph

                if ((m+1) .le. l) then
                    lmindex = PlmIndex(l, m+1) 
                    plm1 = p(lmindex) 
                    if(m .eq. 0) then
                        plm1 = plm1 * sqrt(((l + m + 1) * (l - m)) / 2.0) ! renormalize plm1 to the norm of plm
                    else 
                        plm1 = plm1 * sqrt((l + m + 1) * (l - m) * 1.0)       ! renormalize plm1 to the norm of plm
                    end if
                else
                    plm1 = 0.0_DP  
                end if 
                                                                      
                ! !! Alternative form of dplm

                ! g_sph(1) = g_sph(1) - GMcb * r_0**l / r_mag**(l + 2) * (cssc * m / sin_theta * plm * sin(phi) &
                !                                                        + ccss * cos(phi) * (plm * ((l + m + 1) * sin_theta - m / sin_theta) & 
                !                                                                              + plm1 * cos_theta))
                ! g_sph(2) = g_sph(2) - GMcb * r_0**l / r_mag**(l + 2) * (-cssc * m / sin_theta * plm * cos(phi) &
                !                                                         + ccss * sin(phi) * (plm * ((l + m + 1) * sin_theta - m / sin_theta) & 
                !                                                                              + plm1 * cos_theta)) 
                ! g_sph(3) = g_sph(3) - GMcb * r_0**l / r_mag**(l + 2) * ccss * (plm * (l + m +1) * cos_theta - plm1 * sin_theta)          

                ! Condensed form

                ! fac0 = -(m * cos_tmp * plm / sin_tmp - plm1) / sin_tmp ! dplm
                ! fac3 = plm * (l + m + 1) * cos_theta
                ! fac4 = plm1 * sin_theta
                if(sin_theta .eq. 0) then
                    fac1 = 0.0_DP
                else
                    fac1 = m * plm / sin_theta
                end if

                fac2 = plm * (l + m + 1) * sin_theta + plm1 * cos_theta
                fac3 = fac2 - fac1
                r_fac = -GMcb * r_0**l / r_mag**(l + 2)

                ! g_sph(:) = 0.0_DP
                g_sph(1) = g_sph(1) + r_fac * (cssc * fac1 * sin(phi) + ccss * (fac2 - fac1) * cos(phi))
                g_sph(2) = g_sph(2) + r_fac * (-cssc * fac1 * cos(phi) + ccss * (fac2 - fac1) * sin(phi))
                g_sph(3) = g_sph(3) + r_fac * ccss * (plm * (l + m + 1) * cos_theta - plm1 * sin_theta)

                fac0 = (.mag. g_sph(:))
                ! fac3 = 3 * c_lm(m+1, l+1, 1) / 2 * r_0**2 / r_mag**4 * (3*(cos_theta)**2 - 1) * GMcb ! g_sph for J2
                            
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
        real(DP)        :: theta, phi              !! zenith angle, and azimuthal angle
        real(DP), dimension(NDIM)  :: g_sph        !! Gravitational terms from Spherical Harmonics

        associate(pl => self, npl => self%nbody, cb => nbody_system%cb, rh => self%rh)
            cb%aobl(:) = 0.0_DP

            do i = 1, npl
                if (pl%lmask(i)) then
                    ! theta = atan2(sqrt(rh(1,i)**2 + rh(2,i)**2), rh(3,i))
                    ! phi = atan2(rh(2,i), rh(1,i)) ! - cb%rotphase
    
                    call swiftest_sph_g_acc_one(cb%Gmass, cb%radius, cb%rotphase, rh(:,i), cb%c_lm, g_sph, pl%Gmass(i), cb%aobl)
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
        real(DP)        :: theta, phi              !! zenith angle, and azimuthal angle
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

                    call swiftest_sph_g_acc_one(cb%Gmass, cb%radius, cb%rotphase, rh(:,i), cb%c_lm, g_sph)
                    tp%ah(:, i) = tp%ah(:, i) + g_sph(:) - aoblcb(:)
                    tp%aobl(:, i) = g_sph(:)
                end if
            end do
        end associate
        return
        end subroutine swiftest_sph_g_acc_tp_all
    
end submodule s_swiftest_sph