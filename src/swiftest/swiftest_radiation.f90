!! Copyright 2024 - The Minton Group at Purdue University
!! This file is part of Swiftest.
!! Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
!! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!! Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
!! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!! You should have received a copy of the GNU General Public License along with Swiftest. 
!! If not, see: https://www.gnu.org/licenses. 


!! Swiftest submodule to calculate radiation effects on massive bodies
!! Effects include:
!!  - PR drag + Radiation Pressure
!!  - Yarkovsky effect
!!  - Planetary Yarkovsky + Yarkovsky-Schach effects
!!  - YORP effect

submodule (swiftest) s_swiftest_radiation

contains

    module subroutine swiftest_yarkovsky_getacch_pl_one(lag_angle_constants, mu, mass, radius, r_vec, v_vec, rot, a, emissivity, gamma, albedo, rot_k, L_SUN_sys, inv_c2, a_yark) ! pure module subroutine? 
        !! author: Kaustub P. Anand and David A. Minton
        !! Calculate the Yarkovsky effect on one body 
        !! Based on Ferich, et al, 2022 (https://doi.org/10.3847/1538-4365/ac8d60) and Veras, et al, 2015 (https://doi.org/10.1093/mnras/stv1047)
        !!
        implicit none
        ! Arguments
        real(DP), intent(in)                        :: lag_angle_constants, L_SUN_sys, inv_c2
            !! constants and parameters needed for Yarkovsky calculations
        real(DP), intent(in)                        :: emissivity, gamma, albedo, rot_k
            !! particle characteristics for Yarkovsky calculations
        real(DP), intent(in)                        :: a, mass, radius, mu
            !! semi-major axis, mass, radius, and mu of the particle
        real(DP), dimension(NDIM), intent(in)       :: r_vec, v_vec
            !! position and velocity vectors of the particle
        real(DP), dimension(NDIM), intent(in)       :: rot
            !! rotation vector of the particle
        real(DP), dimension(NDIM), intent(out)      :: a_yark 
            !! Yarkovsky acceleration vector

        ! Internals
        integer(I4B)                    :: j, k
            !! looping index
        real(DP)                        :: phi, zeta
            !! thermal lag angles in the rotational plane and orbital plane respectively
        real(DP)                        :: rmag, vmag, h_mag, s_mag, a_yark_mag
            !! magnitude values for respective vectors
        real(DP)                        :: T_orbit, T_rot
            !! orbital and rotation periods
        real(DP), dimension(NDIM)       :: h
            !! Specific angular momentum vector
        real(DP), dimension(NDIM)       :: i_rad
            !! radiation direction vector
        real(DP), dimension(NDIM, NDIM) :: UM, R_s, R1_s, R2_s, R_h, R1_h, R2_h, Y_dir
            !! rotation matrices

        a_yark(:) = 0.0_DP
        UM(:, :) = 0.0_DP
        UM(1, 1) = 1.0_DP
        UM(2, 2) = 1.0_DP
        UM(3, 3) = 1.0_DP

        rmag = .mag. r_vec(:)
        vmag = .mag. v_vec(:) 
        
        h(:) = r_vec(:) .cross. v_vec(:)
        h_mag = .mag. h(:)
        s_mag = .mag. rot(:)
        T_rot = 2 * PI / s_mag ! TU
        T_orbit = 2*PI*a**(1.5_DP) / sqrt(mu) ! orbital period
        
        ! calculate thermal lag angles from eqn. 19 and 20 in Veras, et. al. (2022)
        phi = atan2(1.0_DP, 1.0_DP + lag_angle_constants * emissivity**(0.25_DP) * T_rot**(0.5_DP) / gamma * (1 - albedo)**(0.75_DP) / rmag**(1.5_DP))
        zeta = atan2(1.0_DP, 1.0_DP + lag_angle_constants * emissivity**(0.25_DP) * T_orbit**(0.5_DP) / gamma * (1 - albedo)**(0.75_DP) / rmag**(1.5_DP))

        ! rotation matrices using MATMUL; left for potential future restructuring
        ! R2_s(:, :) = matmul(rot(:), rot(:)) / s_mag**2! rot(:) .cross. rot(:) / s_mag**2
        ! R2_h(:, :) = matmul(h(:), h(:)) / h_mag**2 !h(:) .cross. h(:) / h_mag**2

        ! Calculate R_1 matrices from eqn. 15 and 17 in Veras, et. al. (2022)
        R1_s(1, :) = [0.0_DP, -rot(3), rot(2)] / s_mag
        R1_s(2, :) = [rot(3), 0.0_DP, -rot(1)] / s_mag
        R1_s(3, :) = [-rot(2), rot(1), 0.0_DP] / s_mag

        R1_h(1, :) = [0.0_DP, -h(3), h(2)] / h_mag
        R1_h(2, :) = [h(3), 0.0_DP, -h(1)] / h_mag
        R1_h(3, :) = [-h(2), h(1), 0.0_DP] / h_mag

        ! Calculate R_2 matrices from eqn. 16 and 18 in Veras, et. al. (2022)
        R2_s(1, :) = [rot(1)**2, rot(1)*rot(2), rot(1)*rot(3)] / s_mag**2
        R2_s(2, :) = [rot(1)*rot(2), rot(2)**2, rot(2)*rot(3)] / s_mag**2
        R2_s(3, :) = [rot(1)*rot(3), rot(2)*rot(3), rot(3)**2] / s_mag**2

        R2_h(1, :) = [h(1)**2, h(1)*h(2), h(1)*h(3)] / h_mag**2
        R2_h(2, :) = [h(1)*h(2), h(2)**2, h(2)*h(3)] / h_mag**2
        R2_h(3, :) = [h(1)*h(3), h(2)*h(3), h(3)**2] / h_mag**2

        ! check for and remove very small numbers to 0 to avoid floating underflow errors in rotation matrix calculations
        do j=1, NDIM
            do k=1, NDIM
                if (abs(R1_s(j, k)) <= EPSILON(0.0_DP)) then
                    R1_s(j, k) = 0.0_DP
                end if
                if (abs(R1_h(j, k)) <= EPSILON(0.0_DP)) then
                    R1_h(j, k) = 0.0_DP
                end if
                if (abs(R2_s(j, k)) <= EPSILON(0.0_DP)) then
                    R2_s(j, k) = 0.0_DP
                end if
                if (abs(R2_h(j, k)) <= EPSILON(0.0_DP)) then
                    R2_h(j, k) = 0.0_DP
                end if
            end do
        end do

        ! Combined rotation matrices
        R_s(:, :) = cos(phi) * UM(:, :) + sin(phi) * R1_s(:, :) + (1.0_DP - cos(phi)) * R2_s(:, :)
        R_h(:, :) = cos(zeta) * UM(:, :) - sin(zeta) * R1_h(:, :) + (1.0_DP - cos(zeta)) * R2_h(:, :)

        !! We will assume that v << c, so radiation direction vector is r_hat. If not:
        ! if vmag**2 * param%inv_c2 > 1e-3 then
        !     i_rad(:) = (1 - dot_product(pl%vh(:, i), pl%rh(:, i)) * sqrt(param%inv_c2) / rmag) * pl%rh(:, i) / rmag - pl%vh(:, i) * sqrt(param%inv_c2) ! radiation direction vector
        ! end if

        i_rad(:) = .unit. r_vec(:)! radiation direction vector

        ! yark acceleration magnitude from eqn. 1 in Ferich, et al (2022) / eqn. 26 in Veras, et al (2015)
        a_yark_mag = rot_k * radius**2 * (1.0_DP - albedo) * L_SUN_sys * sqrt(inv_c2) / (4.0_DP * mass * rmag**2)

        ! calculate yarkovsky direction matrix
        Y_dir = matmul(R_s(:, :), R_h(:, :))

        ! Multiply yarkovsky direction matrix with radiation direction vector
        do j = 1, NDIM
            a_yark(j) = Y_dir(j, 1) * i_rad(1) + Y_dir(j, 2) * i_rad(2) + Y_dir(j, 3) * i_rad(3)
        end do

        a_yark(:) = a_yark_mag * a_yark(:) 

        return

    end subroutine swiftest_yarkovsky_getacch_pl_one

    module subroutine swiftest_yarkovsky_getacch_pl_all(nbody, lmask, mu, mass, radius, r_vec, v_vec, acc, rot, a, emissivity, gamma, albedo, rot_k, L_SUN_sys, inv_c2, sigma_sys)
        !! author: Kaustub P. Anand and David A. Minton
        !! Loop over all bodies to calculate the Yarkovsky effect. 
        !! Based on Ferich, et al, 2022 (https://doi.org/10.3847/1538-4365/ac8d60) and Veras, et al, 2015 (https://doi.org/10.1093/mnras/stv1047)
        !!
        implicit none
        ! Arguments
        integer(I4B), intent(in)                        :: nbody
            !! number of bodies in the system)
        logical, dimension(:), intent(in)          :: lmask
            !! logical mask for active bodies in the system
        real(DP), intent(in)                            :: L_SUN_sys, inv_c2, sigma_sys
            !! constants and parameters needed for Yarkovsky calculations
        real(DP), dimension(:), intent(in)              :: emissivity, gamma, albedo, rot_k
            !! particle characteristics for Yarkovsky calculations
        real(DP), dimension(:), intent(in)              :: a, mass, radius, mu
            !! semi-major axis, mass, radius, and mu of the particle
        real(DP), dimension(:, :), intent(in)           :: r_vec, v_vec
            !! position and velocity vectors of the particle
        real(DP), dimension(:, :), intent(in)           :: rot
            !! rotation vector of the particle
        real(DP), dimension(:, :), intent(inout)        :: acc
            !! Acceleration vector for all bodies

        ! Internals
        integer(I4B)                     :: i
            !! looping index
        real(DP)                         :: lag_angle_constants
            !! constant terms in lag angle calculations
        real(DP), dimension(NDIM)        :: a_yark
            !! Yarkovsky acceleration vector

        ! calculate constants
        lag_angle_constants = 0.5_DP * (sigma_sys / PI**5)**(0.25_DP) * (L_SUN_sys)**(0.75_DP)

        do i=1, nbody
            if (lmask(i) .and. radius(i) * DU2M <= 3e4_DP) then !! check if body radius is <= 30 km for computational efficiency. Yarkovsky effect is negligible for larger bodies (Bottke, et al, 2006; doi:10.1146/annurev.earth.34.031405.125154)
                call swiftest_yarkovsky_getacch_pl_one(lag_angle_constants, mu(i), mass(i), radius(i), r_vec(:, i), v_vec(:, i), rot(:, i) * DEG2RAD, a(i), emissivity(i), gamma(i), albedo(i), rot_k(i), L_SUN_sys, inv_c2, a_yark)
                acc(:, i) = acc(:, i) + a_yark(:)
            end if 
        end do
        
        return

    end subroutine swiftest_yarkovsky_getacch_pl_all

    ! module subroutine swiftest_yarkovsky_schach_getacch_pl(self, nbody_system, param)
    !     !! author: Kaustub P. Anand and David A. Minton
    !     !!
    !     !! Calculate the Yarkovsky-Schach + Planetary Yarkovsky effect on massive bodies. 
    !     !! Based on << >>
    !     implicit none
    !     ! Arguments
    !     class(swiftest_pl),         intent(inout) :: self
    !         !! Swiftest body object
    !     class(swiftest_nbody_system), intent(inout) :: nbody_system
    !         !! Swiftest nbody system object
    !     class(swiftest_parameters),   intent(in)    :: param
    !         !! Current run configuration parameters
    !     ! Internals


    !     associate(pl => self)

    !     end associate

    !     return
    ! end subroutine swiftest_yarkovsky_schach_getacch_pl

    module subroutine swiftest_radiation_getacch_pl(self, nbody_system, param)
        !! author: Kaustub P. Anand and David A. Minton
        !!
        !! Calculate radiation effects (PR-drag + radiation pressure) on massive bodies.
        implicit none
        ! Arguments
        class(swiftest_pl),         intent(inout) :: self
            !! Swiftest body object
        class(swiftest_nbody_system), intent(inout) :: nbody_system
            !! Swiftest nbody system object
        class(swiftest_parameters),   intent(in)    :: param
            !! Current run configuration parameters
        ! ! Internals
        ! integer(I4B)    :: i
        !     !! looping index
        ! real(DP)        :: Q_pr
        !     !! Radiation pressure efficiency factor
        !     !! assumed to be 1.0; Krivov, et al, 1996. http://dx.doi.org/10.1007/BF00692293 gives a good reasoning
        ! real(DP)        :: rmag, vmag
        !     !! magnitude of position and velocity vectors
        ! real(DP), dimension(NDIM) :: S_vec
        !     !! Solar radiation flux vector
        ! real(DP)        :: fac1
        !     !! combined SA/mc factor for acceleration calculation

        ! Q_pr = 1.0_DP ! placeholder in case this needs to changed in the future

        ! associate(pl => self)
        !     do i=1, pl%nbody
        !         if (pl%lmask(i)) then
        !             rmag = sqrt(dot_product(pl%rh(:, i), pl%rh(:, i)))
        !             vmag = sqrt(dot_product(pl%vh(:, i), pl%vh(:, i)))
        !             S_vec(:) = - pl%rh(:, i) / rmag * param%L_SUN_sys / (4.0_DP * PI * rmag**2) ! S_hat = - pl%rh(:, i)
                    
        !             fac1 = param%L_SUN_sys * (param%TU2S)**3 / (param%MU2KG * param%DU2M**2) * sqrt(param%inv_c2) * pl%radius(i)**2 / (4.0_DP * pl%mass(i) * rmag**2) ! SA/mc = L_sun * radius^2 / (4 * c * distance^2 * pl_mass)

        !             pl%ah(:, i) = pl%ah(:, i) + fac1 * Q_pr * ((vmag * param%inv_c - 1.0_DP) * pl%rh(:, i) / rmag - pl%vh(:, i) * param%inv_c) ! eqn. 5 in Burns, et al, 1979. https://doi.org/10.1016/0019-1035(79)90050-2; ICARUS 40, 1 - 48 (1979)

        !         end if
        !     end do
        ! end associate
            

        return
    end subroutine swiftest_radiation_getacch_pl

    ! module subroutine swiftest_yorp_getacch_pl(self, nbody_system, param)
    !     !! author: Kaustub P. Anand and David A. Minton
    !     !!
    !     !! Calculate the Yorp effect on massive bodies. 
    !     !! Based on << >>
    !     implicit none
    !     ! Arguments
    !     class(swiftest_pl),         intent(inout) :: self
    !         !! Swiftest body object
    !     class(swiftest_nbody_system), intent(inout) :: nbody_system
    !         !! Swiftest nbody system object
    !     class(swiftest_parameters),   intent(in)    :: param
    !         !! Current run configuration parameters
    !     ! Internals


    !     associate(pl => self)

    !     end associate

    !     return
    ! end subroutine swiftest_yorp_getacch_pl

end submodule s_swiftest_radiation