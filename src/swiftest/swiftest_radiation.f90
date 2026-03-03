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

    module subroutine swiftest_yarkovsky_getacch_pl(self, nbody_system, param)
        !! author: Kaustub P. Anand and David A. Minton
        !!
        !! Calculate the Yarkovsky effect on massive bodies. 
        !! Based on Ferich, et al, 2022 (https://iopscience.iop.org/article/10.3847/1538-4365/ac8d60) and Veras, et al, 2015 (https://academic.oup.com/mnras/article/451/3/2814/1180328)
        implicit none
        ! Arguments
        class(swiftest_pl),           intent(inout) :: self
            !! Swiftest body object
        class(swiftest_nbody_system), intent(inout) :: nbody_system
            !! Swiftest nbody system object
        class(swiftest_parameters),   intent(in)    :: param
            !! Current run configuration parameters
        ! Internals
        integer(I4B)                    :: i, j, k
            !! looping index
        real(DP)                        :: phi, zeta
            !! thermal lag angles in the rotational plane and orbital plane respectively
        real(DP)                        :: rmag, vmag, h_mag, s_mag, a_yark_mag
            !! magnitude values for respective vectors
        real(DP)                        :: lag_angle_constants
            !! constant terms in lag angle calculations
        real(DP)                        :: T_orbit, T_rot
            !! orbital and rotation periods
        real(DP), dimension(NDIM)       :: h
            !! Specific angular momentum vector
        real(DP), dimension(NDIM)       :: i_rad
            !! radiation direction vector
        real(DP), dimension(NDIM, 1)       :: a_yark
            !! Yarkovsky acceleration vector
        real(DP), dimension(NDIM, NDIM) :: UM, R_s, R1_s, R2_s, R_h, R1_h, R2_h
            !! rotation matrices

        ! calculate constants
        lag_angle_constants = 0.5_DP * (param%sigma_sys / PI**5)**(0.25_DP) * (param%L_SUN_sys)**(0.75_DP)
        UM(:, :) = 0.0_DP
        UM(1, 1) = 1.0_DP
        UM(2, 2) = 1.0_DP
        UM(3, 3) = 1.0_DP

        associate(pl => self)
            do i=1, pl%nbody
                if (pl%lmask(i)) then
                    rmag = .mag. pl%rh(:, i)
                    vmag = .mag. pl%vh(:, i) 
                    
                    !! vb vs vh AND/OR rh vs rb; See 1255-1257 in swiftest_util.f90
                    !! should h be made into a variable to store per body?
                    h(:) = pl%rh(:, i) .cross. pl%vh(:, i) 
                    h_mag = .mag. h(:)
                    s_mag = .mag. pl%rot(:, i) ! DEG/TU
                    T_rot = 360.0_DP / s_mag ! TU
                    T_orbit = 2*PI*pl%a(i)**(1.5_DP) / sqrt(pl%mu(i)) ! orbital period
                    
                    ! calculate thermal lag angles from eqn. 19 and 20 in Veras, et. al. (2022)
                    phi = atan2(1.0_DP, 1.0_DP + lag_angle_constants * pl%emissivity(i)**(0.25_DP) * T_rot**(0.5_DP) / pl%gamma(i) * (1 - pl%albedo(i))**(0.75_DP) / rmag**(1.5_DP))
                    zeta = atan2(1.0_DP, 1.0_DP + lag_angle_constants * pl%emissivity(i)**(0.25_DP) * T_orbit**(0.5_DP) / pl%gamma(i) * (1 - pl%albedo(i))**(0.75_DP) / rmag**(1.5_DP))

                    ! rotation matrices
                    ! R2_s(:, :) = matmul(pl%rot(:, i), pl%rot(:, i)) / s_mag**2! pl%rot(:, i) .cross. pl%rot(:, i) / s_mag**2
                    ! R2_h(:, :) = matmul(h(:), h(:)) / h_mag**2 !h(:) .cross. h(:) / h_mag**2

                    ! Calculate R_1 matrices from eqn. 15 and 17 in Veras, et. al. (2022)
                    R1_s(1, :) = [0.0_DP, -pl%rot(3, i), pl%rot(2, i)] / s_mag !! CHECK row vs column ordering
                    R1_s(2, :) = [pl%rot(3, i), 0.0_DP, -pl%rot(1, i)] / s_mag
                    R1_s(3, :) = [-pl%rot(2, i), pl%rot(1, i), 0.0_DP] / s_mag

                    R1_h(1, :) = [0.0_DP, -h(3), h(2)] / h_mag
                    R1_h(2, :) = [h(3), 0.0_DP, -h(1)] / h_mag
                    R1_h(3, :) = [-h(2), h(1), 0.0_DP] / h_mag

                    ! Calculate R_2 matrices from eqn. 16 and 18 in Veras, et. al. (2022)
                    R2_s(1, :) = [pl%rot(1, i)**2, pl%rot(1, i)*pl%rot(2, i), pl%rot(1, i)*pl%rot(3, i)] / s_mag**2
                    R2_s(2, :) = [pl%rot(1, i)*pl%rot(2, i), pl%rot(2, i)**2, pl%rot(2, i)*pl%rot(3, i)] / s_mag**2
                    R2_s(3, :) = [pl%rot(1, i)*pl%rot(3, i), pl%rot(2, i)*pl%rot(3, i), pl%rot(3, i)**2] / s_mag**2

                    R2_h(1, :) = [h(1)**2, h(1)*h(2), h(1)*h(3)] / h_mag**2
                    R2_h(2, :) = [h(1)*h(2), h(2)**2, h(2)*h(3)] / h_mag**2
                    R2_h(3, :) = [h(1)*h(3), h(2)*h(3), h(3)**2] / h_mag**2

                    ! check for and remove very small numbers to 0 to avoid floating underflow errors in rotation matrix calculations
                    for j=1, NDIM
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

                    !! We will assume that v << c, so radiation direction vector is r_hat
                    ! if vmag**2 * param%inv_c2 > 1e-3 then
                    !     i_rad(:) = (1 - dot_product(pl%vh(:, i), pl%rh(:, i)) * sqrt(param%inv_c2) / rmag) * pl%rh(:, i) / rmag - pl%vh(:, i) * sqrt(param%inv_c2) ! radiation direction vector
                    ! end if

                    i_rad(:) = .unit. pl%rh(:, i)! radiation direction vector

                    ! yark acceleration magnitude from eqn. 1 in Ferich, et al (2022) / eqn. 26 in Veras, et al (2015)
                    a_yark_mag = pl%rot_k(i) * pl%radius(i)**2 * (1.0_DP - pl%albedo(i)) * param%L_SUN_sys * sqrt(param%inv_c2) / (4.0_DP * pl%mass(i) * rmag**2)

                    ! calculate acceleration
                    ! a_yark(:, 1) = a_yark_mag * matmul(matmul(R_s(:, :), R_h(:, :)), i_rad(:))
                    a_yark(:, 1) = matmul(matmul(R_s(:, :), R_h(:, :)), i_rad(:))
                    ! a_yark(:, 1) = matmul(R_s(:, :), i_rad(:))
                    a_yark(:, 1) = a_yark_mag * a_yark(:, 1) 
                    ! a_yark(i) = a_yark_mag * matmul(R_s(:, :), R_h(:, :))

                    ! add to acceleration
                    pl%ah(:, i) = pl%ah(:, i) + a_yark(:, 1)
                    
                end if
            end do

        end associate

        return
    end subroutine swiftest_yarkovsky_getacch_pl

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