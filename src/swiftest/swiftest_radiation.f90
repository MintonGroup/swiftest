!! Copyright 2024 - The Minton Group at Purdue University
!! This file is part of Swiftest.
!! Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
!! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!! Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
!! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!! You should have received a copy of the GNU General Public License along with Swiftest. 
!! If not, see: https://www.gnu.org/licenses. 


!! Swiftest submodule to calculate radiation effects on massive bodies

submodule (swiftest) s_swiftest_radiation
contains

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
        ! Internals
        integer(I4B)    :: i
            !! looping index
        ! real(DP)        :: L_sun
        !     !! Solar luminosity 
        real(DP)        :: Q_pr
            !! Radiation pressure efficiency factor
            !! assumed to be 1.0; Krivov, et al, 1996. http://dx.doi.org/10.1007/BF00692293 gives a good reasoning
        real(DP)        :: rmag, vmag
            !! magnitude of position and velocity vectors
        real(DP), dimension(NDIM) :: S_vec
            !! Solar radiation flux vector
        real(DP)        :: fac1
            !! combined SA/mc factor for acceleration calculation

        Q_pr = 1.0_DP ! placeholder in case this needs to changed in the future
        ! L_sun = L_SUN * (param%TU2S)**3 / (param%MU2KG * param%DU2M**2) ! 3.828e26 W; Mamajek, et al (2015). IAU 2015 Resolution B3. https://doi.org/10.48550/arXiv.1510.07674

        associate(pl => self)
            do i=1, pl%nbody
                if (pl%lmask(i)) then
                    rmag = sqrt(dot_product(pl%rh(:, i), pl%rh(:, i)))
                    vmag = sqrt(dot_product(pl%vh(:, i), pl%vh(:, i)))
                    S_vec(:) = - pl%rh(:, i) / rmag * L_sun / (4.0_DP * PI * rmag**2) ! S_hat = - pl%rh(:, i)
                    
                    fac1 = L_SUN * (param%TU2S)**3 / (param%MU2KG * param%DU2M**2) * sqrt(param%inv_c2) * pl%radius(i)**2 / (4.0_DP * pl%mass(i) * rmag**2) ! SA/mc = L_sun * radius^2 / (4 * c * distance^2 * pl_mass)

                    pl%ah(:, i) = pl%ah(:, i) + fac1 * Q_pr * ((vmag * param%inv_c - 1.0_DP) * pl%rh(:, i) / rmag - pl%vh(:, i) * param%inv_c) ! eqn. 5 in Burns, et al, 1979. https://doi.org/10.1016/0019-1035(79)90050-2; ICARUS 40, 1 - 48 (1979)

                end if
            end do
        end associate
            

        return
    end subroutine swiftest_radiation_getacch_pl

end submodule s_swiftest_radiation