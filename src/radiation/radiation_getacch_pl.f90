!! Copyright 2024 - The Minton Group at Purdue University
!! This file is part of Swiftest.
!! Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
!! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!! Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
!! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!! You should have received a copy of the GNU General Public License along with Swiftest. 
!! If not, see: https://www.gnu.org/licenses. 


!! Swiftest submodule to calculate radiation effects on massive bodies

submodule (radiation) s_radiation
use swiftest

contains

    module subroutine radiation_getacch_pl(self, nbody_system, param)
        !! author: Kaustub P. Anand and David A. Minton
        !!
        !! Calculate radiation effects on massive bodies.
        implicit none
        ! Arguments
        class(base_object),         intent(inout) :: self
            !! Swiftest body object
        class(base_nbody_system), intent(inout) :: nbody_system
            !! Swiftest nbody system object
        class(base_parameters),   intent(in)    :: param
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

        associate(body => self)
            select type(body)
            class is (swiftest_pl)
            select type(nbody_system)
            class is (swiftest_nbody_system)
            select type(param)
            class is (swiftest_parameters)
                do i=1, body%nbody
                    if (body%lmask(i)) then
                        rmag = sqrt(dot_product(body%rh(:, i), body%rh(:, i)))
                        vmag = sqrt(dot_product(body%vh(:, i), body%vh(:, i)))
                        S_vec(:) = - body%rh(:, i) / rmag * L_sun / (4.0_DP * PI * rmag**2) ! S_hat = - body%rh(:, i)
                        
                        fac1 = L_SUN * (param%TU2S)**3 / (param%MU2KG * param%DU2M**2) * sqrt(param%inv_c2) * body%radius(i)**2 / (4.0_DP * body%mass(i) * rmag**2) ! SA/mc = L_sun * radius^2 / (4 * c * distance^2 * pl_mass)

                        body%ah(:, i) = body%ah(:, i) + fac1 * Q_pr * ((vmag * param%inv_c - 1.0_DP) * body%rh(:, i) / rmag - body%vh(:, i) * param%inv_c) ! eqn. 5 in Burns, et al, 1979. https://doi.org/10.1016/0019-1035(79)90050-2; ICARUS 40, 1 - 48 (1979)

                    end if
                end do
            
            class is (swiftest_tp)
                ! Do nothing for test particles
            end select

        return
    end subroutine radiation_getacch_pl

end submodule s_radiation