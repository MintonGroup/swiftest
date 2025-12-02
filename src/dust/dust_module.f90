!! Copyright 2025 - The Minton Group at Purdue University
!! This file is part of Swiftest.
!! Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
!! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!! Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
!! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!! You should have received a copy of the GNU General Public License along with Swiftest. 
!! If not, see: https://www.gnu.org/licenses. 


!! Swiftest submodule to model dust particles affected by radiation forces

module dust
    !! author: Kaustub P Anand and David A. Minton
    !!
    !! This module defines functions used for the computation of radiation effects on massive bodies.
    !! Equations taken from Burns, Lamy & Soter (1979) Icarus 40, 1-48.

    use swiftest
    use whm
    use helio
    use symba
    implicit none
    public

    type, extends()

    type, extends(whm_tp) :: whm_dust
       real(DP),                dimension(:),   allocatable :: mass    
          !! Body mass (units MU)
       real(DP),                dimension(:),   allocatable :: Gmass   
          !! Mass gravitational term G * mass (units GU * MU)
       real(DP),                dimension(:),   allocatable :: density
          !! Body mass density - calculated internally (units MU / DU**3)
    contains
        procedure :: accel_pr => whm_dust_accel_pr
    end type whm_dust


    type, extends(helio_tp) :: helio_dust
       real(DP),                dimension(:),   allocatable :: mass    
          !! Body mass (units MU)
       real(DP),                dimension(:),   allocatable :: Gmass   
          !! Mass gravitational term G * mass (units GU * MU)
       real(DP),                dimension(:),   allocatable :: density
          !! Body mass density - calculated internally (units MU / DU**3)
    contains
        procedure :: accel_pr => helio_dust_accel_pr
    end type helio_dust
    


    type, extends(symba_tp) :: symba_dust
      real(DP),                dimension(:),   allocatable :: mass    
         !! Body mass (units MU)
      real(DP),                dimension(:),   allocatable :: Gmass   
         !! Mass gravitational term G * mass (units GU * MU)
      real(DP),                dimension(:),   allocatable :: density 
         !! Body mass density - calculated internally (units MU / DU**3)
      real(DP),                dimension(:),   allocatable :: radius  
         !! Body radius (units DU)
    contains
        procedure :: accel_pr => symba_dust_accel_pr
    end type symba_dust
    

    interface 
        module subroutine whm_dust_accel_pr(self, nbody_system, param)
            implicit none
            ! Arguments
        class(whm_dust),         intent(inout) :: self
            !! Swiftest body object
        class(whm_nbody_system), intent(inout) :: nbody_system
            !! Swiftest nbody system object
        class(swiftest_parameters),   intent(in)    :: param
            !! Current run configuration parameters
        end subroutine whm_dust_accel_pr

        module subroutine helio_dust_accel_pr(self, nbody_system, param)
            implicit none
            ! Arguments
        class(helio_dust),         intent(inout) :: self
            !! Swiftest body object
        class(helio_nbody_system), intent(inout) :: nbody_system
            !! Swiftest nbody system object
        class(swiftest_parameters),   intent(in)    :: param
            !! Current run configuration parameters
        end subroutine helio_dust_accel_pr


        module subroutine symba_dust_accel_pr(self, nbody_system, param)
            implicit none
            ! Arguments
        class(symba_dust),         intent(inout) :: self
            !! Swiftest body object
        class(symba_nbody_system), intent(inout) :: nbody_system
            !! Swiftest nbody system object
        class(swiftest_parameters),   intent(in)    :: param
            !! Current run configuration parameters
        end subroutine symba_dust_accel_pr

    end interface





end module dust