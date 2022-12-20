!! Copyright 2022 - David Minton, Carlisle Wishard, Jennifer Pouplin, Jake Elliott, & Dana Singh
!! This file is part of Swiftest.
!! Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
!! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!! Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
!! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!! You should have received a copy of the GNU General Public License along with Swiftest. 
!! If not, see: https://www.gnu.org/licenses. 

submodule(symba) s_symba_setup
   use swiftest
contains

   module subroutine symba_setup_initialize_system(self, param)
      !! author: David A. Minton
      !!
      !! Initialize an SyMBA nbody system from files and sets up the planetocentric structures.
      !! This subroutine will also sort the massive bodies in descending order by mass
      !! 
      implicit none
      ! Arguments
      class(symba_nbody_system),  intent(inout) :: self    !! SyMBA system object
      class(swiftest_parameters), intent(inout) :: param  !! Current run configuration parameters 

      ! Call parent method
      associate(system => self)
         call helio_setup_initialize_system(system, param)
         call system%pltp_encounter%setup(0_I8B)
         call system%plpl_encounter%setup(0_I8B)
         call system%plpl_collision%setup(0_I8B)
      end associate

      return
   end subroutine symba_setup_initialize_system


   module subroutine symba_setup_pl(self, n, param)
      !! author: David A. Minton
      !!
      !! Allocate SyMBA test particle structure
      !!
      !! Equivalent in functionality to David E. Kaufmann's Swifter routine symba_setup.f90
      implicit none
      ! Arguments
      class(symba_pl),            intent(inout) :: self  !! SyMBA massive body object
      integer(I4B),               intent(in)    :: n     !! Number of particles to allocate space for
      class(swiftest_parameters), intent(in)    :: param !! Current run configuration parameter
      ! Internals
      integer(I4B) :: i

      !> Call allocation method for parent class. In this case, helio_pl does not have its own setup method so we use the base method for swiftest_pl
      call setup_pl(self, n, param) 
      if (n == 0) return

      allocate(self%levelg(n))
      allocate(self%levelm(n))
      allocate(self%isperi(n))
      allocate(self%peri(n))
      allocate(self%atp(n))
      allocate(self%kin(n))


      self%levelg(:) = -1
      self%levelm(:) = -1
      self%isperi(:) = 0
      self%peri(:) = 0.0_DP
      self%atp(:) = 0.0_DP
      call self%reset_kinship([(i, i=1, n)])
      return
   end subroutine symba_setup_pl


   module subroutine symba_setup_encounter_list(self, n)
      !! author: David A. Minton
      !!
      !! A constructor that sets the number of encounters and allocates and initializes all arrays  
      !!
      implicit none
      ! Arguments
      class(symba_encounter), intent(inout) :: self !! SyMBA pl-tp encounter structure
      integer(I8B),           intent(in)    :: n    !! Number of encounters to allocate space for

      call encounter_setup_list(self, n)
      if (n <= 0_I8B) return

      if (allocated(self%level)) deallocate(self%level)
      if (allocated(self%tcollision)) deallocate(self%tcollision)
      allocate(self%level(n))
      allocate(self%tcollision(n))

      self%level(:) = -1
      self%tcollision(:) = 0.0_DP

      return
   end subroutine symba_setup_encounter_list


   module subroutine symba_setup_tp(self, n, param)
      !! author: David A. Minton
      !!
      !! Allocate WHM test particle structure
      !!
      !! Equivalent in functionality to David E. Kaufmann's Swifter routine whm_setup.f90
      implicit none
      ! Arguments
      class(symba_tp),            intent(inout) :: self  !! SyMBA test particle object
      integer(I4B),               intent(in)    :: n     !! Number of particles to allocate space for
      class(swiftest_parameters), intent(in)    :: param !! Current run configuration parameter

      !> Call allocation method for parent class. In this case, helio_tp does not have its own setup method so we use the base method for swiftest_tp
      call setup_tp(self, n, param) 
      if (n == 0) return

      allocate(self%levelg(n))
      allocate(self%levelm(n))

      self%levelg(:) = -1
      self%levelm(:) = -1
      
      return
   end subroutine symba_setup_tp

end submodule s_symba_setup
