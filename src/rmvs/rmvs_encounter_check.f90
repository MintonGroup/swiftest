!! Copyright 2022 - David Minton, Carlisle Wishard, Jennifer Pouplin, Jake Elliott, & Dana Singh
!! This file is part of Swiftest.
!! Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
!! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!! Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
!! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!! You should have received a copy of the GNU General Public License along with Swiftest. 
!! If not, see: https://www.gnu.org/licenses. 

submodule (rmvs) s_rmvs_encounter_check
   use swiftest
contains

   module function rmvs_encounter_check_tp(self, param, nbody_system, dt) result(lencounter)
      !! author: David A. Minton
      !!
      !! Determine whether a test particle and planet are having or will have an encounter within the next time step
      !!
      !! Adapted from David E. Kaufmann's Swifter routine: rmvs_chk.f90
      !! Adapted from Hal Levison's Swift routine rmvs3_chk.f
      implicit none
      ! Arguments
      class(rmvs_tp),             intent(inout) :: self   !! RMVS test particle object  
      class(swiftest_parameters), intent(inout) :: param  !! Current swiftest run configuration parameters
      class(rmvs_nbody_system),   intent(inout) :: nbody_system !! RMVS nbody system object
      real(DP),                   intent(in)    :: dt     !! step size
      ! Result
      logical                                 :: lencounter  !! Returns true if there is at least one close encounter
      ! Internals
      integer(I4B)                            :: i
      integer(I8B)                            :: nenc
      logical, dimension(:),      allocatable :: lvdotr
      integer(I4B), dimension(:), allocatable :: index1, index2

      ! The minimization and linear solvers can sometimes lead to floating point exceptions. Rather than halting the code entirely if this occurs, we
      ! can simply fail the attempt and try again. So we need to turn off any floating point exception halting modes temporarily 

      lencounter = .false.
      if (self%nbody == 0) return

      select type(pl => nbody_system%pl)
      class is (rmvs_pl)
         associate(tp => self, ntp => self%nbody, npl => pl%nbody, cb => nbody_system%cb)
            tp%plencP(1:ntp) = 0
            call encounter_check_all_pltp(param, npl, ntp, pl%rbeg, pl%vbeg, tp%rh, tp%vh, pl%renc, dt, & 
                                          nenc, index1, index2, lvdotr)

            lencounter = (nenc > 0_I8B)
            if (lencounter) then
               tp%plencP(index2(1_I8B:nenc)) = index1(1_I8B:nenc)
               do i = 1, npl
                  pl%nenc(i) = count(tp%plencP(1:ntp) == i)
               end do
            end if
         end associate
      end select

      return
   end function rmvs_encounter_check_tp


end submodule s_rmvs_encounter_check
