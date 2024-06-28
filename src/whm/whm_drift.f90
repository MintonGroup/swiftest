! Copyight 2022 - David Minton, Carlisle Wishard, Jennifer Pouplin, Jake Elliott, & Dana Singh
! This file is part of Swiftest.
! Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
! Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
! You should have received a copy of the GNU General Public License along with Swiftest. 
! If not, see: https://www.gnu.org/licenses. 

submodule(whm) whm_drift
   use swiftest
contains

   module subroutine whm_drift_pl(self, nbody_system, param, dt)
      !! author: David A. Minton
      !!
      !! Loop through planets and call Danby drift routine
      !!
      !! Adapted from Hal Levison's Swift routine drift.f
      !! Adapted from David E. Kaufmann's Swifter routine whm_drift.f90
      implicit none
      ! Arguments
      class(whm_pl),                 intent(inout) :: self   !! WHM massive body particle data structure
      class(swiftest_nbody_system),  intent(inout) :: nbody_system !! WHM nbody system object
      class(swiftest_parameters),    intent(in)    :: param  !! Current run configuration parameters 
      real(DP),                      intent(in)    :: dt     !! Stepsize
      ! Internals
      integer(I4B)                              :: i
      integer(I4B), dimension(:), allocatable   :: iflag
      character(len=STRMAX) :: message

      if (self%nbody == 0) return

      associate(pl => self, npl => self%nbody)
         allocate(iflag(npl))
         iflag(:) = 0
         call swiftest_drift_all(pl%muj, pl%xj, pl%vj, npl, param, dt, pl%lmask, iflag)
         if (any(iflag(1:npl) /= 0)) then
            where(iflag(1:npl) /= 0) 
               pl%status(1:npl) = DISCARDED_DRIFTERR
               pl%lmask(1:npl) = .false.
            end where
            do i = 1, npl
               if (iflag(i) /= 0) then 
                  write(message, *) " Planet ", pl%id(i), " is lost!!!!!!!!!!!!",new_line('A'), &
                                      pl%muj(i), dt, new_line('A'), &
                                      pl%xj(:,i),new_line('A'), &
                                      pl%vj(:,i),new_line('A'), &
                                    " STOPPING "
                  call swiftest_io_log_one_message(COLLISION_LOG_OUT, message)
               end if
            end do
            call base_util_exit(FAILURE,param%display_unit)
         end if
      end associate

      return
   end subroutine whm_drift_pl

end submodule whm_drift
