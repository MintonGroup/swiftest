! Copyight 2022 - David Minton, Carlisle Wishard, Jennifer Pouplin, Jake Elliott, & Dana Singh
! This file is part of Swiftest.
! Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
! Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
! You should have received a copy of the GNU General Public License along with Swiftest. 
! If not, see: https://www.gnu.org/licenses. 

submodule(rmvs) s_rmvs_discard
   use swiftest
contains  

   module subroutine rmvs_discard_tp(self, nbody_system, param)
      !! author: David A. Minton
      !!
      !! Check to see if test particles should be discarded based on pericenter passage distances with respect to planets encountered
      !!
      !! Adapted from Hal Levison's Swift routine discard_pl.f
      !! Adapted from Hal Levison's Swift routine rmvs_discard_pl.f90
      implicit none
      ! Arguments
      class(rmvs_tp),               intent(inout) :: self   !! RMVS test particle object
      class(swiftest_nbody_system), intent(inout) :: nbody_system !! Swiftest nbody system object
      class(swiftest_parameters),   intent(inout) :: param  !! Current run configuration parameters 
      ! Internals
      integer(I4B)                                 :: i
      character(len=STRMAX) :: timestr, idstri, idstrj, message

      if (self%nbody == 0) return

      associate(tp => self, ntp => self%nbody, pl => nbody_system%pl, t => nbody_system%t)
         do i = 1, ntp
            associate(iplperP => tp%plperP(i))
               if ((tp%status(i) == ACTIVE) .and. (tp%lperi(i))) then 
                  if ((tp%peri(i) < pl%radius(iplperP))) then
                     tp%status(i) = DISCARDED_PLQ
                     write(idstri, *) tp%id(i)
                     write(idstrj, *) pl%id(iplperP)
                     write(timestr, *) t
                     write(message, *) "Particle "  // trim(adjustl(tp%info(i)%name)) // " (" // trim(adjustl(idstri)) &
                              // ") q with respect to massive body " // trim(adjustl(pl%info(iplperP)%name))   &
                              // " (" // trim(adjustl(idstrj)) // ") is too small at t = " // trim(adjustl(timestr))
                     call swiftest_io_log_one_message(COLLISION_LOG_OUT,message)
                     tp%ldiscard(i) = .true.
                     tp%lmask(i) = .false.
                     pl%ldiscard(iplperP) = .true.
                     call tp%info(i)%set_value(status="DISCARDED_PLQ", discard_time=t, discard_rh=tp%rh(:,i), &
                                               discard_vh=tp%vh(:,i), discard_body_id=pl%id(iplperP))
                  end if
               end if
            end associate
         end do
         ! Call the base method that this overrides
         call swiftest_discard_tp(tp, nbody_system, param)
      end associate


      return

   end subroutine rmvs_discard_tp

end submodule s_rmvs_discard