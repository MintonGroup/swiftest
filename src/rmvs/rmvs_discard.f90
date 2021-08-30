submodule(rmvs_classes) s_rmvs_discard
   use swiftest
contains  

   module subroutine rmvs_discard_tp(self, system, param)
      !! author: David A. Minton
      !!
      !! Check to see if test particles should be discarded based on pericenter passage distances with respect to planets encountered
      !!
      !! Adapted from Hal Levison's Swift routine discard_pl.f
      !! Adapted from Hal Levison's Swift routine rmvs_discard_pl.f90
      implicit none
      ! Arguments
      class(rmvs_tp),               intent(inout) :: self   !! RMVS test particle object
      class(swiftest_nbody_system), intent(inout) :: system !! Swiftest nbody system object
      class(swiftest_parameters),   intent(inout) :: param  !! Current run configuration parameters 
      ! Internals
      integer(I4B)                                 :: i
      character(len=STRMAX) :: timestr, idstri, idstrj

      if (self%nbody == 0) return

      associate(tp => self, ntp => self%nbody, pl => system%pl, t => param%t)
         do i = 1, ntp
            associate(iplperP => tp%plperP(i))
               if ((tp%status(i) == ACTIVE) .and. (tp%lperi(i))) then 
                  if ((tp%peri(i) < pl%radius(iplperP))) then
                     tp%status(i) = DISCARDED_PLQ
                     write(idstri, *) tp%id(i)
                     write(idstrj, *) pl%id(iplperP)
                     write(timestr, *) t
                     write(*, *) "Particle "  // trim(adjustl(tp%info(i)%name)) // " (" // trim(adjustl(idstri)) &
                              // ") q with respect to massive body " // trim(adjustl(pl%info(iplperP)%name)) // " (" // trim(adjustl(idstrj)) &
                              // ") is too small at t = " // trim(adjustl(timestr))
                     tp%ldiscard(i) = .true.
                     tp%lmask(i) = .false.
                     call tp%info(i)%set_value(status="DISCARDED_PLQ", discard_time=t, discard_xh=tp%xh(:,i), discard_vh=tp%vh(:,i), discard_body_id=pl%id(iplperP))
                  end if
               end if
            end associate
         end do
         ! Call the base method that this overrides
         call discard_tp(tp, system, param)
      end associate

   end subroutine rmvs_discard_tp

end submodule s_rmvs_discard