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
      class(swiftest_parameters),   intent(in)    :: param  !! Current run configuration parameters 
      ! Internals
      integer(I4B)                                 :: i

      associate(tp => self, ntp => self%nbody, pl => system%pl, t => param%t)
         do i = 1, ntp
            associate(iplperP => tp%plperP(i))
               if ((tp%status(i) == ACTIVE) .and. (tp%lperi(i))) then 
                  if ((tp%peri(i) < pl%radius(iplperP))) then
                     tp%status(i) = DISCARDED_PLQ
                     write(*, *) "Particle ",tp%id(i)," q with respect to Planet ",pl%id(iplperP)," is too small at t = ",t
                     tp%ldiscard(i) = .true.
                  end if
               end if
            end associate
         end do
         ! Call the base method that this overrides
         call discard_tp(tp, system, param)
      end associate

   end subroutine rmvs_discard_tp

end submodule s_rmvs_discard