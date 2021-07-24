submodule (helio_classes) s_helio_getacch
   use swiftest
contains
   module subroutine helio_getacch_pl(self, system, param, t, lbeg)
      !! author: David A. Minton
      !!
      !! Compute heliocentric accelerations of massive bodies
      !!
      !! Adapted from David E. Kaufmann's Swifter routine helio_getacch.f90
      !! Adapted from Hal Levison's Swift routine helio_getacch.f
      implicit none
      ! Arguments
      class(helio_pl),              intent(inout) :: self   !! Helio massive body particle data structure
      class(swiftest_nbody_system), intent(inout) :: system !! WHM nbody system object
      class(swiftest_parameters),   intent(in)    :: param  !! Current run configuration parameters 
      real(DP),                     intent(in)    :: t      !! Current simulation time
      logical, optional,            intent(in)    :: lbeg   !! Optional argument that determines whether or not this is the beginning or end of the step

      associate(cb => system%cb, pl => self, npl => self%nbody)
         pl%ah(:,:) = 0.0_DP
         call pl%accel_int()
         if (param%loblatecb) then 
            cb%aoblbeg = cb%aobl
            call pl%accel_obl(system)
            cb%aoblend = cb%aobl
            if (param%ltides) then
               cb%atidebeg = cb%atide
               call pl%accel_tides(system)
               cb%atideend = cb%atide
            end if
         end if
         if (param%lextra_force) call pl%accel_user(system, param, t)
         !if (param%lgr) call pl%gr_accel(param)
      end associate

      return
      end subroutine helio_getacch_pl

      module subroutine helio_getacch_tp(self, system, param, t, lbeg)
         !! author: David A. Minton
         !!
         !! Compute heliocentric accelerations of test particles
         !!
         !! Adapted from David E. Kaufmann's Swifter routine helio_getacch_tp.f90
         !! Adapted from Hal Levison's Swift routine helio_getacch_tp.f
         implicit none
         ! Arguments
         class(helio_tp),              intent(inout) :: self   !! WHM test particle data structure
         class(swiftest_nbody_system), intent(inout) :: system !! WHM nbody system object
         class(swiftest_parameters),   intent(in)    :: param  !! Current run configuration parameters 
         real(DP),                     intent(in)    :: t      !! Current time
         logical, optional,            intent(in)    :: lbeg   !! Optional argument that determines whether or not this is the beginning or end of the step
      
         associate(tp => self, cb => system%cb, pl => system%pl, npl => system%pl%nbody)
            tp%ah(:,:) = 0.0_DP
            if (present(lbeg)) system%lbeg = lbeg
            if (system%lbeg) then
               call tp%accel_int(pl%Gmass(:), pl%xbeg(:,:), npl)
            else
               call tp%accel_int(pl%Gmass(:), pl%xend(:,:), npl)
            end if
            if (param%loblatecb) call tp%accel_obl(system)
            if (param%lextra_force) call tp%accel_user(system, param, t)
            !if (param%lgr) call tp%gr_accel(param)
         end associate
         return
      end subroutine helio_getacch_tp

end submodule s_helio_getacch
