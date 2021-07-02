submodule(helio_classes) s_helio_step
   use swiftest
contains
   module subroutine helio_step_system(self, param, t, dt)
      !! author: David A. Minton
      !!
      !! Step massive bodies and and active test particles ahead in heliocentric coordinates
      !! 
      !! Adapted from Hal Levison's Swift routine step_kdk.f
      !! Adapted from David E. Kaufmann's Swifter routine helio_step.f90
      implicit none
      ! Arguments
      class(helio_nbody_system),  intent(inout) :: self   !! Helio nbody system object
      class(swiftest_parameters), intent(inout) :: param  !! Current run configuration parameters
      real(DP),                   intent(in)    :: t      !! Simulation time
      real(DP),                   intent(in)    :: dt     !! Current stepsize

      select type(system => self)
      class is (helio_nbody_system)
      select type(cb => self%cb)
      class is (helio_cb)
      select type(pl => self%pl)
      class is (helio_pl)
      select type(tp => self%tp)
      class is (helio_tp)
         call pl%set_rhill(cb)
         call tp%set_beg_end(xbeg = pl%xh)
         call pl%step(system, param, t, dt)
         if (ntp > 0) then
            call tp%set_beg_end(xend = pl%xh)
            call tp%step(system, param, t, dt)
         end if
      end select
      end select
      end select
      end select
      return
   end subroutine helio_step_system 

   module subroutine helio_step_pl(self, system, param, t, dt)
      !! author: David A. Minton
      !!
      !! Step massive bodies ahead Democratic Heliocentric method
      !!
      !! Adapted from David E. Kaufmann's Swifter helio_step_pl.f90
      !! Adapted from Hal Levison's Swift routine helio_step_pl.f
      implicit none
      ! Arguments
      class(helio_pl),              intent(inout) :: self   !! Helio massive body particle data structure
      class(swiftest_nbody_system), intent(inout) :: system !! Swiftest nboody system
      class(swiftest_parameters),   intent(inout) :: param  !! Current run configuration parameters 
      real(DP),                     intent(in)    :: t      !! Current simulation time
      real(DP),                     intent(in)    :: dt     !! Stepsize
      ! Internals 
      integer(I4B)     :: i
      real(DP)         :: dth, msys
      real(DP), dimension(NDIM) :: ptbeg, ptend !! TODO: Incorporate these into the tp structure
      logical, save :: lfirst = .true.
  
      associate(cb => system%cb)
         dth = 0.5_DP * dt
         if (lfirst) then
            call self%vh2vb(cb)
            lfirst = .false.
         end if
         call self%lindrift(cb, dth, ptbeg)
         call self%getacch(system, param, t)
         call self%kickvb(dth)

         call self%drift(system, param, dt)
         call self%getacch(system, param, t + dt)
         call self%kickvb(dth)
         call self%lindrift(cb, dth, ptend)
         call self%vb2vh(cb)
      end associate
   
      return
   
   end subroutine helio_step_pl

   module subroutine helio_step_tp(self, system, param, t, dt)
      !! author: David A. Minton
      !!
      !! Step active test particles ahead using Democratic Heliocentric method
      !!
      !! Adapted from David E. Kaufmann's Swifter routine helio_step_tp.f90
      !! Adapted from Hal Levison's Swift routine helio_step_tp.f
      implicit none
      ! Arguments
      class(helio_tp),              intent(inout) :: self    !! Helio test particle data structure
      class(swiftest_nbody_system), intent(inout) :: system !! Swiftest nboody system
      class(swiftest_parameters),   intent(inout) :: param  !! Current run configuration parameters 
      real(DP),                     intent(in)    :: t      !! Current simulation time
      real(DP),                     intent(in)    :: dt     !! Stepsiz
      ! Internals
      logical, save  :: lfirst = .true. !! Flag to indicate that this is the first call
      real(DP) :: dth   !! Half step size
   
   ! executable code
      associate(cb => system%cb, pl => system%pl)
         dth = 0.5_DP * dt
         if (lfirst) then
            call self%vh2vb(vbcb = -self%ptbeg)
            lfirst = .false.
         end if
         call self%lindrift(dth, self%ptbeg)
         call self%getacch(system, param, t, self%xbeg)
         call self%kickvb(dth)
         call self%drift(system, param, dt)
         call self%getacch(system, param, t + dt, self%xend)
         call self%kickvb(dth)
         call self%lindrift(dth, self%ptend)
         call self%vb2vh(vbcb = -self%ptend)
      end associate
   
      return
   
   end subroutine helio_step_tp

end submodule s_helio_step
