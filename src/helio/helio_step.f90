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
                  call system%set_beg_end(xbeg = pl%xh)
                  call pl%step(system, param, t, dt)
                  if (tp%nbody > 0) then
                     call system%set_beg_end(xend = pl%xh)
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
  
      associate(pl => self, cb => system%cb)
         dth = 0.5_DP * dt
         if (lfirst) then
            call pl%vh2vb(cb)
            lfirst = .false.
         end if
         call pl%lindrift(cb, dth, ptbeg)
         call pl%getacch(system, param, t)
         call pl%kickvb(dth)

         call pl%drift(system, param, dt)
         call pl%getacch(system, param, t + dt)
         call pl%kickvb(dth)
         call pl%lindrift(cb, dth, ptend)
         call pl%vb2vh(cb)
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
   
      select type(system)
      class is (helio_nbody_system)
         associate(tp => self, cb => system%cb, pl => system%pl, xbeg => system%xbeg, xend => system%xend)
            dth = 0.5_DP * dt
            if (lfirst) then
               call tp%vh2vb(vbcb = -tp%ptbeg)
               lfirst = .false.
            end if
            call tp%lindrift(dth, tp%ptbeg)
            call tp%getacch(system, param, t, xbeg)
            call tp%kickvb(dth)
            call tp%drift(system, param, dt)
            call tp%getacch(system, param, t + dt, xend)
            call tp%kickvb(dth)
            call tp%lindrift(dth, tp%ptend)
            call tp%vb2vh(vbcb = -tp%ptend)
         end associate
      end select
   
      return
   
   end subroutine helio_step_tp

end submodule s_helio_step
