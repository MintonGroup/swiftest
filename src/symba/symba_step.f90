submodule (symba_classes) s_symba_step
   use swiftest
contains
   module subroutine symba_step_system(self, param, t, dt)
      !! author: David A. Minton
      !!
      !! Step planets and active test particles ahead in democratic heliocentric coordinates, descending the recursive
      !!   branch if necessary to handle possible close encounters
      !!
      !! Adapted from David E. Kaufmann's Swifter routine: symba_step.f90
      !! Adapted from Hal Levison's Swift routine symba5_step_pl.f
      implicit none
      ! Arguments
      class(symba_nbody_system),  intent(inout) :: self   !! SyMBA nbody system object
      class(swiftest_parameters), intent(inout) :: param  !! Current run configuration parameters
      real(DP),                   intent(in)    :: t      !! Simulation time
      real(DP),                   intent(in)    :: dt     !! Current stepsize
      ! Internals
      logical           :: lencounter_pl, lencounter_tp, lencounter
     
      select type(pl => self%pl)
      class is (symba_pl)
         select type(tp => self%tp)
         class is (symba_tp)
            lencounter = pl%encounter_check(self, dt) .or. tp%encounter_check(self, dt)
            if (lencounter) then
               call self%interp(param, t, dt)
            else
               call helio_step_system(self, param, t, dt)
            end if
         end select
      end select

      return

   end subroutine symba_step_system

   module subroutine symba_step_interp_system(self, param, t, dt)
      !! author: David A. Minton
      !!
      !! Step planets and active test particles ahead in democratic heliocentric coordinates, calling the recursive
      !!         subroutine to descend to the appropriate level to handle close encounters
      !!
      !! Adapted from David E. Kaufmann's Swifter routine: symba_step_interp.f90
      !! Adapted from Hal Levison's Swift routine symba5_step_interp.f
      implicit none
      ! Arguments
      class(symba_nbody_system),  intent(inout) :: self   !! SyMBA nbody system object
      class(swiftest_parameters), intent(inout) :: param  !! Current run configuration parameters 
      real(DP),                   intent(in)    :: t      !! Simulation time
      real(DP),                   intent(in)    :: dt     !! Current stepsize
      ! Internals 
      real(DP)                                  :: dth

      dth = 0.5_DP * dt
      associate(system => self)
         select type(pl => system%pl)
         class is (symba_pl)
            select type(tp => system%tp)
            class is (symba_tp)
               select type(cb => system%cb)
               class is (symba_cb)
                  call pl%vh2vb(cb)
                  call pl%lindrift(cb, dth, lbeg=.true.)
                  call tp%vh2vb(vbcb = -cb%ptbeg)
                  call tp%lindrift(cb, dth, lbeg=.true.)

                  call pl%set_beg_end(xbeg = pl%xh)
                  call pl%accel(system, param, t)
                  call tp%accel(system, param, t, lbeg=.true.)

                  call pl%kick(dth)
                  call tp%kick(dth)

                  system%irec = -1
                  call pl%drift(system, param, dt)
                  call tp%drift(system, param, dt)
                  system%irec = 0

                  call system%recursive_step(param, t, dt)

                  call pl%set_beg_end(xend = pl%xh)
                  call pl%accel(system, param, t + dt)
                  call tp%accel(system, param, t + dt, lbeg=.false.)

                  call pl%kick(dth)
                  call tp%kick(dth)

                  call pl%vh2vb(cb)
                  call pl%lindrift(cb, dth, lbeg=.false.)
                  call tp%vh2vb(vbcb = -cb%ptend)
                  call tp%lindrift(cb, dth, lbeg=.false.)
               end select
            end select
         end select
      end associate
      return
   end subroutine symba_step_interp_system

   module recursive subroutine symba_step_recur_system(self, param, t, dt)
      !! author: David A. Minton
      !!
      !! Step interacting planets and active test particles ahead in democratic heliocentric coordinates at the current
      !!         recursion level, if applicable, and descend to the next deeper level if necessarys
      !!
      !! Adapted from David E. Kaufmann's Swifter routine: symba_step_recur.f90
      !! Adapted from Hal Levison's Swift routine symba5_step_recur.f
      implicit none
      ! Arguments
      class(symba_nbody_system),  intent(inout) :: self  !! SyMBA nbody system object
      class(swiftest_parameters), intent(inout) :: param !! Current run configuration parameters 
      real(DP),                   intent(in)    :: t     !! Simulation time
      real(DP),                   intent(in)    :: dt    !! Current stepsize
      ! Internals
   end subroutine symba_step_recur_system

end submodule s_symba_step
