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
     
      call self%reset()
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
      class(symba_nbody_system),  intent(inout) :: self  !! SyMBA nbody system object
      class(swiftest_parameters), intent(inout) :: param !! Current run configuration parameters 
      real(DP),                   intent(in)    :: t     !! Simulation time
      real(DP),                   intent(in)    :: dt    !! Current stepsize
      ! Internals
      real(DP)                                  :: dth   !! Half step size
      integer(I4B)                              :: irec  !! Recursion level                 

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

                  call pl%drift(system, param, dt, pl%status(:) == ACTIVE)
                  call tp%drift(system, param, dt, tp%status(:) == ACTIVE)
                  irec = 0
                  call system%recursive_step(param, t, dt, irec)

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

   module recursive subroutine symba_step_recur_system(self, param, t, dt, ireci)
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
      integer(I4B), value,        intent(in)    :: ireci !! input recursion level
      ! Internals
      integer(I4B) :: i, j, irecp, icflg, index_i, index_j, index_pl, index_tp
      real(DP) :: dtl, dth,sgn

      associate(plplenc_list => self%plplenc_list, pltpenc_list => self%pltpenc_list)
         dtl = param%dt / (NTENC**ireci)
         dth = 0.5_DP * dtl
         IF (dtl / param%dt < VSMALL) THEN
            WRITE(*, *) "SWIFTEST Warning:"
            WRITE(*, *) "   In symba_step_recur_system, local time step is too small"
            WRITE(*, *) "   Roundoff error will be important!"
            call util_exit(FAILURE)
         END IF
         irecp = ireci + 1
         if (ireci == 0) then
            icflg = 0
            
         end if
      end associate

   end subroutine symba_step_recur_system

   module subroutine symba_step_reset_system(self)
      !! author: David A. Minton
      !!
      !! Resets pl, tp,and encounter structures at the start of a new step
      !! 
      implicit none
      ! Arguments
      class(symba_nbody_system),  intent(inout) :: self !! SyMBA nbody system object
      ! Internals
      integer(I4B) :: i

      associate(system => self, pltpenc_list => self%pltpenc_list, plplenc_list => self%plplenc_list, mergeadd_list => self%mergeadd_list, mergesub_list => self%mergesub_list)
         select type(pl => system%pl)
         class is (symba_pl)
            select type(tp => system%tp)
            class is (symba_tp)
               pl%lcollision(:) = .false.
               pl%kin(:)%parent = [(i, i=1, pl%nbody)]
               pl%kin(:)%nchild = 0
               do i = 1, pl%nbody
                  if (allocated(pl%kin(i)%child)) deallocate(pl%kin(i)%child)
               end do
               pl%nplenc(:) = 0
               pl%ntpenc(:) = 0
               pl%levelg(:) = 0
               pl%levelm(:) = 0
               pl%lencounter = .false.
               pl%lcollision = .false.
            
               tp%nplenc(:) = 0 
               tp%levelg(:) = 0
               tp%levelm(:) = 0
            
               plplenc_list%nenc = 0
               pltpenc_list%nenc = 0

               mergeadd_list%nbody = 0
               mergesub_list%nbody = 0
            end select
         end select
      end associate

   
   end subroutine symba_step_reset_system



end submodule s_symba_step
