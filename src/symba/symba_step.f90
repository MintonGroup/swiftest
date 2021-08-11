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
      logical :: lencounter
     
      call self%reset()
      select type(pl => self%pl)
      class is (symba_pl)
         select type(tp => self%tp)
         class is (symba_tp)
            lencounter = pl%encounter_check(self, dt, 0) .or. tp%encounter_check(self, dt, 0)
            if (lencounter) then
               tp%lfirst = pl%lfirst
               call self%interp(param, t, dt)
               pl%lfirst = .true.
               tp%lfirst = .true.
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

      dth = 0.5_DP * dt
      associate(system => self)
         select type(pl => system%pl)
         class is (symba_pl)
            select type(tp => system%tp)
            class is (symba_tp)
               select type(cb => system%cb)
               class is (symba_cb)
                  system%irec = -1
                  call pl%vh2vb(cb)
                  call pl%lindrift(cb, dth, lbeg=.true.)
                  call pl%kick(system, param, t, dth, lbeg=.true.)
                  call pl%drift(system, param, dt)

                  call tp%vh2vb(vbcb = -cb%ptbeg)
                  call tp%lindrift(cb, dth, lbeg=.true.)
                  call tp%kick(system, param, t, dth, lbeg=.true.)
                  call tp%drift(system, param, dt)

                  call system%recursive_step(param, t, 0)

                  call pl%kick(system, param, t, dth, lbeg=.false.)
                  call pl%vb2vh(cb)
                  call pl%lindrift(cb, dth, lbeg=.false.)

                  call tp%kick(system, param, t, dth, lbeg=.false.)
                  call tp%vb2vh(vbcb = -cb%ptend)
                  call tp%lindrift(cb, dth, lbeg=.false.)
               end select
            end select
         end select
      end associate

      return
   end subroutine symba_step_interp_system


   module subroutine symba_step_set_recur_levels_system(self, ireci)
      !! author: David A. Minton
      !!
      !! Resets pl, tp,and encounter structures at the start of a new step
      !! 
      !! Adapted from David E. Kaufmann's Swifter routine: symba_step_recur.f90
      !! Adapted from Hal Levison's Swift routine symba5_step_recur.f
      implicit none
      ! Arguments
      class(symba_nbody_system),  intent(inout) :: self  !! SyMBA nbody system object
      integer(I4B),               intent(in)    :: ireci !! Input recursion level 
      ! Internals
      integer(I4B) :: k, irecp

      associate(system => self, plplenc_list => self%plplenc_list, pltpenc_list => self%pltpenc_list)
         select type(pl => self%pl)
         class is (symba_pl)
            select type(tp => self%tp)
            class is (symba_tp)
               irecp = ireci + 1

               if (plplenc_list%nenc > 0) then
                  do k = 1, plplenc_list%nenc
                     associate(i => plplenc_list%index1(k), j => plplenc_list%index2(k))
                        if (pl%levelg(i) == irecp) pl%levelg(i) = ireci
                        if (pl%levelg(j) == irecp) pl%levelg(j) = ireci
                     end associate
                  end do
                  where(plplenc_list%level(1:plplenc_list%nenc) == irecp) plplenc_list%level(1:plplenc_list%nenc) = ireci
               end if

               if (pltpenc_list%nenc > 0) then
                  do k = 1, pltpenc_list%nenc
                     associate(i => pltpenc_list%index1(k), j => pltpenc_list%index2(k))
                        if (pl%levelg(i) == irecp) pl%levelg(i) = ireci
                        if (tp%levelg(j) == irecp) tp%levelg(j) = ireci
                     end associate
                  end do
                  where(pltpenc_list%level(1:pltpenc_list%nenc) == irecp) pltpenc_list%level(1:pltpenc_list%nenc) = ireci
               end if

               system%irec = ireci

            end select
         end select
      end associate

      return
   end subroutine symba_step_set_recur_levels_system


   module recursive subroutine symba_step_recur_system(self, param, t, ireci)
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
      real(DP),                   value         :: t
      integer(I4B),               value         :: ireci !! input recursion level
      ! Internals
      integer(I4B) :: i, j, irecp, nloops
      real(DP) :: dtl, dth
      real(DP), dimension(NDIM) :: xr, vr
      logical :: lencounter

      associate(system => self, plplenc_list => self%plplenc_list, pltpenc_list => self%pltpenc_list)
         select type(pl => self%pl)
         class is (symba_pl)
            select type(tp => self%tp)
            class is (symba_tp)
               system%irec = ireci
               dtl = param%dt / (NTENC**ireci)
               dth = 0.5_DP * dtl
               IF (dtl / param%dt < VSMALL) THEN
                  write(*, *) "SWIFTEST Warning:"
                  write(*, *) "   In symba_step_recur_system, local time step is too small"
                  write(*, *) "   Roundoff error will be important!"
                  call util_exit(FAILURE)
               END IF
               irecp = ireci + 1
               if (ireci == 0) then
                  nloops = 1
               else
                  nloops = NTENC
               end if
               do j = 1, nloops
                  lencounter = plplenc_list%encounter_check(system, dtl, irecp) .or. pltpenc_list%encounter_check(system, dtl, irecp)
                   
                  call plplenc_list%kick(system, dth, irecp, 1)
                  call pltpenc_list%kick(system, dth, irecp, 1)
                  if (ireci /= 0) then
                     call plplenc_list%kick(system, dth, irecp, -1)
                     call pltpenc_list%kick(system, dth, irecp, -1)
                  end if

                  call pl%drift(system, param, dtl)
                  call tp%drift(system, param, dtl)

                  if (lencounter) call system%recursive_step(param, t+dth,irecp)
                  system%irec = ireci

                  call plplenc_list%kick(system, dth, irecp, 1)
                  call pltpenc_list%kick(system, dth, irecp, 1)
                  if (ireci /= 0) then
                     call plplenc_list%kick(system, dth, irecp, -1)
                     call pltpenc_list%kick(system, dth, irecp, -1)
                  end if

                  if (param%lclose) then
                     call plplenc_list%collision_check(system, param, t+dtl, dtl, ireci) 
                     call pltpenc_list%collision_check(system, param, t+dtl, dtl, ireci) 
                  end if

                  call self%set_recur_levels(ireci)

               end do
            end select
         end select
      end associate

      return
   end subroutine symba_step_recur_system


   module subroutine symba_step_reset_system(self)
      !! author: David A. Minton
      !!
      !! Resets pl, tp,and encounter structures at the start of a new step
      !! 
      !! Adapted from David E. Kaufmann's Swifter routine: symba_step.f90
      !! Adapted from Hal Levison's Swift routine symba5_step.f
      implicit none
      ! Arguments
      class(symba_nbody_system),  intent(inout) :: self !! SyMBA nbody system object
      ! Internals
      integer(I4B) :: i

      associate(system => self)
         select type(pl => system%pl)
         class is (symba_pl)
            select type(tp => system%tp)
            class is (symba_tp)
               if (pl%nbody > 0) then
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
                  system%plplenc_list%nenc = 0
                  system%plplcollision_list%nenc = 0
               end if
           
               if (tp%nbody > 0) then
                  tp%nplenc(:) = 0 
                  tp%levelg(:) = 0
                  tp%levelm(:) = 0
               system%pltpenc_list%nenc = 0
               end if

               call system%pl_adds%resize(0)
               call system%pl_discards%resize(0)
            end select
         end select
      end associate

      return
   end subroutine symba_step_reset_system

end submodule s_symba_step
