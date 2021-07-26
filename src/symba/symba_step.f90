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
                  irec = -1
                  call pl%vh2vb(cb)
                  call pl%lindrift(cb, dth, mask=(pl%status(:) == ACTIVE), lbeg=.true.)
                  call pl%kick(system, param, t, dth, mask=(pl%status(:) == ACTIVE), lbeg=.true.)
                  call pl%drift(system, param, dt, mask=(pl%status(:) == ACTIVE .and. pl%levelg(:) == irec))

                  call tp%vh2vb(vbcb = -cb%ptbeg)
                  call tp%lindrift(cb, dth, mask=(tp%status(:) == ACTIVE), lbeg=.true.)
                  call tp%kick(system, param, t, dth, mask=(tp%status(:) == ACTIVE), lbeg=.true.)
                  call tp%drift(system, param, dt, mask=(tp%status(:) == ACTIVE .and. tp%levelg(:) == irec))

                  irec = 0
                  call system%recursive_step(param, irec)

                  call pl%kick(system, param, t, dth, mask=(pl%status(:) == ACTIVE), lbeg=.false.)
                  call pl%vb2vh(cb)
                  call pl%lindrift(cb, dth, mask=(pl%status(:) == ACTIVE), lbeg=.false.)

                  call tp%kick(system, param, t, dth, mask=(tp%status(:) == ACTIVE), lbeg=.true.)
                  call tp%vb2vh(vbcb = -cb%ptend)
                  call tp%lindrift(cb, dth, mask=(tp%status(:) == ACTIVE), lbeg=.false.)
               end select
            end select
         end select
      end associate
      return
   end subroutine symba_step_interp_system

   module recursive subroutine symba_step_recur_system(self, param, ireci)
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
      integer(I4B), value,        intent(in)    :: ireci !! input recursion level
      ! Internals
      integer(I4B) :: i, j, irecp, nloops, sgn
      real(DP) :: dtl, dth
      real(DP), dimension(NDIM) :: xr, vr
      logical :: lencounter

      associate(system => self, plplenc_list => self%plplenc_list, pltpenc_list => self%pltpenc_list)
         select type(pl => self%pl)
         class is (symba_pl)
            select type(tp => self%tp)
            class is (symba_tp)
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
                  sgn = 1
                  call plplenc_list%kick(system, dth, irecp, sgn)
                  call pltpenc_list%kick(system, dth, irecp, sgn)
                  if (ireci /= 0) then
                     sgn = -1
                     call plplenc_list%kick(system, dth, irecp, sgn)
                     call pltpenc_list%kick(system, dth, irecp, sgn)
                  end if

                  call pl%drift(system, param, dtl, mask=(pl%status(:) == ACTIVE .and. pl%levelg(:) == ireci))
                  call tp%drift(system, param, dtl, mask=(tp%status(:) == ACTIVE .and. tp%levelg(:) == ireci))

                  if (lencounter) call system%recursive_step(param, irecp)

                  sgn = 1
                  call plplenc_list%kick(system, dth, irecp, sgn)
                  call pltpenc_list%kick(system, dth, irecp, sgn)
                  if (ireci /= 0) then
                     sgn = -1
                     call plplenc_list%kick(system, dth, irecp, sgn)
                     call pltpenc_list%kick(system, dth, irecp, sgn)
                  end if
                  associate (plind1 => plplenc_list%index1(1:plplenc_list%nenc), &
                             plind2 => plplenc_list%index2(1:plplenc_list%nenc), &
                             plind3 => pltpenc_list%index1(1:pltpenc_list%nenc), &
                             tpind  => pltpenc_list%index2(1:pltpenc_list%nenc))
                     where(pl%levelg([plind1,plind2,plind3]) == irecp) pl%levelg(:) = ireci
                     where(tp%levelg(tpind) == irecp) tp%levelg(:) = ireci
                  end associate
               end do
            end select
         end select
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
