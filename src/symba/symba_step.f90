!! Copyright 2022 - David Minton, Carlisle Wishard, Jennifer Pouplin, Jake Elliott, & Dana Singh
!! This file is part of Swiftest.
!! Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
!! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!! Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
!! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!! You should have received a copy of the GNU General Public License along with Swiftest. 
!! If not, see: https://www.gnu.org/licenses. 

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
     
      select type(pl => self%pl)
      class is (symba_pl)
         select type(tp => self%tp)
         class is (symba_tp)
            select type(param)
            class is (symba_parameters)
               call self%reset(param)
               lencounter = pl%encounter_check(param, self, dt, 0) .or. tp%encounter_check(param, self, dt, 0)
               if (lencounter) then
                  call self%interp(param, t, dt)
                  !call self%encounter_history%dump(param) 
               else
                  self%irec = -1
                  call helio_step_system(self, param, t, dt)
               end if
               param%lfirstkick = pl%lfirst
            end select
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
                  if (pl%lfirst) call pl%vh2vb(cb)
                  call pl%lindrift(cb, dth, lbeg=.true.)
                  call pl%kick(system, param, t, dth, lbeg=.true.)
                  if (param%lgr) call pl%gr_pos_kick(system, param, dth)
                  call pl%drift(system, param, dt)

                  if (tp%nbody > 0) then
                     if (tp%lfirst) call tp%vh2vb(vbcb = -cb%ptbeg)
                     call tp%lindrift(cb, dth, lbeg=.true.)
                     call tp%kick(system, param, t, dth, lbeg=.true.)
                     if (param%lgr) call tp%gr_pos_kick(system, param, dth)
                     call tp%drift(system, param, dt)
                  end if

                  call system%recursive_step(param, t, 0)
                  system%irec = -1

                  if (param%lgr) call pl%gr_pos_kick(system, param, dth)
                  call pl%kick(system, param, t, dth, lbeg=.false.)
                  call pl%lindrift(cb, dth, lbeg=.false.)
                  call pl%vb2vh(cb)

                  if (tp%nbody > 0) then
                     if (param%lgr) call tp%gr_pos_kick(system, param, dth)
                     call tp%kick(system, param, t, dth, lbeg=.false.)
                     call tp%lindrift(cb, dth, lbeg=.false.)
                     call tp%vb2vh(vbcb = -cb%ptend)
                  end if
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
      integer(I4B) :: irecp

      associate(system => self, plplenc_list => self%plplenc_list, pltpenc_list => self%pltpenc_list, &
                npl => self%pl%nbody, ntp => self%tp%nbody)
         select type(pl => self%pl)
         class is (symba_pl)
            select type(tp => self%tp)
            class is (symba_tp)
               irecp = ireci + 1

               if (npl >0) where(pl%levelg(1:npl) == irecp) pl%levelg(1:npl) = ireci
               if (ntp > 0) where(tp%levelg(1:ntp) == irecp) tp%levelg(1:ntp) = ireci
               if (plplenc_list%nenc > 0) then
                  where(plplenc_list%level(1:plplenc_list%nenc) == irecp) 
                     plplenc_list%level(1:plplenc_list%nenc) = ireci
                  endwhere
               end if
               if (pltpenc_list%nenc > 0) then
                  where(pltpenc_list%level(1:pltpenc_list%nenc) == irecp) 
                     pltpenc_list%level(1:pltpenc_list%nenc) = ireci
                  endwhere
               end if

               system%irec = ireci

            end select
         end select
      end associate

      return
   end subroutine symba_step_set_recur_levels_system


   recursive module subroutine symba_step_recur_system(self, param, t, ireci)
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
      real(DP),                   intent(in)    :: t
      integer(I4B),               intent(in)    :: ireci !! input recursion level
      ! Internals
      integer(I4B) :: j, irecp, nloops
      real(DP) :: dtl, dth
      logical :: lencounter, lplpl_collision, lpltp_collision

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
                  lencounter = plplenc_list%encounter_check(param, system, dtl, irecp) &
                          .or. pltpenc_list%encounter_check(param, system, dtl, irecp)
                   
                  call plplenc_list%kick(system, dth, irecp, 1)
                  call pltpenc_list%kick(system, dth, irecp, 1)
                  if (ireci /= 0) then
                     call plplenc_list%kick(system, dth, irecp, -1)
                     call pltpenc_list%kick(system, dth, irecp, -1)
                  end if

                  if (param%lgr) then
                     call pl%gr_pos_kick(system, param, dth)
                     call tp%gr_pos_kick(system, param, dth)
                  end if

                  call pl%drift(system, param, dtl)
                  call tp%drift(system, param, dtl)

                  !call system%

                  if (lencounter) call system%recursive_step(param, t+dth,irecp)
                  system%irec = ireci

                  if (param%lgr) then
                     call pl%gr_pos_kick(system, param, dth)
                     call tp%gr_pos_kick(system, param, dth)
                  end if

                  call plplenc_list%kick(system, dth, irecp, 1)
                  call pltpenc_list%kick(system, dth, irecp, 1)
                  if (ireci /= 0) then
                     call plplenc_list%kick(system, dth, irecp, -1)
                     call pltpenc_list%kick(system, dth, irecp, -1)
                  end if

                  if (param%lclose) then
                     lplpl_collision = plplenc_list%collision_check(system, param, t+dtl, dtl, ireci) 
                     lpltp_collision = pltpenc_list%collision_check(system, param, t+dtl, dtl, ireci) 

                     if (lplpl_collision) call plplenc_list%resolve_collision(system, param, t+dtl, dtl, ireci)
                     if (lpltp_collision) call pltpenc_list%resolve_collision(system, param, t+dtl, dtl, ireci)
                  end if

                  call self%set_recur_levels(ireci)

               end do
            end select
         end select
      end associate

      return
   end subroutine symba_step_recur_system


   module subroutine symba_step_reset_system(self, param)
      !! author: David A. Minton
      !!
      !! Resets pl, tp,and encounter structures at the start of a new step
      !! 
      !! Adapted from David E. Kaufmann's Swifter routine: symba_step.f90
      !! Adapted from Hal Levison's Swift routine symba5_step.f
      implicit none
      ! Arguments
      class(symba_nbody_system), intent(inout) :: self  !! SyMBA nbody system object
      class(symba_parameters),   intent(in)    :: param !! Current run configuration parameters with SyMBA additions
      ! Internals
      integer(I4B) :: i
      integer(I8B) :: nenc_old

      associate(system => self)
         select type(pl => system%pl)
         class is (symba_pl)
            select type(tp => system%tp)
            class is (symba_tp)
               associate(npl => pl%nbody, ntp => tp%nbody)
                  nenc_old = system%plplenc_list%nenc
                  call system%plplenc_list%setup(0_I8B)
                  call system%plplcollision_list%setup(0_I8B)
                  system%ienc_frame = 0
                  if (allocated(system%encounter_history)) deallocate(system%encounter_history)
                  if (npl > 0) then
                     pl%lcollision(1:npl) = .false.
                     call pl%reset_kinship([(i, i=1, npl)])
                     pl%nplenc(1:npl) = 0
                     pl%ntpenc(1:npl) = 0
                     pl%levelg(1:npl) = -1
                     pl%levelm(1:npl) = -1
                     pl%lencounter(1:npl) = .false.
                     pl%lcollision(1:npl) = .false.
                     pl%ldiscard(1:npl) = .false.
                     pl%lmask(1:npl) = .true.
                     call pl%set_renc(0)
                     call system%plplenc_list%setup(nenc_old) ! This resizes the pl-pl encounter list to be the same size as it was the last step, to decrease the number of potential resize operations that have to be one inside the step
                     system%plplenc_list%nenc = 0 ! Sets the true number of encounters back to 0 after resizing
                  end if
            
                  nenc_old = system%pltpenc_list%nenc
                  call system%pltpenc_list%setup(0_I8B)
                  if (ntp > 0) then
                     tp%nplenc(1:ntp) = 0 
                     tp%levelg(1:ntp) = -1
                     tp%levelm(1:ntp) = -1
                     tp%lmask(1:ntp) = .true.
                     tp%ldiscard(1:ntp) = .false.
                     call system%pltpenc_list%setup(nenc_old)! This resizes the pl-tp encounter list to be the same size as it was the last step, to decrease the number of potential resize operations that have to be one inside the step
                     system%pltpenc_list%nenc = 0 ! Sets the true number of encounters back to 0 after resizing
                  end if

                  call system%pl_adds%setup(0, param)
                  call system%pl_discards%setup(0, param)

                  tp%lfirst = param%lfirstkick
                  pl%lfirst = param%lfirstkick
               end associate
            end select
         end select
      end associate

      return
   end subroutine symba_step_reset_system

end submodule s_symba_step
