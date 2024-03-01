! Copyight 2022 - David Minton, Carlisle Wishard, Jennifer Pouplin, Jake Elliott, & Dana Singh
! This file is part of Swiftest.
! Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
! Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
! You should have received a copy of the GNU General Public License along with Swiftest. 
! If not, see: https://www.gnu.org/licenses. 

submodule (symba) s_symba_step
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
      class is (swiftest_parameters)
         call self%reset(param)
         lencounter = pl%encounter_check(param, self, dt, 0) .or. tp%encounter_check(param, self, dt, 0)
         if (lencounter) then
            if (param%lenc_save_trajectory) call self%encounter_history%take_snapshot(param, self, t, "trajectory") 
            call self%interp(param, t, dt)
            if (param%lenc_save_trajectory) call self%encounter_history%take_snapshot(param, self, t+dt, "trajectory") 
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

      select type(pl => self%pl)
      class is (symba_pl)
      select type(tp => self%tp)
      class is (symba_tp)
      select type(cb => self%cb)
      class is (symba_cb)
         associate(nbody_system => self)
            dth = 0.5_DP * dt
            nbody_system%irec = -1
            if (pl%lfirst) call pl%vh2vb(cb)
            call pl%lindrift(cb, dth, lbeg=.true.)
            call pl%kick(nbody_system, param, t, dth, lbeg=.true.)
            if (param%lgr) call pl%gr_pos_kick(nbody_system, param, dth)
            call pl%drift(nbody_system, param, dt)

            if (tp%nbody > 0) then
               if (tp%lfirst) call tp%vh2vb(vbcb = -cb%ptbeg)
               call tp%lindrift(cb, dth, lbeg=.true.)
               call tp%kick(nbody_system, param, t, dth, lbeg=.true.)
               if (param%lgr) call tp%gr_pos_kick(nbody_system, param, dth)
               call tp%drift(nbody_system, param, dt)
            end if

            call nbody_system%recursive_step(param, t, 0)
            nbody_system%irec = -1

            if (param%lgr) call pl%gr_pos_kick(nbody_system, param, dth)
            call pl%kick(nbody_system, param, t, dth, lbeg=.false.)
            call pl%lindrift(cb, dth, lbeg=.false.)
            call pl%vb2vh(cb)

            if (tp%nbody > 0) then
               if (param%lgr) call tp%gr_pos_kick(nbody_system, param, dth)
               call tp%kick(nbody_system, param, t, dth, lbeg=.false.)
               call tp%lindrift(cb, dth, lbeg=.false.)
               call tp%vb2vh(vbcb = -cb%ptend)
            end if
         end associate
      end select
      end select
      end select

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

      select type(pl => self%pl)
      class is (symba_pl)
      select type(tp => self%tp)
      class is (symba_tp)
         associate(nbody_system => self, plpl_encounter => self%plpl_encounter, pltp_encounter => self%pltp_encounter, npl => self%pl%nbody, ntp => self%tp%nbody)

            irecp = ireci + 1

            if (npl >0) where(pl%levelg(1:npl) == irecp) pl%levelg(1:npl) = ireci
            if (ntp > 0) where(tp%levelg(1:ntp) == irecp) tp%levelg(1:ntp) = ireci
            if (plpl_encounter%nenc > 0) then
               where(plpl_encounter%level(1:plpl_encounter%nenc) == irecp) 
                  plpl_encounter%level(1:plpl_encounter%nenc) = ireci
               endwhere
            end if
            if (pltp_encounter%nenc > 0) then
               where(pltp_encounter%level(1:pltp_encounter%nenc) == irecp) 
                  pltp_encounter%level(1:pltp_encounter%nenc) = ireci
               endwhere
            end if

            nbody_system%irec = ireci

         end associate
      end select
      end select

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
      logical :: lencounter

      select type(param)
      class is (swiftest_parameters)
      select type(pl => self%pl)
      class is (symba_pl)
      select type(tp => self%tp)
      class is (symba_tp)
      select type(plpl_encounter => self%plpl_encounter)
      class is (symba_list_plpl)
      select type(pltp_encounter => self%pltp_encounter)
      class is (symba_list_pltp)
         associate(nbody_system => self, lplpl_collision => plpl_encounter%lcollision, lpltp_collision => pltp_encounter%lcollision)
            nbody_system%irec = ireci
            dtl = param%dt / (NTENC**ireci)
            dth = 0.5_DP * dtl
            IF (dtl / param%dt < epsilon(1.0_DP)) THEN
               write(*, *) "SWIFTEST Warning:"
               write(*, *) "   In symba_step_recur_system, local time step is too small"
               write(*, *) "   Roundoff error will be important!"
               call base_util_exit(FAILURE,param%display_unit)
            END IF
            irecp = ireci + 1
            if (ireci == 0) then
               nloops = 1
            else
               nloops = NTENC
            end if
            do j = 1, nloops
               lencounter = plpl_encounter%encounter_check(param, nbody_system, dtl, irecp) &
                     .or. pltp_encounter%encounter_check(param, nbody_system, dtl, irecp)
               
               call plpl_encounter%kick(nbody_system, dth, irecp, 1)
               call pltp_encounter%kick(nbody_system, dth, irecp, 1)
               if (ireci /= 0) then
                  call plpl_encounter%kick(nbody_system, dth, irecp, -1)
                  call pltp_encounter%kick(nbody_system, dth, irecp, -1)
               end if

               if (param%lgr) then
                  call pl%gr_pos_kick(nbody_system, param, dth)
                  call tp%gr_pos_kick(nbody_system, param, dth)
               end if

               call pl%drift(nbody_system, param, dtl)
               call tp%drift(nbody_system, param, dtl)

               if (lencounter) call nbody_system%recursive_step(param, t+(j-1)*dtl, irecp)
               nbody_system%irec = ireci

               if (param%lgr) then
                  call pl%gr_pos_kick(nbody_system, param, dth)
                  call tp%gr_pos_kick(nbody_system, param, dth)
               end if

               call plpl_encounter%kick(nbody_system, dth, irecp, 1)
               call pltp_encounter%kick(nbody_system, dth, irecp, 1)
               if (ireci /= 0) then
                  call plpl_encounter%kick(nbody_system, dth, irecp, -1)
                  call pltp_encounter%kick(nbody_system, dth, irecp, -1)
               end if

               if (param%lclose) then
                  call plpl_encounter%collision_check(nbody_system, param, t+j*dtl, dtl, ireci, lplpl_collision) 
                  call pltp_encounter%collision_check(nbody_system, param, t+j*dtl, dtl, ireci, lpltp_collision) 

                  if (lplpl_collision) call plpl_encounter%resolve_collision(nbody_system, param, t+j*dtl, dtl, ireci)
                  if (lpltp_collision) call pltp_encounter%resolve_collision(nbody_system, param, t+j*dtl, dtl, ireci)
               end if
               if (param%lenc_save_trajectory) call self%encounter_history%take_snapshot(param, self, t+j*dtl, "trajectory") 

               call self%set_recur_levels(ireci)

            end do
         end associate
      end select
      end select
      end select
      end select
      end select

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
      class(swiftest_parameters),   intent(in)    :: param !! Current run configuration parameters with SyMBA additions
      ! Internals
      integer(I4B) :: i
      integer(I8B) :: nenc_old

      associate(nbody_system => self)
      select type(pl => nbody_system%pl)
      class is (symba_pl)
      select type(tp => nbody_system%tp)
      class is (symba_tp)
         associate(npl => pl%nbody, ntp => tp%nbody)
            nenc_old = nbody_system%plpl_encounter%nenc
            call nbody_system%plpl_encounter%setup(0_I8B)
            call nbody_system%plpl_collision%setup(0_I8B)
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
               call nbody_system%plpl_encounter%setup(nenc_old) ! This resizes the pl-pl encounter list to be the same size as it was the last step, to decrease the number of potential resize operations that have to be one inside the step
               nbody_system%plpl_encounter%nenc = 0 ! Sets the true number of encounters back to 0 after resizing
               nbody_system%plpl_encounter%lcollision = .false.
            end if
      
            nenc_old = nbody_system%pltp_encounter%nenc
            call nbody_system%pltp_encounter%setup(0_I8B)
            call nbody_system%pltp_collision%setup(0_I8B)
            if (ntp > 0) then
               tp%lcollision(1:ntp) = .false.
               tp%nplenc(1:ntp) = 0 
               tp%levelg(1:ntp) = -1
               tp%levelm(1:ntp) = -1
               tp%lmask(1:ntp) = .true.
               tp%ldiscard(1:ntp) = .false.
               call nbody_system%pltp_encounter%setup(nenc_old)! This resizes the pl-tp encounter list to be the same size as it was the last step, to decrease the number of potential resize operations that have to be one inside the step
               nbody_system%pltp_encounter%nenc = 0 ! Sets the true number of encounters back to 0 after resizing
               nbody_system%pltp_encounter%lcollision = .false.
            end if

            call nbody_system%pl_adds%setup(0, param)
            call nbody_system%pl_discards%setup(0, param)

            tp%lfirst = param%lfirstkick
            pl%lfirst = param%lfirstkick

         end associate
      end select
      end select
      end associate

      return
   end subroutine symba_step_reset_system

end submodule s_symba_step
