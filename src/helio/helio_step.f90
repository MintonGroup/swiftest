!! Copyright 2022 - David Minton, Carlisle Wishard, Jennifer Pouplin, Jake Elliott, & Dana Singh
!! This file is part of Swiftest.
!! Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
!! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!! Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
!! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!! You should have received a copy of the GNU General Public License along with Swiftest. 
!! If not, see: https://www.gnu.org/licenses. 

submodule(helio) s_helio_step
   use swiftest
contains

   module subroutine helio_step_system(self, param, t, dt)
      !! author: David A. Minton
      !!
      !! Step massive bodies and and active test particles ahead in heliocentric coordinates.
      !!
      !! Currently there's no difference between this and the WHM nbody_system stepper, so this is just
      !! a wrapper function to keep the method calls consistent for inherited types.
      !! 
      !! Adapted from Hal Levison's Swift routine step_kdk.f
      !! Adapted from David E. Kaufmann's Swifter routine helio_step.f90
      implicit none
      ! Arguments
      class(helio_nbody_system),  intent(inout) :: self   !! Helio nbody system object
      class(swiftest_parameters), intent(inout) :: param  !! Current run configuration parameters
      real(DP),                   intent(in)    :: t      !! Simulation time
      real(DP),                   intent(in)    :: dt     !! Current stepsize

      call whm_step_system(self, param, t, dt)

      return
   end subroutine helio_step_system 


   module subroutine helio_step_pl(self, nbody_system, param, t, dt)
      !! author: David A. Minton
      !!
      !! Step massive bodies ahead Democratic Heliocentric method
      !!
      !! Adapted from David E. Kaufmann's Swifter helio_step_pl.f90
      !! Adapted from Hal Levison's Swift routine helio_step_pl.f
      implicit none
      ! Arguments
      class(helio_pl),              intent(inout) :: self   !! Helio massive body particle data structure
      class(swiftest_nbody_system), intent(inout) :: nbody_system !! Swiftest nbody system
      class(swiftest_parameters),   intent(inout) :: param  !! Current run configuration parameters 
      real(DP),                     intent(in)    :: t      !! Current simulation time
      real(DP),                     intent(in)    :: dt     !! Stepsize
      ! Internals 
      real(DP) :: dth   !! Half step size 

      if (self%nbody == 0) return

      associate(pl => self)
         select type(cb => nbody_system%cb)
         class is (helio_cb)
            dth = 0.5_DP * dt
            if (pl%lfirst) then
               call pl%vh2vb(cb)
               pl%lfirst = .false.
            end if
            call pl%lindrift(cb, dth, lbeg=.true.)
            call pl%kick(nbody_system, param, t, dth, lbeg=.true.)
            if (param%lgr) call pl%gr_pos_kick(nbody_system, param, dth)
            call pl%drift(nbody_system, param, dt)
            if (param%lgr) call pl%gr_pos_kick(nbody_system, param, dth)
            call pl%kick(nbody_system, param, t + dt, dth, lbeg=.false.)
            call pl%lindrift(cb, dth, lbeg=.false.)
            call pl%vb2vh(cb)
         end select
      end associate
   
      return
   end subroutine helio_step_pl


   module subroutine helio_step_tp(self, nbody_system, param, t, dt)

      !! author: David A. Minton
      !!
      !! Step active test particles ahead using Democratic Heliocentric method
      !!
      !! Adapted from David E. Kaufmann's Swifter routine helio_step_tp.f90
      !! Adapted from Hal Levison's Swift routine helio_step_tp.f
      implicit none
      ! Arguments
      class(helio_tp),              intent(inout) :: self    !! Helio test particle data structure
      class(swiftest_nbody_system), intent(inout) :: nbody_system !! Swiftest nbody system
      class(swiftest_parameters),   intent(inout) :: param  !! Current run configuration parameters 
      real(DP),                     intent(in)    :: t      !! Current simulation time
      real(DP),                     intent(in)    :: dt     !! Stepsize
      ! Internals
      real(DP) :: dth   !! Half step size 
   
      if (self%nbody == 0) return

      associate(tp => self)
         select type(cb => nbody_system%cb)
         class is (helio_cb)
            dth = 0.5_DP * dt
            if (tp%lfirst) then
               call tp%vh2vb(vbcb = -cb%ptbeg)
               tp%lfirst = .false.
            end if
            call tp%lindrift(cb, dth, lbeg=.true.)
            call tp%kick(nbody_system, param, t, dth, lbeg=.true.)
            if (param%lgr) call tp%gr_pos_kick(nbody_system, param, dth)
            call tp%drift(nbody_system, param, dt)
            if (param%lgr) call tp%gr_pos_kick(nbody_system, param, dth)
            call tp%kick(nbody_system, param, t + dt, dth, lbeg=.false.)
            call tp%lindrift(cb, dth, lbeg=.false.)
            call tp%vb2vh(vbcb = -cb%ptend)
         end select
      end associate
   
      return
   end subroutine helio_step_tp

end submodule s_helio_step
