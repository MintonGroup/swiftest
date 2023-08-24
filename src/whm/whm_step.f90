!! Copyright 2022 - David Minton, Carlisle Wishard, Jennifer Pouplin, Jake Elliott, & Dana Singh
!! This file is part of Swiftest.
!! Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
!! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!! Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
!! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!! You should have received a copy of the GNU General Public License along with Swiftest. 
!! If not, see: https://www.gnu.org/licenses. 

submodule(whm) s_whm_step
   use swiftest
contains

   module subroutine whm_step_system(self, param, t, dt)
      !! author: David A. Minton
      !!
      !! Step massive bodies and and active test particles ahead in heliocentric coordinates
      !! 
      !! Adapted from Hal Levison's Swift routine step_kdk.f
      !! Adapted from David E. Kaufmann's Swifter routine whm_step.f90
      implicit none
      ! Arguments
      class(whm_nbody_system),    intent(inout) :: self  !! WHM nbody system object
      class(swiftest_parameters), intent(inout) :: param !! Current run configuration parameters 
      real(DP),                   intent(in)    :: t     !! Current simulation time
      real(DP),                   intent(in)    :: dt    !! Current stepsize

      associate(nbody_system => self, cb => self%cb, pl => self%pl, tp => self%tp)
         tp%lfirst = pl%lfirst
         call pl%step(nbody_system, param, t, dt)
         call tp%step(nbody_system, param, t, dt)
         ! if (param%ltides) call nbody_system%step_spin(param, t, dt)
      end associate
      return
   end subroutine whm_step_system 


   module subroutine whm_step_pl(self, nbody_system, param, t, dt)
      !! author: David A. Minton
      !!
      !! Step planets ahead using kick-drift-kick algorithm
      !! 
      !! Adapted from Hal Levison's Swift routine step_kdk_pl.f
      !! Adapted from David E. Kaufmann's Swifter routine whm_step_pl.f90
      implicit none
      ! Arguments
      class(whm_pl),                intent(inout) :: self   !! WHM massive body particle data structure
      class(swiftest_nbody_system), intent(inout) :: nbody_system !! Swiftest nbody_system object
      class(swiftest_parameters),   intent(inout) :: param  !! Current run configuration parameters 
      real(DP),                     intent(in)    :: t     !! Current simulation time
      real(DP),                     intent(in)    :: dt    !! Current stepsize
      ! Internals
      real(DP)                                    :: dth
      
      if (self%nbody == 0) return

      associate(pl => self, cb => nbody_system%cb)
         dth = 0.5_DP * dt
         call pl%kick(nbody_system, param, t, dth,lbeg=.true.)
         call pl%vh2vj(cb) 
         if (param%lgr) call pl%gr_pos_kick(nbody_system, param, dth)
         call pl%drift(nbody_system, param, dt)
         if (param%lgr) call pl%gr_pos_kick(nbody_system, param, dth)
         call pl%j2h(cb)
         call pl%kick(nbody_system, param, t + dt, dth, lbeg=.false.)
      end associate

      return
   end subroutine whm_step_pl


   module subroutine whm_step_tp(self, nbody_system, param, t, dt)
      !! author: David A. Minton
      !!
      !! Step active test particles ahead using kick-drift-kick algorithm
      !! 
      !! Adapted from Hal Levison's Swift routine step_kdk_tp.f
      !! Adapted from David E. Kaufmann's Swifter routine whm_step_tp.f90
      implicit none
      ! Arguments
      class(whm_tp),                intent(inout) :: self   !! WHM test particle data structure
      class(swiftest_nbody_system), intent(inout) :: nbody_system !! Swiftest nbody_system object
      class(swiftest_parameters),   intent(inout) :: param  !! Current run configuration parameters 
      real(DP),                     intent(in)    :: t     !! Current simulation time
      real(DP),                     intent(in)    :: dt    !! Current stepsize
      ! Internals
      real(DP)                                    :: dth

      if (self%nbody == 0) return

      select type(nbody_system)
      class is (whm_nbody_system)
         associate(tp => self, cb => nbody_system%cb, pl => nbody_system%pl)
            dth = 0.5_DP * dt
            call tp%kick(nbody_system, param, t, dth, lbeg=.true.)
            if (param%lgr) call tp%gr_pos_kick(nbody_system, param, dth)
            call tp%drift(nbody_system, param, dt)
            if (param%lgr) call tp%gr_pos_kick(nbody_system, param, dth)
            call tp%kick(nbody_system, param, t + dt, dth, lbeg=.false.)
         end associate
      end select

      return
   end subroutine whm_step_tp   

end submodule s_whm_step
