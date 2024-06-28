! Copyight 2022 - David Minton, Carlisle Wishard, Jennifer Pouplin, Jake Elliott, & Dana Singh
! This file is part of Swiftest.
! Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
! Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
! You should have received a copy of the GNU General Public License along with Swiftest. 
! If not, see: https://www.gnu.org/licenses. 

submodule(rmvs) s_rmvs_kick
   use swiftest
contains  

   module subroutine rmvs_kick_getacch_tp(self, nbody_system, param, t, lbeg)
      !! author: David A. Minton
      !!
      !! Compute the oblateness acceleration in the inner encounter region with planets 
      !! 
      !! Performs a similar task as David E. Kaufmann's Swifter routine rmvs_kick_getacch_tp.f90, but 
      !! uses object polymorphism, and so is not directly adapted.
      implicit none
      ! Arguments
      class(rmvs_tp),               intent(inout) :: self   !! RMVS test particle data structure
      class(swiftest_nbody_system), intent(inout) :: nbody_system !! Swiftest central body particle data structuree 
      class(swiftest_parameters),   intent(inout) :: param  !! Current run configuration parameters
      real(DP),                     intent(in)    :: t      !! Current time
      logical,                      intent(in)    :: lbeg   !! Logical flag that determines whether or not this is the beginning or end of the step
      ! Internals
      class(swiftest_parameters), allocatable   :: param_planetocen
      real(DP), dimension(:, :), allocatable    :: rh_original
      real(DP)                                  :: GMcb_original
      integer(I4B)                              :: i, ntp, inner_index

      if (self%nbody == 0) return

      associate(tp => self, ipleP => self%ipleP)
         ntp = self%nbody
         inner_index = self%index
         select type(nbody_system)
         class is (rmvs_nbody_system)
            if (nbody_system%lplanetocentric) then  ! This is a close encounter step, so any accelerations requiring heliocentric position values
                                              ! must be handeled outside the normal WHM method call
               select type(pl => nbody_system%pl)
               class is (rmvs_pl)
                  select type (cb => nbody_system%cb)
                  class is (rmvs_cb)
                     associate(xpc => pl%rh, xpct => self%rh, apct => self%ah, system_planetocen => nbody_system)
                        system_planetocen%lbeg = lbeg

                        ! Save the original heliocentric position for later
                        allocate(rh_original, source=tp%rh)

                        ! Temporarily turn off the heliocentric-dependent acceleration terms during an inner encounter using a copy of the parameter list with all of the heliocentric-specific acceleration terms turned off
                        allocate(param_planetocen, source=param)
                        param_planetocen%lnon_spherical_cb = .false.
                        param_planetocen%lextra_force = .false.
                        param_planetocen%lgr = .false.

                        ! Compute the planetocentric values of acceleration
                        call whm_kick_getacch_tp(tp, system_planetocen, param_planetocen, t, lbeg)

                        ! Now compute any heliocentric values of acceleration 
                        if (tp%lfirst) then
#ifdef DOCONLOC
                           do concurrent(i = 1:ntp, tp%lmask(i)) shared(tp)
#else
                           do concurrent(i = 1:ntp, tp%lmask(i))
#endif
                              tp%rheliocentric(:,i) = tp%rh(:,i) + cb%inner(inner_index - 1)%x(:,1)
                           end do
                        else
#ifdef DOCONLOC
                           do concurrent(i = 1:ntp, tp%lmask(i)) shared(tp)
#else
                           do concurrent(i = 1:ntp, tp%lmask(i))
#endif
                              tp%rheliocentric(:,i) = tp%rh(:,i) + cb%inner(inner_index    )%x(:,1)
                           end do
                        end if

                        ! Swap the planetocentric and heliocentric position vectors and central body masses
#ifdef DOCONLOC
                        do concurrent(i = 1:ntp, tp%lmask(i)) shared(tp)
#else
                        do concurrent(i = 1:ntp, tp%lmask(i))
#endif
                           tp%rh(:, i) = tp%rheliocentric(:, i)
                        end do
                        GMcb_original = cb%Gmass
                        cb%Gmass = tp%cb_heliocentric%Gmass

                        ! If the heliocentric-specifc acceleration terms are requested, compute those now
                        if (param%lnon_spherical_cb) call tp%accel_non_spherical_cb(system_planetocen)
                        if (param%lextra_force) call tp%accel_user(system_planetocen, param, t, lbeg)
                        if (param%lgr) call tp%accel_gr(param)

                        ! Put everything back the way we found it
                        call move_alloc(rh_original, tp%rh)
                        cb%Gmass = GMcb_original

                     end associate
                  end select
               end select
            else ! Not a close encounter, so just proceded with the standard WHM method
               call whm_kick_getacch_tp(tp, nbody_system, param, t, lbeg)
            end if
         end select
      end associate

      return
   end subroutine rmvs_kick_getacch_tp

end submodule s_rmvs_kick