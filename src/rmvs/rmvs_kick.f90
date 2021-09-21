submodule(rmvs_classes) s_rmvs_kick
   use swiftest
contains  

   module subroutine rmvs_kick_getacch_tp(self, system, param, t, lbeg)
      !! author: David A. Minton
      !!
      !! Compute the oblateness acceleration in the inner encounter region with planets 
      !! 
      !! Performs a similar task as David E. Kaufmann's Swifter routine rmvs_kick_getacch_tp.f90, but 
      !! uses object polymorphism, and so is not directly adapted.
      implicit none
      ! Arguments
      class(rmvs_tp),               intent(inout) :: self   !! RMVS test particle data structure
      class(swiftest_nbody_system), intent(inout) :: system !! Swiftest central body particle data structuree 
      class(swiftest_parameters),   intent(inout) :: param  !! Current run configuration parameters
      real(DP),                     intent(in)    :: t      !! Current time
      logical,                      intent(in)    :: lbeg   !! Logical flag that determines whether or not this is the beginning or end of the step
      ! Internals
      class(swiftest_parameters), allocatable   :: param_planetocen
      real(DP), dimension(:, :), allocatable    :: xh_original
      real(DP)                                  :: GMcb_original
      integer(I4B)                              :: i

      if (self%nbody == 0) return

      associate(tp => self, ntp => self%nbody, ipleP => self%ipleP, inner_index => self%index)
         select type(system)
         class is (rmvs_nbody_system)
            if (system%lplanetocentric) then  ! This is a close encounter step, so any accelerations requiring heliocentric position values
                                              ! must be handeled outside the normal WHM method call
               select type(pl => system%pl)
               class is (rmvs_pl)
                  select type (cb => system%cb)
                  class is (rmvs_cb)
                     associate(xpc => pl%xh, xpct => self%xh, apct => self%ah, system_planetocen => system)
                        system_planetocen%lbeg = lbeg

                        ! Save the original heliocentric position for later
                        allocate(xh_original, source=tp%xh)

                        ! Temporarily turn off the heliocentric-dependent acceleration terms during an inner encounter using a copy of the parameter list with all of the heliocentric-specific acceleration terms turned off
                        allocate(param_planetocen, source=param)
                        param_planetocen%loblatecb = .false.
                        param_planetocen%lextra_force = .false.
                        param_planetocen%lgr = .false.

                        ! Compute the planetocentric values of acceleration
                        call whm_kick_getacch_tp(tp, system_planetocen, param_planetocen, t, lbeg)

                        ! Now compute any heliocentric values of acceleration 
                        if (tp%lfirst) then
                           do concurrent(i = 1:ntp, tp%lmask(i))
                              tp%xheliocentric(:,i) = tp%xh(:,i) + cb%inner(inner_index - 1)%x(:,1)
                           end do
                        else
                           do concurrent(i = 1:ntp, tp%lmask(i))
                              tp%xheliocentric(:,i) = tp%xh(:,i) + cb%inner(inner_index    )%x(:,1)
                           end do
                        end if

                        ! Swap the planetocentric and heliocentric position vectors and central body masses
                        do concurrent(i = 1:ntp, tp%lmask(i))
                           tp%xh(:, i) = tp%xheliocentric(:, i)
                        end do
                        GMcb_original = cb%Gmass
                        cb%Gmass = tp%cb_heliocentric%Gmass

                        ! If the heliocentric-specifc acceleration terms are requested, compute those now
                        if (param%loblatecb) call tp%accel_obl(system_planetocen)
                        if (param%lextra_force) call tp%accel_user(system_planetocen, param, t, lbeg)
                        if (param%lgr) call tp%accel_gr(param)

                        ! Put everything back the way we found it
                        call move_alloc(xh_original, tp%xh)
                        cb%Gmass = GMcb_original

                     end associate
                  end select
               end select
            else ! Not a close encounter, so just proceded with the standard WHM method
               call whm_kick_getacch_tp(tp, system, param, t, lbeg)
            end if
         end select
      end associate

      return
   end subroutine rmvs_kick_getacch_tp

end submodule s_rmvs_kick