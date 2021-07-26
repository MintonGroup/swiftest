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
      class(swiftest_parameters),   intent(in)    :: param  !! Current run configuration parameters
      real(DP),                     intent(in)    :: t      !! Current time
      logical,                      intent(in)    :: lbeg   !! Logical flag that determines whether or not this is the beginning or end of the step
      ! Internals
      type(swiftest_parameters)                 :: param_planetocen
      real(DP), dimension(:, :), allocatable    :: xh_original
      real(DP)                                  :: GMcb_original
      integer(I4B)                              :: i
      real(DP), dimension(:, :), allocatable    :: xhp

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

                        if (system_planetocen%lbeg) then
                           allocate(xhp, source=pl%xbeg)
                        else
                           allocate(xhp, source=pl%xend)
                        end if
               
                        allocate(xh_original, source=tp%xh)
                        param_planetocen = param
                        ! Temporarily turn off the heliocentric-dependent acceleration terms during an inner encounter
                        param_planetocen%loblatecb = .false.
                        param_planetocen%lextra_force = .false.
                        param_planetocen%lgr = .false.
                        ! Now compute the planetocentric values of acceleration
                        call whm_kick_getacch_tp(tp, system_planetocen, param_planetocen, t, lbeg)

                        ! Now compute any heliocentric values of acceleration 
                        if (tp%lfirst) then
                           do i = 1, ntp
                              tp%xheliocentric(:,i) = tp%xh(:,i) + cb%inner(inner_index - 1)%x(:,1)
                           end do
                        else
                           do i = 1, ntp
                              tp%xheliocentric(:,i) = tp%xh(:,i) + cb%inner(inner_index    )%x(:,1)
                           end do
                        end if
                        ! Swap the planetocentric and heliocentric position vectors and central body masses
                        tp%xh(:,:) = tp%xheliocentric(:,:)
                        GMcb_original = cb%Gmass
                        cb%Gmass = tp%cb_heliocentric%Gmass
                        if (param%loblatecb) call tp%accel_obl(system_planetocen)
                        if (param%lextra_force) call tp%accel_user(system_planetocen, param, t, lbeg)
                        if (param%lgr) call tp%accel_gr(param)
                        tp%xh(:,:) = xh_original(:,:)
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