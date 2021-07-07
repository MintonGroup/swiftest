submodule(rmvs_classes) s_rmvs_getacch
   use swiftest
contains  
   module subroutine rmvs_getacch_tp(self, system, param, t, xhp)

      !! author: David A. Minton
      !!
      !! Compute the oblateness acceleration in the inner encounter region with planets 
      !! 
      !! Performs a similar task as David E. Kaufmann's Swifter routine rmvs_getacch_tp.f90, but 
      !! uses object polymorphism, and so is not directly adapted.
      implicit none
      ! Arguments
      class(rmvs_tp),               intent(inout) :: self   !! RMVS test particle data structure
      class(swiftest_nbody_system), intent(inout) :: system !! Swiftest central body particle data structuree 
      class(swiftest_parameters),   intent(in)    :: param  !! Current run configuration parameters
      real(DP),                     intent(in)    :: t      !! Current time
      real(DP), dimension(:,:),     intent(in)    :: xhp    !! Heliocentric positions of planets at current substep
      ! Internals
      type(swiftest_parameters)                 :: param_planetocen
      real(DP), dimension(:, :), allocatable       :: xh_original
      integer(I4B)                                 :: i


      associate(tp => self, ntp => self%nbody, ipleP => self%ipleP, inner_index => self%index, cb_heliocentric => self%cb_heliocentric)
         select type(system)
         class is (rmvs_nbody_system)
            if (system%lplanetocentric) then  ! This is a close encounter step, so any accelerations requiring heliocentric position values
                                          ! must be handeled outside the normal WHM method call
               select type(pl => system%pl)
               class is (rmvs_pl)
               select type (cb => system%cb)
               class is (rmvs_cb)
                  associate(xpc => pl%xh, xpct => self%xh, apct => self%ah)
                     allocate(xh_original, source=tp%xh)
                     param_planetocen = param
                     ! Temporarily turn off the heliocentric-dependent acceleration terms during an inner encounter
                     param_planetocen%loblatecb = .false.
                     param_planetocen%lextra_force = .false.
                     param_planetocen%lgr = .false.
                     ! Now compute the planetocentric values of acceleration
                     call whm_getacch_tp(tp, system, param_planetocen, t, xhp)

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
                     ! Swap the planetocentric and heliocentric position vectors
                     tp%xh(:,:) = tp%xheliocentric(:,:)
                     if (param%loblatecb) then
                        ! Put in the current encountering planet's oblateness acceleration as the central body's
                        if (tp%lfirst) then
                           cb_heliocentric%aobl(:) = cb%inner(inner_index - 1)%aobl(:,1)
                        else
                           cb_heliocentric%aobl(:) = cb%inner(inner_index    )%aobl(:,1)
                        end if
                        call tp%obl_acc(cb_heliocentric)
                     end if

                     if (param%lextra_force) call tp%user_getacch(system, param, t)
                     if (param%lgr) call tp%gr_get_accel(param)
                     
                     tp%xh(:,:) = xh_original(:,:)
                  end associate
               end select
               end select
            else ! Not a close encounter, so just proceded with the standard WHM method
               call whm_getacch_tp(tp, system, param, t, xhp)
            end if
         end select

      end associate

      return

   end subroutine rmvs_getacch_tp

end submodule s_rmvs_getacch