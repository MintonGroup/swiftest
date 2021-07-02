submodule(rmvs_classes) s_rmvs_getacch
   use swiftest
contains  
   module subroutine rmvs_getacch_tp(self, cb, pl, param, t, xh)
      !! author: David A. Minton
      !!
      !! Compute the oblateness acceleration in the inner encounter region with planets 
      !! 
      !! Performs a similar task as David E. Kaufmann's Swifter routine rmvs_getacch_tp.f90, but 
      !! uses object polymorphism, and so is not directly adapted.
      implicit none
      ! Arguments
      class(rmvs_tp),                intent(inout) :: self   !! RMVS test particle data structure
      class(swiftest_cb),            intent(inout) :: cb     !! Swiftest central body particle data structuree 
      class(whm_pl),                 intent(inout) :: pl     !! WHM massive body particle data structure. 
      class(swiftest_parameters), intent(in)    :: param !! Current run configuration parameters of  parameter
      real(DP),                      intent(in)    :: t      !! Current time
      real(DP), dimension(:,:),      intent(in)    :: xh     !! Heliocentric positions of planets
      ! Internals
      type(swiftest_parameters)                 :: param_planetocen
      real(DP), dimension(:, :), allocatable       :: xh_original
      integer(I4B)                                 :: i


      associate(tp => self, ntp => self%nbody, ipleP => self%ipleP, &
                inner_index => self%index, cb_heliocentric => self%cb_heliocentric)
         
         if (tp%lplanetocentric) then  ! This is a close encounter step, so any accelerations requiring heliocentric position values
                                       ! must be handeled outside the normal WHM method call
            select type(pl)
            class is (rmvs_pl)
            select type (cb)
            class is (rmvs_cb)
               associate(xpc => pl%xh, xpct => self%xh)
                  allocate(xh_original, source=tp%xh)
                  param_planetocen = param
                  ! Temporarily turn off the heliocentric-dependent acceleration terms during an inner encounter
                  param_planetocen%loblatecb = .false.
                  param_planetocen%lextra_force = .false.
                  param_planetocen%lgr = .false.
                  ! Now compute the planetocentric values of acceleration
                  call whm_getacch_tp(tp, cb, pl, param_planetocen, t, xh)

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

                  if (param%lextra_force) call tp%user_getacch(cb_heliocentric, param, t)
                  if (param%lgr) call tp%gr_getacch(cb_heliocentric, param)
                  
                  tp%xh(:,:) = xh_original(:,:)
               end associate
            end select
            end select
         else ! Not a close encounter, so just proceded with the standard WHM method
            call whm_getacch_tp(tp, cb, pl, param, t, xh)
         end if

      end associate

      return

   end subroutine rmvs_getacch_tp

end submodule s_rmvs_getacch