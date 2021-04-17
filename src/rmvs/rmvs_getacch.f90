submodule(rmvs_classes) s_rmvs_getacch
   use swiftest
contains  
   module subroutine rmvs_getacch_tp(self, cb, pl, config, t, xh)
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
      class(swiftest_configuration), intent(in)    :: config !! Input collection of  parameter
      real(DP),                      intent(in)    :: t      !! Current time
      real(DP), dimension(:,:),      intent(in)    :: xh     !! Heliocentric positions of planets
      ! Internals
      type(swiftest_configuration)                 :: config_planetocen
      real(DP), dimension(:, :), allocatable       :: xh_original
      integer(I4B)                                 :: i


      associate(tp => self, ntp => self%nbody, ipleP => self%ipleP, &
                enc_index => self%index, cb_heliocentric => self%cb_heliocentric)
         
         if (tp%lplanetocentric) then  ! This is a close encounter step, so any accelerations requiring heliocentric position values
                                       ! must be handeled outside the normal WHM method call
            select type(pl)
               class is (rmvs_pl)
                  allocate(xh_original, source=tp%xh)
                  config_planetocen = config
                  ! Temporarily turn off the heliocentric-dependent acceleration terms during an inner encounter
                  config_planetocen%loblatecb = .false.
                  config_planetocen%lextra_force = .false.
                  config_planetocen%lgr = .false.
                  ! Now compute the planetocentric values of acceleration
                  call whm_getacch_tp(tp, cb, pl%planetocentric(ipleP)%pl, config_planetocen, t, pl%xh)

                  ! Now compute any heliocentric values of acceleration 
                  if (tp%lfirst) then
                     do i = 1, ntp
                        tp%xheliocentric(:,i) = tp%xh(:,i) + pl%inner(enc_index - 1)%x(:,ipleP)
                     end do
                  else
                     do i = 1, ntp
                        tp%xheliocentric(:,i) = tp%xh(:,i) + pl%inner(enc_index    )%x(:,ipleP)
                     end do
                  end if
                  ! Swap the planetocentric and heliocentric position vectors
                  tp%xh(:,:) = tp%xheliocentric(:,:)
                  if (config%loblatecb) then
                     ! Put in the current encountering planet's oblateness acceleration as the central body's
                     if (tp%lfirst) then
                        cb_heliocentric%aobl(:) = pl%inner(enc_index - 1)%aobl(:, ipleP)
                     else
                        cb_heliocentric%aobl(:) = pl%inner(enc_index    )%aobl(:, ipleP)
                     end if
                     call tp%obl_acc(cb_heliocentric)
                  end if

                  if (config%lextra_force) call tp%user_getacch(cb_heliocentric, config, t)
                  if (config%lgr) call tp%gr_getacch(cb_heliocentric, config)
                  
                  call move_alloc(xh_original, tp%xh)
            end select
         else ! Not a close encounter, so just proceded with the standard WHM method
            call whm_getacch_tp(tp, cb, pl, config, t, pl%xh)
         end if

      end associate

      return

   end subroutine rmvs_getacch_tp

end submodule s_rmvs_getacch