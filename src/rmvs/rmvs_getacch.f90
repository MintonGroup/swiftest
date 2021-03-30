submodule(rmvs_classes) s_rmvs_getacch
contains  

   module subroutine rmvs_getacch_in_tp(self, cb, pl, config, t, xh)
      !! author: David A. Minton
      !!
      !! Compute the oblateness acceleration in the inner encounter region with planets 
      !! 
      !! Performs a similar task as David E. Kaufmann's Swifter routine rmvs_getacch_tp.f90, but 
      !! uses object polymorphism, and so is not directly adapted.
      use swiftest
      implicit none
      !! Arguments
      class(rmvs_tp),                intent(inout) :: self   !! RMVS test particle data structure
      class(swiftest_cb),            intent(inout) :: cb     !! Swiftest central body particle data structuree 
      class(whm_pl),                 intent(inout) :: pl     !! WHM massive body particle data structure. 
      class(swiftest_configuration), intent(in)    :: config !! Input collection of  parameter
      real(DP),                      intent(in)    :: t      !! Current time
      real(DP), dimension(:,:),      intent(in)    :: xh     !! Heliocentric positions of planets
      !! Internals
      type(swiftest_configuration)                 :: config_planetocen
      real(DP), dimension(:, :), allocatable       :: xh_original
      integer(I4B)                                 :: i


      associate(tp => self, cb_heliocen => self%cb, ipleP => self%ipleP, index => self%index, &
         xht => self%xh, vht => self%vh, aht => self%ah)
         if (tp%lplanetocentric) then  ! This is a close encounter step, so any accelerations requiring heliocentric position values
                                       ! must be handeled outside the normal WHM method call
            allocate(xh_original, source=tp%xh)
            config_planetocen = config
            ! Temporarily turn off the heliocentric-dependent acceleration terms during an inner encounter
            config_planetocen%loblatecb = .false.
            config_planetocen%lextra_force = .false.
            config_planetocen%lgr = .false.
            call whm_getacch_tp(tp, cb, pl, config_planetocen, t, xh)

            ! Now compute the heliocentric values of acceleration
            if (tp%lfirst) then
               do i = 1, NDIM
                  tp%xheliocen(i,:) = tp%xh(i,:) + tp%xh_pl(i, index - 1)
               end do
            else
               do i = 1, NDIM
                  tp%xheliocen(i,:) = tp%xh(i,:) + tp%xh_pl(i, index)
               end do
            end if

            tp%xh(:,:) = tp%xheliocen(:,:)
            if (config%loblatecb) then
               ! Put in the current encountering planet's oblateness acceleration as the central body's
               if (tp%lfirst) then
                  cb_heliocen%aobl(:) = tp%aoblin_pl(:, index - 1)
               else
                  cb_heliocen%aobl(:) = tp%aoblin_pl(:, index)
               end if
               call tp%obl_acc(cb_heliocen)
            end if

            if (config%lextra_force) call tp%user_getacch(cb_heliocen, config, t)
            if (config%lgr) call tp%gr_getacch(cb_heliocen, config)
            tp%xh(:,:) = xh_original(:,:)
         else ! Not a close encounter, so just proceded with the standard WHM method
            call whm_getacch_tp(tp, cb, pl, config, t, xh)
         end if

      end associate

      return

   end subroutine rmvs_getacch_in_tp

end submodule s_rmvs_getacch