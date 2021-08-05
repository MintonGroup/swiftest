  submodule (symba_classes) s_symba_drift
   use swiftest
contains

   module subroutine symba_drift_pl(self, system, param, dt)
      !! author: David A. Minton
      !!
      !! Wrapper function used to call the body drift routine from a symba_pl structure
      implicit none
      ! Arguments
      class(symba_pl),              intent(inout) :: self   !! Helio massive body object
      class(swiftest_nbody_system), intent(inout) :: system !! Swiftest nbody system object
      class(swiftest_parameters),   intent(in)    :: param  !! Current run configuration parameters 
      real(DP),                     intent(in)    :: dt     !! Stepsize

      if (self%nbody == 0) return

      select type(system)
      class is (symba_nbody_system)
         self%lmask(:) = self%status(:) /= INACTIVE .and. self%levelg(:) == system%irec
         call helio_drift_body(self, system, param, dt)
         self%lmask(:) = self%status(:) /= INACTIVE 
      end select

      return
   end subroutine symba_drift_pl


   module subroutine symba_drift_tp(self, system, param, dt)
      !! author: David A. Minton
      !!
      !! Wrapper function used to call the body drift routine from a symba_pl structure
      implicit none
      ! Arguments
      class(symba_tp),              intent(inout) :: self   !! Helio massive body object
      class(swiftest_nbody_system), intent(inout) :: system !! Swiftest nbody system object
      class(swiftest_parameters),   intent(in)    :: param  !! Current run configuration parameters 
      real(DP),                     intent(in)    :: dt     !! Stepsize

      if (self%nbody == 0) return

      select type(system)
      class is (symba_nbody_system)
         self%lmask(:) = self%status(:) /= INACTIVE .and. self%levelg(:) == system%irec
         call helio_drift_body(self, system, param, dt)
         self%lmask(:) = self%status(:) /= INACTIVE 
      end select

      return
   end subroutine symba_drift_tp

end submodule s_symba_drift
