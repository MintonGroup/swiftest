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
      associate(pl => self, npl => self%nbody)
         select type(system)
         class is (symba_nbody_system)
            pl%lmask(1:npl) = pl%status(1:npl) /= INACTIVE .and. pl%levelg(1:npl) == system%irec
            call helio_drift_body(pl, system, param, dt)
            pl%lmask(1:npl) = pl%status(1:npl) /= INACTIVE 
         end select
      end associate

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
      associate (tp => self, ntp => self%nbody)
         select type(system)
         class is (symba_nbody_system)
            tp%lmask(1:ntp) = tp%status(1:ntp) /= INACTIVE .and. tp%levelg(1:ntp) == system%irec
            call helio_drift_body(tp, system, param, dt)
            tp%lmask(1:ntp) = tp%status(1:ntp) /= INACTIVE 
         end select
      end associate

      return
   end subroutine symba_drift_tp

end submodule s_symba_drift
