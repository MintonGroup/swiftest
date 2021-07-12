submodule (symba_classes) s_symba_encounter_check
   use swiftest
contains
   module function symba_encounter_check_pl(self, system, dt) result(lencounter)
      implicit none
      ! Arguments
      class(symba_pl),           intent(inout) :: self       !! SyMBA test particle object  
      class(symba_nbody_system), intent(inout) :: system     !! SyMBA nbody system object
      real(DP),                  intent(in)    :: dt         !! step size
      ! Result
      logical                                  :: lencounter !! Returns true if there is at least one close encounter      

      lencounter = .false.
      return
   end function symba_encounter_check_pl

   module function symba_encounter_check_tp(self, system, dt) result(lencounter)
      implicit none
      ! Arguments
      class(symba_tp),           intent(inout) :: self       !! SyMBA test particle object  
      class(symba_nbody_system), intent(inout) :: system     !! SyMBA nbody system object
      real(DP),                  intent(in)    :: dt         !! step size
      ! Result
      logical                                  :: lencounter !! Returns true if there is at least one close encounter      

      lencounter = .false.
      return
   end function symba_encounter_check_tp

end submodule s_symba_encounter_check