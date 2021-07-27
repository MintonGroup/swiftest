submodule (symba_classes) s_symba_discard
   use swiftest
contains

   module subroutine symba_discard_pl(self, system, param)
      implicit none
      ! Arguments
      class(symba_pl),              intent(inout) :: self   !! SyMBA test particle object
      class(swiftest_nbody_system), intent(inout) :: system !! Swiftest nbody system object
      class(swiftest_parameters),   intent(in)    :: param  !! Current run configuration parameters 
      
      return
   end subroutine symba_discard_pl

end submodule s_symba_discard