submodule (symba_classes) s_symba_collision
   use swiftest
contains
   module subroutine symba_collision_check_plplenc(self, system, param, t, dt, irec)
      !! author: Jennifer L.L. Pouplin, Carlisle A. wishard, and David A. Minton
      !!
      !! Check for merger between massive bodies in SyMBA. If the user has turned on the FRAGMENTATION feature, it will call the 
      !! symba_regime subroutine to determine what kind of collision will occur.
      !! 
      !! Adapted from David E. Kaufmann's Swifter routine symba_merge_pl.f90
      !!
      !! Adapted from Hal Levison's Swift routine symba5_merge.f
      implicit none
      class(symba_plplenc),       intent(inout) :: self   !! SyMBA pl-tp encounter list object
      class(symba_nbody_system),  intent(inout) :: system !! SyMBA nbody system object
      class(swiftest_parameters), intent(in)    :: param  !! Current run configuration parameters 
      real(DP),                   intent(in)    :: t      !! current time
      real(DP),                   intent(in)    :: dt     !! step size
      integer(I4B),               intent(in)    :: irec   !! Current recursion level
   end subroutine symba_collision_check_plplenc


   module subroutine symba_collision_check_pltpenc(self, system, param, t, dt, irec)
      !! author: David A. Minton
      !!
      !! Check for merger between massive bodies and test particles in SyMBA
      !! 
      !! Adapted from David E. Kaufmann's Swifter routine symba_merge_tp.f90
      !!
      !! Adapted from Hal Levison's Swift routine symba5_merge.f
      implicit none
      class(symba_pltpenc),       intent(inout) :: self   !! SyMBA pl-tp encounter list object
      class(symba_nbody_system),  intent(inout) :: system !! SyMBA nbody system object
      class(swiftest_parameters), intent(in)    :: param  !! Current run configuration parameters 
      real(DP),                   intent(in)    :: t      !! current time
      real(DP),                   intent(in)    :: dt     !! step size
      integer(I4B),               intent(in)    :: irec   !! Current recursion level
   end subroutine symba_collision_check_pltpenc

end submodule s_symba_collision