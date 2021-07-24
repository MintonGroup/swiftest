submodule(symba_classes) s_symba_kick
   use swiftest
contains

   module subroutine symba_kick_plplenc(self, system, dt, irec, sgn)
      !! author: David A. Minton
      !!
      !! !! Kick barycentric velocities of massive bodies within SyMBA recursion.
      !!
      !! Adapted from David E. Kaufmann's Swifter routine: symba_kick.f90
      !! Adapted from Hal Levison's Swift routine symba5_kick.f
      implicit none
      class(symba_plplenc),      intent(in)    :: self   !! SyMBA pl-pl encounter list object
      class(symba_nbody_system), intent(inout) :: system !! SyMBA nbody system object
      real(DP),                  intent(in)    :: dt     !! step size
      integer(I4B),              intent(in)    :: irec   !! Current recursion level
      integer(I4B),              intent(in)    :: sgn    !! sign to be applied to acceleration
   end subroutine symba_kick_plplenc

   module subroutine symba_kick_pltpenc(self, system, dt, irec, sgn)
      !! author: David A. Minton
      !!
      !! !! Kick barycentric velocities of active test particles within SyMBA recursion.
      !!
      !! Adapted from David E. Kaufmann's Swifter routine: symba_kick.f90
      !! Adapted from Hal Levison's Swift routine symba5_kick.f
      implicit none
      class(symba_pltpenc),      intent(in)    :: self   !! SyMBA pl-tp encounter list object
      class(symba_nbody_system), intent(inout) :: system !! SyMBA nbody system object
      real(DP),                  intent(in)    :: dt     !! step size
      integer(I4B),              intent(in)    :: irec   !! Current recursion level
      integer(I4B),              intent(in)    :: sgn    !! sign to be applied to acceleration
   end subroutine symba_kick_pltpenc

end submodule s_symba_kick