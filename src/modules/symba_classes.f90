module symba_classes
   !! author: The Purdue Swiftest Team - David A. Minton, Carlisle A. Wishard, Jennifer L.L. Pouplin, and Jacob R. Elliott
   !!
   !! Definition of classes and methods specific to the Democratic Heliocentric Method
   !! Adapted from David E. Kaufmann's Swifter routine: helio.f90
   use swiftest_globals
   use helio_classes, only : helio_cb, helio_pl, helio_tp, helio_nbody_system
   implicit none


   !********************************************************************************************************************************
   !  symba_nbody_system class definitions and method interfaces
   !********************************************************************************************************************************
   type, public, extends(helio_nbody_system) :: symba_nbody_system
   contains
      private
      procedure, public :: step          => symba_step_system         !! Advance the SyMBA nbody system forward in time by one step
      procedure, public :: interp        => symba_step_interp_system  !! Perform an interpolation step on the SymBA nbody system 
   end type symba_nbody_system

   !********************************************************************************************************************************
   ! symba_cb class definitions and method interfaces
   !*******************************************************************************************************************************
   !> Helio central body particle class
   type, public, extends(helio_cb) :: symba_cb
   contains
   end type symba_cb

   !********************************************************************************************************************************
   !                                    symba_pl class definitions and method interfaces
   !*******************************************************************************************************************************

   !! Helio massive body particle class
   type, public, extends(helio_pl) :: symba_pl
   contains
   end type symba_pl

   !********************************************************************************************************************************
   !                                    symba_tp class definitions and method interfaces
   !*******************************************************************************************************************************

   !! Helio test particle class
   type, public, extends(helio_tp) :: symba_tp
   contains
   end type symba_tp

   interface
      module subroutine symba_step_system(self, param, t, dt)
         use swiftest_classes, only : swiftest_parameters
         implicit none
         class(symba_nbody_system),  intent(inout) :: self    !! RMVS nbody system object
         class(swiftest_parameters), intent(inout) :: param  !! Current run configuration parameters 
         real(DP),                   intent(in)    :: t      !! Simulation time
         real(DP),                   intent(in)    :: dt     !! Current stepsize
      end subroutine symba_step_system

      module subroutine symba_step_interp_system(self, param, t, dt)
         use swiftest_classes, only : swiftest_parameters
         implicit none
         class(symba_nbody_system),  intent(inout) :: self    !! RMVS nbody system object
         class(swiftest_parameters), intent(inout) :: param  !! Current run configuration parameters 
         real(DP),                   intent(in)    :: t      !! Simulation time
         real(DP),                   intent(in)    :: dt     !! Current stepsize
      end subroutine symba_step_interp_system
   end interface
end module symba_classes