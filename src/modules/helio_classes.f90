module helio_classes
   !! author: The Purdue Swiftest Team - David A. Minton, Carlisle A. Wishard, Jennifer L.L. Pouplin, and Jacob R. Elliott
   !!
   !! Definition of classes and methods specific to the Democratic Heliocentric Method
   !! Adapted from David E. Kaufmann's Swifter routine: helio.f90
   use swiftest_globals
   use rmvs_classes, only : rmvs_cb, rmvs_pl, rmvs_tp, rmvs_nbody_system
   implicit none


   !********************************************************************************************************************************
   ! helio_cb class definitions and method interfaces
   !*******************************************************************************************************************************
   !> RMVS central body particle class
   type, public, extends(rmvs_cb) :: helio_cb
   contains
   end type helio_cb

   !********************************************************************************************************************************
   !                                    helio_tp class definitions and method interfaces
   !*******************************************************************************************************************************

   !! Helio test particle class
   type, public, extends(rmvs_tp) :: helio_tp
   contains
   end type helio_tp

   interface
   end interface

   !********************************************************************************************************************************
   !                                    helio_pl class definitions and method interfaces
   !*******************************************************************************************************************************

   !! Helio massive body particle class
   type, public, extends(rmvs_pl) :: helio_pl
   contains
   end type helio_pl

   !********************************************************************************************************************************
   !  rmvs_nbody_system class definitions and method interfaces
   !********************************************************************************************************************************
   type, public, extends(rmvs_nbody_system) :: helio_nbody_system
      !> In the RMVS integrator, only test particles are discarded
   contains
      private
      !> Replace the abstract procedures with concrete ones
      procedure, public :: initialize    => helio_setup_system  !! Performs RMVS-specific initilization steps, like calculating the Jacobi masses
   end type helio_nbody_system

   interface
   module subroutine helio_setup_system(self, config)
      use swiftest_classes
      implicit none
      class(helio_nbody_system),      intent(inout) :: self    !! RMVS system object
      class(swiftest_configuration), intent(inout) :: config  !! Input collection of  configuration parameters 
   end subroutine helio_setup_system
   end interface



end module helio_classes
