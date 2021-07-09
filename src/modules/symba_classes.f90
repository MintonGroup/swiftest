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
      private
   end type symba_cb

   !********************************************************************************************************************************
   !                                    symba_pl class definitions and method interfaces
   !*******************************************************************************************************************************

   !! Helio massive body particle class
   type, public, extends(helio_pl) :: symba_pl
   contains
      private
      procedure, public :: discard         => symba_discard_pl         !! Process massive body discards
      procedure, public :: encounter_check => symba_encounter_check_pl !! Checks if massive bodies are going through close encounters with each other
   end type symba_pl

   !********************************************************************************************************************************
   !                                    symba_tp class definitions and method interfaces
   !*******************************************************************************************************************************

   !! Helio test particle class
   type, public, extends(helio_tp) :: symba_tp
   contains
      private
      procedure, public :: discard         => symba_discard_tp         !! process test particle discards
      procedure, public :: encounter_check => symba_encounter_check_tp !! Checks if any test particles are undergoing a close encounter with a massive body
   end type symba_tp

   interface
      module subroutine symba_discard_pl(self, system, param)
         use swiftest_classes, only : swiftest_nbody_system, swiftest_parameters
         implicit none
         class(symba_pl),              intent(inout) :: self   !! RMVS test particle object
         class(swiftest_nbody_system), intent(inout) :: system !! Swiftest nbody system object
         class(swiftest_parameters),   intent(in)    :: param  !! Current run configuration parameters 
      end subroutine symba_discard_pl

      module subroutine symba_discard_tp(self, system, param)
         use swiftest_classes, only : swiftest_nbody_system, swiftest_parameters
         implicit none
         class(symba_tp),              intent(inout) :: self   !! RMVS test particle object
         class(swiftest_nbody_system), intent(inout) :: system !! Swiftest nbody system object
         class(swiftest_parameters),   intent(in)    :: param  !! Current run configuration parameters 
      end subroutine symba_discard_tp

      module function symba_encounter_check_pl(self, system, dt) result(lencounter)
         implicit none
         class(symba_pl),           intent(inout) :: self       !! RMVS test particle object  
         class(symba_nbody_system), intent(inout) :: system     !! RMVS nbody system object
         real(DP),                  intent(in)    :: dt         !! step size
         logical                                  :: lencounter !! Returns true if there is at least one close encounter      
      end function symba_encounter_check_pl

      module function symba_encounter_check_tp(self, system, dt) result(lencounter)
         implicit none
         class(symba_tp),           intent(inout) :: self       !! RMVS test particle object  
         class(symba_nbody_system), intent(inout) :: system     !! RMVS nbody system object
         real(DP),                  intent(in)    :: dt         !! step size
         logical                                  :: lencounter !! Returns true if there is at least one close encounter      
      end function symba_encounter_check_tp

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