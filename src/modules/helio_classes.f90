module helio_classes
   !! author: The Purdue Swiftest Team - David A. Minton, Carlisle A. Wishard, Jennifer L.L. Pouplin, and Jacob R. Elliott
   !!
   !! Definition of classes and methods specific to the Democratic Heliocentric Method
   !! Adapted from David E. Kaufmann's Swifter routine: helio.f90
   use swiftest_globals
   use swiftest_classes, only : swiftest_cb, swiftest_pl, swiftest_tp, swiftest_nbody_system
   use whm_classes, only : whm_nbody_system
   implicit none


   !********************************************************************************************************************************
   !  helio_nbody_system class definitions and method interfaces
   !********************************************************************************************************************************
   type, public, extends(whm_nbody_system) :: helio_nbody_system
   contains
   end type helio_nbody_system

   !********************************************************************************************************************************
   ! helio_cb class definitions and method interfaces
   !*******************************************************************************************************************************
   !> Helio central body particle class
   type, public, extends(swiftest_cb) :: helio_cb
      real(DP), dimension(NDIM)  :: ptbeg !! negative barycentric velocity of the central body at the beginning of time step
      real(DP), dimension(NDIM)  :: ptend !! negative barycentric velocity of the central body at the end of time step
   contains
   end type helio_cb

   !********************************************************************************************************************************
   !                                    helio_pl class definitions and method interfaces
   !*******************************************************************************************************************************

   !! Helio massive body particle class
   type, public, extends(swiftest_pl) :: helio_pl
   contains
      procedure, public :: vh2vb    => helio_coord_vh2vb_pl  !! Convert massive bodies from heliocentric to barycentric coordinates (velocity only)
      procedure, public :: vb2vh    => helio_coord_vb2vh_pl  !! Convert massive bodies from barycentric to heliocentric coordinates (velocity only)
      procedure, public :: drift    => helio_drift_pl        !! Method for Danby drift in Democratic Heliocentric coordinates 
      procedure, public :: lindrift => helio_drift_linear_pl !! Method for linear drift of massive bodies due to barycentric momentum of Sun
      procedure, public :: accel    => helio_getacch_pl      !! Compute heliocentric accelerations of massive bodies
      procedure, public :: kick     => helio_kickvb_pl       !! Kicks the barycentric velocities
      procedure, public :: step     => helio_step_pl         !! Steps the body forward one stepsize
   end type helio_pl

   !********************************************************************************************************************************
   !                                    helio_tp class definitions and method interfaces
   !*******************************************************************************************************************************

   !! Helio test particle class
   type, public, extends(swiftest_tp) :: helio_tp
   contains
      procedure, public :: vh2vb    => helio_coord_vh2vb_tp  !! Convert test particles from heliocentric to barycentric coordinates (velocity only)
      procedure, public :: vb2vh    => helio_coord_vb2vh_tp  !! Convert test particles from barycentric to heliocentric coordinates (velocity only)
      procedure, public :: drift    => helio_drift_tp        !! Method for Danby drift in Democratic Heliocentric coordinates 
      procedure, public :: lindrift => helio_drift_linear_tp !! Method for linear drift of massive bodies due to barycentric momentum of Sun
      procedure, public :: accel    => helio_getacch_tp      !! Compute heliocentric accelerations of massive bodies
      procedure, public :: kick     => helio_kickvb_tp       !! Kicks the barycentric velocities
      procedure, public :: step     => helio_step_tp         !! Steps the body forward one stepsize
   end type helio_tp

   interface
      module subroutine helio_coord_vb2vh_pl(self, cb)
         use swiftest_classes, only : swiftest_cb
         implicit none
         class(helio_pl),              intent(inout) :: self !! Helio massive body object
         class(swiftest_cb),           intent(inout) :: cb   !! Swiftest central body object
      end subroutine helio_coord_vb2vh_pl
   
      module subroutine helio_coord_vb2vh_tp(self, vbcb)
         implicit none
         class(helio_tp),              intent(inout) :: self !! Helio massive body object
         real(DP), dimension(:),       intent(in)    :: vbcb  !! Barycentric velocity of the central body
      end subroutine helio_coord_vb2vh_tp
   
      module subroutine helio_coord_vh2vb_pl(self, cb)
         use swiftest_classes, only : swiftest_cb
         implicit none
         class(helio_pl),        intent(inout) :: self !! Helio massive body object
         class(swiftest_cb),     intent(inout) :: cb   !! Swiftest central body object
      end subroutine helio_coord_vh2vb_pl
   
      module subroutine helio_coord_vh2vb_tp(self, vbcb)
         implicit none
         class(helio_tp),        intent(inout) :: self !! Helio massive body object
         real(DP), dimension(:), intent(in)    :: vbcb !! Barycentric velocity of the central body
      end subroutine helio_coord_vh2vb_tp
   
      module subroutine helio_drift_pl(self, system, param, dt)
         use swiftest_classes, only : swiftest_nbody_system, swiftest_parameters
         implicit none
         class(helio_pl),              intent(inout) :: self   !! Helio massive body object
         class(swiftest_nbody_system), intent(inout) :: system !! Swiftest nbody system object
         class(swiftest_parameters),   intent(in)    :: param  !! Current run configuration parameters of 
         real(DP),                     intent(in)    :: dt     !! Stepsize
      end subroutine helio_drift_pl

      module subroutine helio_drift_tp(self, system, param, dt)
         use swiftest_classes, only : swiftest_nbody_system, swiftest_parameters
         implicit none
         class(helio_tp),              intent(inout) :: self   !! Helio test particle object
         class(swiftest_nbody_system), intent(inout) :: system !! Swiftest nbody system object
         class(swiftest_parameters),   intent(in)    :: param  !! Current run configuration parameters of 
         real(DP),                     intent(in)    :: dt     !! Stepsize
      end subroutine helio_drift_tp
   
      module subroutine helio_drift_linear_pl(self, cb, dt, lbeg)
         implicit none
         class(helio_pl), intent(inout) :: self  !! Helio massive body object
         class(helio_cb), intent(inout) :: cb    !! Helio central body object
         real(DP),        intent(in)    :: dt    !! Stepsize
         logical,         intent(in)    :: lbeg  !! Argument that determines whether or not this is the beginning or end of the step
      end subroutine helio_drift_linear_pl 

      module subroutine helio_drift_linear_tp(self, cb, dt, lbeg)
         implicit none
         class(helio_tp), intent(inout) :: self !! Helio test particle object
         class(helio_cb), intent(in)    :: cb   !! Helio nbody system object
         real(DP),        intent(in)    :: dt   !! Stepsize
         logical,         intent(in)    :: lbeg !! Argument that determines whether or not this is the beginning or end of the step
      end subroutine helio_drift_linear_tp

      module subroutine helio_getacch_pl(self, system, param, t, lbeg)
         use swiftest_classes, only : swiftest_parameters, swiftest_nbody_system
         implicit none
         class(helio_pl),              intent(inout) :: self   !! Helio massive body object
         class(swiftest_nbody_system), intent(inout) :: system !! Swiftest nbody system object
         class(swiftest_parameters),   intent(in)    :: param  !! Current run configuration parameters of 
         real(DP),                     intent(in)    :: t      !! Current simulation time
         logical, optional,            intent(in)    :: lbeg   !! Optional argument that determines whether or not this is the beginning or end of the step
      end subroutine helio_getacch_pl

      module subroutine helio_getacch_tp(self, system, param, t, lbeg)
         use swiftest_classes, only : swiftest_parameters, swiftest_nbody_system
         implicit none
         class(helio_tp),              intent(inout) :: self   !! Helio test particle object
         class(swiftest_nbody_system), intent(inout) :: system !! Swiftest nbody system object
         class(swiftest_parameters),   intent(in)    :: param  !! Current run configuration parameters
         real(DP),                     intent(in)    :: t      !! Current time
         logical, optional,            intent(in)    :: lbeg   !! Optional argument that determines whether or not this is the beginning or end of the step
      end subroutine helio_getacch_tp

      module subroutine helio_kickvb_pl(self, dt)
         implicit none
         class(helio_pl), intent(inout) :: self !! Helio massive body object
         real(DP),        intent(in)    :: dt   !! Stepsize
      end subroutine helio_kickvb_pl

      module subroutine helio_kickvb_tp(self, dt)
         implicit none
         class(helio_tp), intent(inout) :: self !! Helio test particle object
         real(DP),        intent(in)    :: dt   !! Stepsize
      end subroutine helio_kickvb_tp

      module subroutine helio_step_system(self, param, t, dt)
         use swiftest_classes, only : swiftest_parameters
         implicit none
         class(helio_nbody_system),  intent(inout) :: self   !! Helio nbody system object
         class(swiftest_parameters), intent(inout) :: param  !! Current run configuration parameters
         real(DP),                   intent(in)    :: t      !! Simulation time
         real(DP),                   intent(in)    :: dt     !! Current stepsize 
      end subroutine helio_step_system

      module subroutine helio_step_pl(self, system, param, t, dt)
         use swiftest_classes, only : swiftest_nbody_system, swiftest_parameters
         implicit none
         class(helio_pl),              intent(inout) :: self   !! Helio massive body particle object
         class(swiftest_nbody_system), intent(inout) :: system !! Swiftest nboody system
         class(swiftest_parameters),   intent(inout) :: param  !! Current run configuration parameters 
         real(DP),                     intent(in)    :: t      !! Current simulation time
         real(DP),                     intent(in)    :: dt     !! Stepsize
      end subroutine helio_step_pl

      module subroutine helio_step_tp(self, system, param, t, dt)
         use swiftest_classes, only : swiftest_cb, swiftest_parameters, swiftest_nbody_system
         implicit none
         class(helio_tp),              intent(inout) :: self   !! Helio test particle object
         class(swiftest_nbody_system), intent(inout) :: system !! Swiftest nbody system object
         class(swiftest_parameters),   intent(inout) :: param  !! Current run configuration parameters 
         real(DP),                     intent(in)    :: t      !! Current simulation time
         real(DP),                     intent(in)    :: dt     !! Stepsizee
      end subroutine helio_step_tp
   end interface
end module helio_classes
