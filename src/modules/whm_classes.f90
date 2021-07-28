module whm_classes
   !! author: David A. Minton
   !!
   !! Definition of classes and methods specific to the Democratic Heliocentric Method
   !! Partially adapted from David E. Kaufmann's Swifter module: module_whm.f90
   use swiftest_globals
   use swiftest_classes, only : swiftest_cb, swiftest_pl, swiftest_tp, swiftest_nbody_system
   implicit none
   public

   !********************************************************************************************************************************
   ! whm_cb class definitions and method interfaces
   !*******************************************************************************************************************************
   !> Swiftest central body particle class
   type, extends(swiftest_cb) :: whm_cb
   contains
   end type whm_cb

   !********************************************************************************************************************************
   !                                    whm_pl class definitions and method interfaces
   !*******************************************************************************************************************************

   !> WHM massive body particle class
   type, extends(swiftest_pl) :: whm_pl
      real(DP), dimension(:),   allocatable :: eta    !! Jacobi mass
      real(DP), dimension(:,:), allocatable :: xj     !! Jacobi position
      real(DP), dimension(:,:), allocatable :: vj     !! Jacobi velocity
      real(DP), dimension(:),   allocatable :: muj    !! Jacobi mu: GMcb * eta(i) / eta(i - 1) 
      real(DP), dimension(:),   allocatable :: ir3j    !! Third term of heliocentric acceleration
      !! Note to developers: If you add componenets to this class, be sure to update methods and subroutines that traverse the
      !!    component list, such as whm_setup_pl and whm_util_spill_pl
   contains
      procedure :: h2j         => whm_coord_h2j_pl       !! Convert position and velcoity vectors from heliocentric to Jacobi coordinates 
      procedure :: j2h         => whm_coord_j2h_pl       !! Convert position and velcoity vectors from Jacobi to helliocentric coordinates 
      procedure :: vh2vj       => whm_coord_vh2vj_pl     !! Convert velocity vectors from heliocentric to Jacobi coordinates 
      procedure :: drift       => whm_drift_pl           !! Loop through massive bodies and call Danby drift routine to jacobi coordinates
      procedure :: fill        => whm_util_fill_pl       !! "Fills" bodies from one object into another depending on the results of a mask (uses the MERGE intrinsic)
      procedure :: accel       => whm_kick_getacch_pl    !! Compute heliocentric accelerations of massive bodies
      procedure :: kick        => whm_kick_vh_pl         !! Kick heliocentric velocities of massive bodies
      procedure :: accel_gr    => whm_gr_kick_getacch_pl !! Acceleration term arising from the post-Newtonian correction
      procedure :: gr_pos_kick => whm_gr_p4_pl           !! Position kick due to p**4 term in the post-Newtonian correction
      procedure :: setup       => whm_setup_pl           !! Constructor method - Allocates space for number of particles
      procedure :: set_mu      => whm_util_set_mu_eta_pl !! Sets the Jacobi mass value for all massive bodies.
      procedure :: set_ir3     => whm_util_set_ir3j     !! Sets both the heliocentric and jacobi inverse radius terms (1/rj**3 and 1/rh**3)
      procedure :: step        => whm_step_pl            !! Steps the body forward one stepsize
      procedure :: spill       => whm_util_spill_pl      !!"Spills" bodies from one object to another depending on the results of a mask (uses the PACK intrinsic)
   end type whm_pl

   !********************************************************************************************************************************
   !  whm_tp class definitions and method interfaces
   !*******************************************************************************************************************************

   !! WHM test particle class
   type, extends(swiftest_tp) :: whm_tp
      !! Note to developers: If you add componenets to this class, be sure to update methods and subroutines that traverse the
      !!    component list, such as whm_setup_tp and whm_util_spill_tp
   contains
      procedure :: accel       => whm_kick_getacch_tp    !! Compute heliocentric accelerations of test particles
      procedure :: kick        => whm_kick_vh_tp         !! Kick heliocentric velocities of test particles
      procedure :: accel_gr    => whm_gr_kick_getacch_tp !! Acceleration term arising from the post-Newtonian correction
      procedure :: gr_pos_kick => whm_gr_p4_tp           !! Position kick due to p**4 term in the post-Newtonian correction
      procedure :: setup       => whm_setup_tp           !! Allocates new components of the whm class and recursively calls parent allocations
      procedure :: step        => whm_step_tp            !! Steps the particle forward one stepsize
   end type whm_tp

   !********************************************************************************************************************************
   !  whm_nbody_system class definitions and method interfaces
   !********************************************************************************************************************************
   !> An abstract class for the WHM integrator nbody system 
   type, extends(swiftest_nbody_system) :: whm_nbody_system
   contains
      !> Replace the abstract procedures with concrete ones
      procedure :: initialize   => whm_setup_initialize_system      !! Performs WHM-specific initilization steps, like calculating the Jacobi masses
      procedure :: step         => whm_step_system       !! Advance the WHM nbody system forward in time by one step
   end type whm_nbody_system

   interface
      module subroutine whm_coord_h2j_pl(self, cb)
         use swiftest_classes, only : swiftest_cb
         implicit none
         class(whm_pl),      intent(inout) :: self   !! WHM massive body particle data structure
         class(swiftest_cb), intent(inout) :: cb     !! Swiftest central body particle data structuree
      end subroutine whm_coord_h2j_pl

      module subroutine whm_coord_j2h_pl(self, cb)
         use swiftest_classes, only : swiftest_cb
         implicit none
         class(whm_pl),      intent(inout) :: self   !! WHM massive body particle data structure
         class(swiftest_cb), intent(inout) :: cb     !! Swiftest central body particle data structuree
      end subroutine whm_coord_j2h_pl

      module subroutine whm_coord_vh2vj_pl(self, cb)
         use swiftest_classes, only : swiftest_cb
         implicit none
         class(whm_pl),      intent(inout) :: self   !! WHM massive body particle data structure
         class(swiftest_cb), intent(inout) :: cb     !! Swiftest central body particle data structuree
      end subroutine whm_coord_vh2vj_pl

      module subroutine whm_drift_pl(self, system, param, dt, mask)
         use swiftest_classes, only : swiftest_nbody_system, swiftest_parameters
         implicit none
         class(whm_pl),                intent(inout) :: self   !! WHM massive body particle data structure
         class(swiftest_nbody_system), intent(inout) :: system !! WHM nbody system object
         class(swiftest_parameters),   intent(in)    :: param  !! Current run configuration parameters 
         real(DP),                     intent(in)    :: dt     !! Stepsize
         logical, dimension(:),        intent(in)    :: mask   !! Logical mask of size self%nbody that determines which bodies to drift
      end subroutine whm_drift_pl

      module subroutine whm_util_fill_pl(self, inserts, lfill_list)
         use swiftest_classes, only : swiftest_body
         implicit none
         class(whm_pl),         intent(inout) :: self       !! WHM massive body object
         class(swiftest_body),  intent(inout) :: inserts    !! inserted object 
         logical, dimension(:), intent(in)    :: lfill_list !! Logical array of bodies to merge into the keeps
      end subroutine whm_util_fill_pl

      !> Get heliocentric accelration of massive bodies
      module subroutine whm_kick_getacch_pl(self, system, param, t, lbeg)
         use swiftest_classes, only : swiftest_cb, swiftest_parameters
         implicit none
         class(whm_pl),                intent(inout) :: self   !! WHM massive body particle data structure
         class(swiftest_nbody_system), intent(inout) :: system !! WHM nbody system object
         class(swiftest_parameters),   intent(in)    :: param  !! Current run configuration parameters 
         real(DP),                     intent(in)    :: t      !! Current simulation time
         logical,                      intent(in)    :: lbeg   !! Logical flag that determines whether or not this is the beginning or end of the step
      end subroutine whm_kick_getacch_pl

      !> Get heliocentric accelration of the test particle
      module subroutine whm_kick_getacch_tp(self, system, param, t, lbeg)
         use swiftest_classes, only : swiftest_cb, swiftest_parameters
         implicit none
         class(whm_tp),                intent(inout) :: self   !! WHM test particle data structure
         class(swiftest_nbody_system), intent(inout) :: system !! WHM nbody system object
         class(swiftest_parameters),   intent(in)    :: param  !! Current run configuration parameters
         real(DP),                     intent(in)    :: t      !! Current time
         logical,                      intent(in)    :: lbeg   !! Logical flag that determines whether or not this is the beginning or end of the step
      end subroutine whm_kick_getacch_tp

      module subroutine whm_kick_vh_pl(self, system, param, t, dt, mask, lbeg)
         use swiftest_classes, only : swiftest_nbody_system, swiftest_parameters
         implicit none
         class(whm_pl),                intent(inout) :: self   !! WHM massive body object
         class(swiftest_nbody_system), intent(inout) :: system !! Swiftest nbody system object
         class(swiftest_parameters),   intent(in)    :: param  !! Current run configuration parameters 
         real(DP),                     intent(in)    :: t      !! Current time
         real(DP),                     intent(in)    :: dt     !! Stepsize
         logical, dimension(:),        intent(in)    :: mask   !! Mask that determines which bodies to kick
         logical,                      intent(in)    :: lbeg   !! Logical flag indicating whether this is the beginning of the half step or not. 
      end subroutine whm_kick_vh_pl

      module subroutine whm_kick_vh_tp(self, system, param, t, dt, mask, lbeg)
         use swiftest_classes, only : swiftest_nbody_system, swiftest_parameters
         implicit none
         class(whm_tp),                intent(inout) :: self   !! WHM test particle object
         class(swiftest_nbody_system), intent(inout) :: system !! Swiftest nbody system object
         class(swiftest_parameters),   intent(in)    :: param  !! Current run configuration parameters 
         real(DP),                     intent(in)    :: t      !! Current time
         real(DP),                     intent(in)    :: dt     !! Stepsize
         logical, dimension(:),        intent(in)    :: mask   !! Mask that determines which bodies to kick
         logical,                      intent(in)    :: lbeg   !! Logical flag indicating whether this is the beginning of the half step or not. 
      end subroutine whm_kick_vh_tp

      module subroutine whm_gr_kick_getacch_pl(self, param)
         use swiftest_classes, only : swiftest_cb, swiftest_parameters
         implicit none
         class(whm_pl),              intent(inout) :: self  !! WHM massive body particle data structure
         class(swiftest_parameters), intent(in)    :: param !! Current run configuration parameters 
      end subroutine whm_gr_kick_getacch_pl

      module subroutine whm_gr_kick_getacch_tp(self, param)
         use swiftest_classes, only : swiftest_cb, swiftest_parameters
         implicit none
         class(whm_tp),              intent(inout) :: self  !! WHM test particle data structure
         class(swiftest_parameters), intent(in)    :: param !! Current run configuration parameters
      end subroutine whm_gr_kick_getacch_tp

      module pure subroutine whm_gr_p4_pl(self, param, dt)
         use swiftest_classes, only : swiftest_parameters
         implicit none
         class(whm_pl),              intent(inout) :: self  !! WHM massive body object
         class(swiftest_parameters), intent(in)    :: param !! Current run configuration parameters 
         real(DP),                   intent(in)    :: dt    !! Step size
      end subroutine whm_gr_p4_pl

      module pure subroutine whm_gr_p4_tp(self, param, dt)
         use swiftest_classes, only : swiftest_parameters
         implicit none
         class(whm_tp),              intent(inout) :: self  !! WHM test particle object
         class(swiftest_parameters), intent(in)    :: param !! Current run configuration parameters 
         real(DP),                    intent(in)   :: dt    !! Step size
      end subroutine whm_gr_p4_tp

      !> Reads WHM massive body object in from file
      module subroutine whm_setup_pl(self,n)
         implicit none
         class(whm_pl), intent(inout)    :: self !! WHM massive body objectobject
         integer(I4B),  intent(in)       :: n    !! Number of test particles to allocate
      end subroutine whm_setup_pl

      module subroutine whm_util_set_ir3j(self)
         implicit none
         class(whm_pl),                intent(inout) :: self    !! WHM massive body object
      end subroutine whm_util_set_ir3j

      module subroutine whm_util_set_mu_eta_pl(self, cb)
         use swiftest_classes, only : swiftest_cb
         implicit none
         class(whm_pl),                intent(inout) :: self    !! WHM massive body object
         class(swiftest_cb),           intent(inout) :: cb     !! Swiftest central body object
      end subroutine whm_util_set_mu_eta_pl

      module subroutine whm_setup_initialize_system(self, param)
         use swiftest_classes, only : swiftest_parameters
         implicit none
         class(whm_nbody_system),    intent(inout) :: self   !! WHM nbody system object
         class(swiftest_parameters), intent(inout) :: param  !! Current run configuration parameters 
      end subroutine whm_setup_initialize_system

      !> Reads WHM test particle object in from file
      module subroutine whm_setup_tp(self,n)
         implicit none
         class(whm_tp), intent(inout) :: self   !! WHM test particle data structure
         integer,       intent(in)    :: n      !! Number of test particles to allocate
      end subroutine whm_setup_tp

      module subroutine whm_step_pl(self, system, param, t, dt)
         use swiftest_classes, only : swiftest_nbody_system, swiftest_parameters
         implicit none
         class(whm_pl),                intent(inout) :: self   !! WHM massive body object
         class(swiftest_nbody_system), intent(inout) :: system !! Swiftest system object
         class(swiftest_parameters),   intent(inout) :: param  !! Current run configuration parameters 
         real(DP),                     intent(in)    :: t       !! Simulation time
         real(DP),                     intent(in)    :: dt     !! Current stepsize
      end subroutine whm_step_pl

      module subroutine whm_step_tp(self, system, param, t, dt)
         use swiftest_classes, only : swiftest_nbody_system, swiftest_parameters
         implicit none
         class(whm_tp),                intent(inout) :: self   !! WHM test particle data structure
         class(swiftest_nbody_system), intent(inout) :: system !! Swiftest nbody system object
         class(swiftest_parameters),   intent(inout) :: param  !! Current run configuration parameters
         real(DP),                     intent(in)    :: t      !! Current simulation time
         real(DP),                     intent(in)    :: dt     !! Stepsize
      end subroutine whm_step_tp

      module subroutine whm_util_spill_pl(self, discards, lspill_list)
         use swiftest_classes, only : swiftest_body
         implicit none
         class(whm_pl),         intent(inout) :: self        !! WHM massive body object
         class(swiftest_body),  intent(inout) :: discards    !! Discarded object 
         logical, dimension(:), intent(in)    :: lspill_list !! Logical array of bodies to spill into the discards
      end subroutine whm_util_spill_pl

      !> Steps the Swiftest nbody system forward in time one stepsize
      module subroutine whm_step_system(self, param, t, dt)
         use swiftest_classes, only : swiftest_parameters
         implicit none
         class(whm_nbody_system),    intent(inout) :: self    !! WHM system object
         class(swiftest_parameters), intent(inout) :: param  !! Current run configuration parameters 
         real(DP),                   intent(in)    :: t      !! Simulation time
         real(DP),                   intent(in)    :: dt     !! Current stepsize
      end subroutine whm_step_system
   end interface

end module whm_classes
