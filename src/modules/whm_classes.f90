module whm_classes
   !! author: David A. Minton
   !!
   !! Definition of classes and methods specific to the Democratic Heliocentric Method
   !! Partially adapted from David E. Kaufmann's Swifter module: module_whm.f90
   use swiftest_globals
   use swiftest_classes
   implicit none

   !********************************************************************************************************************************
   !                                    whm_pl class definitions and method interfaces
   !*******************************************************************************************************************************

   !! WHM massive body particle class
   type, public, extends(swiftest_pl) :: whm_pl
      real(DP), dimension(:),   allocatable :: eta     ! Jacobi mass
      real(DP), dimension(:,:), allocatable :: xj      ! Jacobi position
      real(DP), dimension(:,:), allocatable :: vj      ! Jacobi velocity
      real(DP), dimension(:,:), allocatable :: ah1     ! First term of heliocentric acceleration
      real(DP), dimension(:,:), allocatable :: ah2     ! Second term of heliocentric acceleration
      real(DP), dimension(:,:), allocatable :: ah3     ! Third term of heliocentric acceleration
      real(DP), dimension(:,:), allocatable :: ah      ! Total heliocentric acceleration
   contains
      procedure, public :: alloc => whm_allocate_pl     !! Constructor method - Allocates space for number of particles
      procedure, public :: getacch => whm_getacch_pl    !! Compute heliocentric accelerations of massive bodies
      procedure, public :: step => whm_step_pl          !! Step massive bodies ahead Democratic Heliocentric method
      procedure, public :: drift => whm_drift_pl        !! Loop through massive bodies and call Danby drift routine
      procedure, public :: kick => whm_kickvh_pl        !! Kick barycentric velocities of active massive bodies
   end type whm_pl

   interface
      !> WHM massive body constructor method
      module subroutine whm_allocate_pl(self,n)
         implicit none
         class(whm_pl), intent(inout)    :: self !! Swiftest test particle object
         integer, intent(in)               :: n    !! Number of test particles to allocate
      end subroutine whm_allocate_pl

      !> Get heliocentric accelration of massive bodies
      module subroutine whm_getacch_pl(self, config, t, lflag, whm_plA, xh)
         implicit none
         class(whm_pl), intent(inout)                  :: self      !! WHM massive body particle data structure. 
         type(swiftest_configuration),intent(in)       :: config    !! Input collection of user-defined parameter
         real(DP), intent(in)                          :: t         !! Current time. This is passed to the user-defined acceleration function.
         logical, intent(in)                           :: lflag     !! Logical flag indicating whether to recompute direct cross term accelrations
         class(whm_pl), optional, intent(in)           :: whm_plA !! Dummy argument used to make this a polymorphic method for both pl and tp classes
         real(DP), dimension(:,:), optional,intent(in) :: xh        !! Dummy argument to make this method polymorphic with the tp class
      end subroutine whm_getacch_pl

      module subroutine whm_drift_pl(self, dt, mu)
         implicit none
         class(whm_pl), intent(inout)                  :: self   !! WHM massive body particle data structure
         real(DP), intent(in)                          :: dt     !! Stepsize
         real(DP), optional, intent(in)                :: mu     !! G * (m1 + m2), G = gravitational constant, m1 = mass of central body, m2 = mass of body to drift
      end subroutine whm_drift_pl

      module subroutine whm_getacch_int_pl(self)
         implicit none
         class(whm_pl), intent(inout)                  :: self   !! WHM massive body particle data structure.
      end subroutine whm_getacch_int_pl

      module subroutine whm_kickvh_pl(self, dt)
         implicit none
         class(whm_pl), intent(inout)                   :: self  !! WHM massive body particle data structure.
         real(DP), intent(in)                           :: dt    !! Stepsize
      end subroutine whm_kickvh_pl

      module subroutine whm_step_pl(self, config, t, dt)
         implicit none
         class(whm_pl), intent(inout)                   :: self   !! WHM massive body particle data structure.
         type(swiftest_configuration),intent(in)        :: config !! Input collection of user-defined parameter
         real(DP), intent(in)                           :: t      !! Current time
         real(DP), intent(in)                           :: dt     !! Stepsize
      end subroutine whm_step_pl

      module subroutine whm_user_getacch_pl(self, t)
         implicit none
         class(whm_pl), intent(inout)                   :: self   !! WHM massive body particle data structure
         real(DP), intent(in)                           :: t      !! Current time
      end subroutine whm_user_getacch_pl
   end interface

   !********************************************************************************************************************************
   !                                    whm_tp class definitions and method interfaces
   !*******************************************************************************************************************************

   !! WHM test particle class
   type, public, extends(swiftest_tp) :: whm_tp
      real(DP), dimension(:,:), allocatable :: ah  !! Total heliocentric acceleration
   contains
      private
      procedure, public     :: alloc    => whm_allocate_tp  !! Allocates new components of the whm class and recursively calls parent allocations
      procedure, public     :: getacch  => whm_getacch_tp   !! Compute heliocentric accelerations of test particles
      procedure, public     :: step => whm_step_tp          !! Step active test particles ahead using Democratic Heliocentric method
      procedure, public     :: drift => whm_drift_tp        !! Loop through test particles and call Danby drift routine
      procedure, public     :: kick => whm_kickvh_tp        !! Kick barycentric velocities of active test particles
   end type whm_tp

   interface
      !> WHM test particle constructor 
      module subroutine whm_allocate_tp(self,n)
         implicit none
         class(whm_tp), intent(inout)           :: self !! Swiftest test particle object
         integer, intent(in)                    :: n    !! Number of test particles to allocate
      end subroutine whm_allocate_tp

      module subroutine whm_drift_tp(self, dt, mu)
         implicit none
         class(whm_tp), intent(inout)           :: self !! WHM test particle data structure
         real(DP), intent(in)                   :: dt   !! Stepsize
         real(DP), optional, intent(in)         :: mu   !! G * m1, G = gravitational constant, m1 = mass of central body
      end subroutine whm_drift_tp

      !> Get heliocentric accelration of the test particle
      module subroutine whm_getacch_tp(self, pl, config, t, xh)
         implicit none
         class(whm_tp), intent(inout)            :: self   !! WHM test particle data structure
         type(whm_pl),  intent(in)               :: pl     !! WHM massive body particle data structure. 
         type(swiftest_configuration),intent(in) :: config !! Input collection of user-defined parameter
         real(DP), intent(in)                    :: t      !! Current time. This is passed to the user-defined acceleration function.
         real(DP), dimension(:,:), intent(in)    :: xh     !! Heliocentric positions of massive bodies at time t
      end subroutine whm_getacch_tp

      module subroutine whm_kickvh_tp(self, dt)
         implicit none
         class(whm_tp), intent(inout)            :: self !! WHM test particle data structure
         real(DP), intent(in)                    :: dt   !! Stepsize
      end subroutine whm_kickvh_tp

      module subroutine whm_step_tp(self, pl, config, t, dt,  xbeg, xend)
         implicit none
         class(whm_tp), intent(inout)              :: self      !! WHM test particle data structure.
         class(whm_pl), intent(in)                 :: pl        !! WHM massive body particle data structure.
         type(swiftest_configuration),intent(in)   :: config    !! Input collection of user-defined parameter
         real(DP), intent(in)                      :: t         !! Current time
         real(DP), intent(in)                      :: dt        !! Stepsize
         real(DP), dimension(:, :), intent(inout)  :: xbeg      !! Heliocentric massive body positions at beginning of time step
         real(DP), dimension(:, :), intent(inout)  :: xend      !! Heliocentric massive body positions at end of time step
      end subroutine whm_step_tp

      module subroutine whm_user_getacch_tp(self,t)
         implicit none
         class(whm_tp), intent(inout)              :: self      !! WHM massive body particle data structure
         real(DP), intent(in)                      :: t          !! Current time
      end subroutine whm_user_getacch_tp

   end interface

   !********************************************************************************************************************************
   !                            whm_nbody_system class definitions and method interfaces
   !********************************************************************************************************************************
   !> An abstract class for the WHM integrator nbody system 
   type, public, extends(swiftest_nbody_system) :: whm_nbody_system
      !> In the WHM integrator, only test particles are discarded
      class(swiftest_tp), allocatable :: tp_discards
   contains
      private
      !> Replace the abstract procedures with concrete ones
      procedure, public :: construct     => whm_construct_system   !! Perform a discard operation and spill any discarded bodies to list for output.  
      procedure, public :: step          => whm_step               !! Method to advance the system one step in time given by the step size dt
   end type whm_nbody_system

!> Interfaces for all non-type bound whm methods that are implemented in separate submodules 
interface
   !> Constructs a WHM nbody system
   module subroutine whm_construct_system(self, config, integrator)
      implicit none
      class(whm_nbody_system),       intent(inout) :: self       !! Swiftest system object
      class(swiftest_configuration), intent(out)   :: config     !! Input collection of user-defined configuration parameters
      integer, intent(in)                          :: integrator !! Integrator type code
   end subroutine whm_construct_system

   !> Steps the Swiftest nbody system forward in time one stepsize
   module subroutine whm_step(self, config, t, dt) 
      implicit none
      class(whm_nbody_system),       intent(inout) :: self    !! Swiftest system object
      class(swiftest_configuration), intent(in)    :: config  !! Input collection of user-defined configuration parameters
      real(DP),                      intent(in)    :: t       !! Current simulation time
      real(DP),                      intent(in)    :: dt      !! Stepsize
   end subroutine whm_step
end interface


end module whm_classes
