module helio
   !! author: The Purdue Swiftest Team - David A. Minton, Carlisle A. Wishard, Jennifer L.L. Pouplin, and Jacob R. Elliott
   !!
   !! Definition of classes and methods specific to the Democratic Heliocentric Method
   !! Adapted from David E. Kaufmann's Swifter modules: helio.f90
   use swiftest_globals
   use swiftest_data_structures
   implicit none

   !********************************************************************************************************************************
   !                                    helio_tp class definitions and method interfaces
   !*******************************************************************************************************************************

   !! Helio test particle class
   type, public, extends(swiftest_pl) :: helio_tp
      real(DP), dimension(:,:), allocatable :: ah  !! Total heliocentric acceleration
      real(DP), dimension(:,:), allocatable :: ahi !! Heliocentric acceleration due to interactions
   contains
      private
      procedure, public     :: alloc => helio_allocate_tp     !! Allocates new components of the helio class and recursively calls parent allocations
      procedure, public     :: set_vec => swiftest_set_vec    !! Method used to construct the vectorized form of the central body mass
      !procedure, public     :: spill => helio_spill_tp        
      procedure, public     :: getacch => helio_getacch_tp    !! Compute heliocentric accelerations of test particles
      procedure, public     :: step => helio_step_tp          !! Step active test particles ahead using Democratic Heliocentric method
      procedure, public     :: drift => helio_drift_tp        !! Loop through test particles and call Danby drift routine
      procedure, public     :: lindrift => helio_lindrift_tp  !! Perform linear drift of test particles due to barycentric momentum of Sun
      procedure, public     :: kick => helio_kickvb_tp        !! Kick barycentric velocities of active test particles
      final :: helio_deallocate_tp
   end type helio_tp

   interface
      !> Helio constructor and desctructor methods
      module subroutine helio_allocate_tp(self,n)
         implicit none
         class(helio_tp), intent(inout) :: self !! Swiftest test particle object
         integer, intent(in)            :: n    !! Number of test particles to allocate
      end subroutine helio_allocate_tp

      !> Helio test particle destructor/finalizer
      module subroutine helio_deallocate_tp(self)
         implicit none
         type(helio_tp), intent(inout)  :: self
      end subroutine helio_deallocate_tp

      !> Method to remove the inactive bodies and spill them to a discard object
      module subroutine helio_spill_tp(self,discard)
         implicit none
         class(helio_tp), intent(inout) :: self    !! Swiftest test particle object to input
         class(helio_tp), intent(inout) :: discard !! Discarded body list
      end subroutine helio_spill_tp
  
   end interface

   !********************************************************************************************************************************
   !                                    helio_pl class definitions and method interfaces
   !*******************************************************************************************************************************

   !! Helio massive body particle class
   type, public, extends(helio_tp) :: helio_pl
   contains
      procedure, public :: alloc => helio_allocate_pl     !! Constructor method - Allocates space for number of particles
      !procedure, public :: set_from_file => io_read_pl_in !! Override helio_tp io reader with the pl reader
      procedure, public :: getacch => helio_getacch_pl    !! Compute heliocentric accelerations of massive bodies
      procedure, public :: step => helio_step_pl          !! Step massive bodies ahead Democratic Heliocentric method
      procedure, public :: drift => helio_drift_pl        !! Loop through massive bodies and call Danby drift routine
      procedure, public :: lindrift => helio_lindrift_pl  !! Perform linear drift of massive bodies due to barycentric momentum of Sun
      procedure, public :: kick => helio_kickvb_pl        !! Kick barycentric velocities of active massive bodies
      final :: helio_deallocate_pl                        !! Finalizer method - Deallocates all allocatables 
   end type helio_pl

   interface
      !> Helio massive body constructor method
      module subroutine helio_allocate_pl(self,n)
         implicit none
         class(helio_pl), intent(inout)    :: self !! Swiftest test particle object
         integer, intent(in)               :: n    !! Number of test particles to allocate
      end subroutine helio_allocate_pl

      !> Helio massive body destructor/finalizer
      module subroutine helio_deallocate_pl(self)
         implicit none
         type(helio_pl), intent(inout)    :: self
      end subroutine helio_deallocate_pl
   end interface

!> Interfaces for all non-type bound helio methods that are implemented in separate submodules 
interface
   !> Call discard routine to determine spilled test particles, then remove them from active list
   module function helio_discard(helio_plA, helio_tpA, config, t, dt) result(helio_tp_discard)
      implicit none
      type(helio_pl), intent(inout)           :: helio_plA !! Helio massive body particle data structure. 
      type(helio_tp), intent(inout)           :: helio_tpA !! Helio test particle data structure
      type(swiftest_configuration),intent(in) :: config    !! Input collection of user-defined parameter
      real(DP), intent(in)                    :: t         !! Current time
      real(DP), intent(in)                    :: dt        !! Stepsize
      type(helio_tp)                          :: helio_tp_discard !! Discarded Helio test particles
   end function helio_discard

   !> Get helioctric accelration of massive bodies
   module subroutine helio_getacch_pl(self, config, t, lflag, helio_plA, xh)
      implicit none
      class(helio_pl), intent(inout)                :: self      !! Helio massive body particle data structure. 
      type(swiftest_configuration),intent(in)       :: config    !! Input collection of user-defined parameter
      real(DP), intent(in)                          :: t         !! Current time. This is passed to the user-defined acceleration function.
      logical, intent(in)                           :: lflag     !! Logical flag indicating whether to recompute direct cross term accelrations
      class(helio_pl), optional, intent(inout)      :: helio_plA !! Dummy argument used to make this a polymorphic method for both pl and tp classes
      real(DP), dimension(:,:), optional,intent(in) :: xh        !! Dummy argument to make this method polymorphic with the tp class
   end subroutine helio_getacch_pl

   !> Get heliocentric accelration of the test particle
   module subroutine helio_getacch_tp(self, config, t, lflag, helio_plA, xh)
      implicit none
      class(helio_tp), intent(inout)                 :: self      !! Helio test particle data structure
      type(swiftest_configuration),intent(in)        :: config    !! Input collection of user-defined parameter
      real(DP), intent(in)                           :: t         !! Current time. This is passed to the user-defined acceleration function.
      logical, intent(in)                            :: lflag     !! Logical flag indicating whether to recompute direct cross term accelrations
      class(helio_pl), optional, intent(inout)       :: helio_plA !! Helio massive body particle data structure. Optional to allow this method to be polymorphic for pl and tp classes
      real(DP), dimension(:,:), optional, intent(in) :: xh        !! Heliocentric positions of massive bodies at time t
   end subroutine helio_getacch_tp

   module subroutine helio_drift_pl(self, dt, mu)
      implicit none
      class(helio_pl), intent(inout) :: self !! Helio massive body particle data structure
      real(DP), intent(in)           :: dt   !! Stepsize
      real(DP), optional, intent(in) :: mu   !! G * (m1 + m2), G = gravitational constant, m1 = mass of central body, m2 = mass of body to drift
   end subroutine helio_drift_pl

   module subroutine helio_drift_tp(self, dt, mu)
      implicit none
      class(helio_tp), intent(inout) :: self !! Helio test particle data structure
      real(DP), intent(in)           :: dt   !! Stepsize
      real(DP), optional, intent(in) :: mu   !! G * m1, G = gravitational constant, m1 = mass of central body
   end subroutine helio_drift_tp

   module subroutine helio_getacch_int_pl(self, helio_plA)
      implicit none
      class(helio_pl), intent(inout)           :: self      !! Helio massive body particle data structure.
      class(helio_pl), optional, intent(inout) :: helio_plA !! Dummy argument used to make this a polymorphic method for both pl and tp classes
   end subroutine helio_getacch_int_pl

   module subroutine helio_getacch_int_tp(self, helio_plA)
      implicit none
      class(helio_tp), intent(inout)           :: self      !! Helio test particle data structure
      class(helio_pl), optional, intent(inout) :: helio_plA !! Helio massive body particle data structure. Optional to allow this method to be polymorphic for pl and tp classes
   end subroutine helio_getacch_int_tp

   module subroutine helio_kickvb_pl(self, dt)
      implicit none
      class(helio_pl), intent(inout) :: self !! Helio massive body particle data structure.
      real(DP), intent(in)           :: dt   !! Stepsize
   end subroutine helio_kickvb_pl

   module subroutine helio_kickvb_tp(self, dt)
      implicit none
      class(helio_tp), intent(inout) :: self !! Helio test particle data structure
      real(DP), intent(in)           :: dt   !! Stepsize
   end subroutine helio_kickvb_tp

   module subroutine helio_lindrift_pl(self, dt, pt)
      implicit none
      class(helio_pl), intent(inout)           :: self !! Helio massive body particle data structur
      real(DP), intent(in)                     :: dt   !! Stepsize
      real(DP), dimension(NDIM), intent(inout) :: pt   !! Negative barycentric velocity of the Sun
   end subroutine helio_lindrift_pl

   module subroutine helio_lindrift_tp(self, dt, pt)
      implicit none
      class(helio_tp), intent(inout)           :: self !! Helio test particle data structure
      real(DP), intent(in)                     :: dt   !! Stepsize
      real(DP), dimension(NDIM), intent(inout) :: pt   !! Negative barycentric velocity of the Sun
   end subroutine helio_lindrift_tp

   module subroutine helio_step(helio_plA, helio_tpA, config, t, dt, lfirst)
      implicit none
      type(helio_pl), intent(inout)           :: helio_plA !! Helio massive body particle data structure. 
      type(helio_tp), intent(inout)           :: helio_tpA !! Helio test particle data structure
      type(swiftest_configuration),intent(in) :: config    !! Input collection of user-defined parameter
      real(DP), intent(in)                    :: t         !! Current time
      real(DP), intent(in)                    :: dt        !! Stepsize
      logical, intent(inout)                  :: lfirst    !! Logical flag indicating whether current invocation is the first 
                                                           !!     TODO: Replace lfirst with internal flag with save attribute
   end subroutine helio_step

   module subroutine helio_step_pl(self, helio_plA, config, t, dt, lfirst, xbeg, xend, ptb, pte)
      implicit none
      class(helio_pl), intent(inout)           :: self            !! Helio massive body particle data structure.
      class(helio_pl), optional, intent(in)    :: helio_plA       !! Dummy argument used to make this a polymorphic method for both pl and tp classes
      type(swiftest_configuration),intent(in)  :: config          !! Input collection of user-defined parameter
      real(DP), intent(in)                     :: t               !! Current time
      real(DP), intent(in)                     :: dt              !! Stepsize
      logical, intent(inout)                   :: lfirst          !! Logical flag indicating whether current invocation is the first 
                                                                  !!     TODO: Replace lfirst with internal flag with save attribute
      real(DP), dimension(:), intent(inout)    :: ptb             !! Negative barycentric velocity of the Sun at beginning of time step
      real(DP), dimension(:), intent(inout)    :: pte             !! Negative barycentric velocity of the Sun at end of time step
      real(DP), dimension(:, :), intent(inout) :: xbeg            !! Heliocentric massive body positions at beginning of time step
      real(DP), dimension(:, :), intent(inout) :: xend            !! Heliocentric massive body positions at end of time step
   end subroutine helio_step_pl

   module subroutine helio_step_tp(self, helio_plA, config, t, dt, lfirst, xbeg, xend, ptb, pte)
      implicit none
      class(helio_tp), intent(inout)           :: self      !! Helio test particle data structure.
      class(helio_pl), optional, intent(in)    :: helio_plA !! Helio massive body particle data structure.
      type(swiftest_configuration),intent(in)  :: config    !! Input collection of user-defined parameter
      real(DP), intent(in)                     :: t         !! Current time
      real(DP), intent(in)                     :: dt        !! Stepsize
      logical, intent(inout)                   :: lfirst    !! Logical flag indicating whether current invocation is the first 
                                                            !!     TODO: Replace lfirst with internal flag with save attribute
      real(DP), dimension(:), intent(inout)    :: ptb       !! Negative barycentric velocity of the Sun at beginning of time step
      real(DP), dimension(:), intent(inout)    :: pte       !! Negative barycentric velocity of the Sun at end of time step
      real(DP), dimension(:, :), intent(inout) :: xbeg      !! Heliocentric massive body positions at beginning of time step
      real(DP), dimension(:, :), intent(inout) :: xend      !! Heliocentric massive body positions at end of time step
   end subroutine helio_step_tp

   module subroutine helio_user_getacch_pl(helio_plA, t)
      implicit none
      class(helio_pl), intent(inout) :: helio_plA  !! Helio massive body particle data structure
      real(DP), intent(in)           :: t          !! Current time
   end subroutine helio_user_getacch_pl

   module subroutine helio_user_getacch_tp(helio_tpA, t)
      implicit none
      class(helio_tp), intent(inout) :: helio_tpA  !! Helio massive body particle data structure
      real(DP), intent(in)           :: t          !! Current tim
   end subroutine helio_user_getacch_tp
end interface

!> Only the constructor and destructor method implementations are listed here. All other methods are implemented in the helio submodules.
   contains





      

end module helio
