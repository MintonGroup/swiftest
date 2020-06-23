module whm
   !! author: The Purdue Swiftest Team - David A. Minton, Carlisle A. Wishard, Jennifer L.L. Pouplin, and Jacob R. Elliott
   !!
   !! Definition of classes and methods specific to the Democratic Heliocentric Method
   !! Partially adapted from David E. Kaufmann's Swifter module: module_whm.f90
   use swiftest_globals
   use swiftest_classes
   implicit none

   !********************************************************************************************************************************
   !                                    whm_tp class definitions and method interfaces
   !*******************************************************************************************************************************

   !! Helio test particle class
   type, public, extends(swiftest_tp) :: whm_tp
      real(DP), dimension(:,:), allocatable :: ah  !! Total whmcentric acceleration
   contains
      private
      procedure, public     :: alloc    => whm_allocate_tp     !! Allocates new components of the whm class and recursively calls parent allocations
      procedure, public     :: getacch => whm_getacch_tp    !! Compute whmcentric accelerations of test particles
      procedure, public     :: step => whm_step_tp          !! Step active test particles ahead using Democratic Heliocentric method
      procedure, public     :: drift => whm_drift_tp        !! Loop through test particles and call Danby drift routine
      procedure, public     :: lindrift => whm_lindrift_tp  !! Perform linear drift of test particles due to barycentric momentum of Sun
      procedure, public     :: kick => whm_kickvb_tp        !! Kick barycentric velocities of active test particles
      final :: whm_deallocate_tp
   end type whm_tp

   interface
      !> Helio constructor and desctructor methods
      module subroutine whm_allocate_tp(self,n)
         implicit none
         class(whm_tp), intent(inout) :: self !! Swiftest test particle object
         integer, intent(in)            :: n    !! Number of test particles to allocate
      end subroutine whm_allocate_tp

      !> Helio test particle destructor/finalizer
      module subroutine whm_deallocate_tp(self)
         implicit none
         type(whm_tp), intent(inout)  :: self
      end subroutine whm_deallocate_tp

      !> Method to remove the inactive bodies and spill them to a discard object
      module subroutine whm_spill_tp(self,discard)
         implicit none
         class(whm_tp), intent(inout) :: self    !! Swiftest test particle object to input
         class(whm_tp), intent(inout) :: discard !! Discarded body list
      end subroutine whm_spill_tp
  
   end interface

   !********************************************************************************************************************************
   !                                    whm_pl class definitions and method interfaces
   !*******************************************************************************************************************************

   !! Helio massive body particle class
   type, public, extends(whm_tp) :: whm_pl
   contains
      procedure, public :: alloc => whm_allocate_pl     !! Constructor method - Allocates space for number of particles
      procedure, public :: getacch => whm_getacch_pl    !! Compute whmcentric accelerations of massive bodies
      procedure, public :: step => whm_step_pl          !! Step massive bodies ahead Democratic Heliocentric method
      procedure, public :: drift => whm_drift_pl        !! Loop through massive bodies and call Danby drift routine
      procedure, public :: lindrift => whm_lindrift_pl  !! Perform linear drift of massive bodies due to barycentric momentum of Sun
      procedure, public :: kick => whm_kickvb_pl        !! Kick barycentric velocities of active massive bodies
      final :: whm_deallocate_pl                        !! Finalizer method - Deallocates all allocatables 
   end type whm_pl

   interface
      !> Helio massive body constructor method
      module subroutine whm_allocate_pl(self,n)
         implicit none
         class(whm_pl), intent(inout)    :: self !! Swiftest test particle object
         integer, intent(in)               :: n    !! Number of test particles to allocate
      end subroutine whm_allocate_pl

      !> Helio massive body destructor/finalizer
      module subroutine whm_deallocate_pl(self)
         implicit none
         type(whm_pl), intent(inout)    :: self
      end subroutine whm_deallocate_pl
   end interface

!> Interfaces for all non-type bound whm methods that are implemented in separate submodules 
interface
   !> Call discard routine to determine spilled test particles, then remove them from active list
   module function whm_discard(whm_plA, whm_tpA, config, t, dt) result(whm_tp_discard)
      implicit none
      type(whm_pl), intent(inout)           :: whm_plA !! Helio massive body particle data structure. 
      type(whm_tp), intent(inout)           :: whm_tpA !! Helio test particle data structure
      type(swiftest_configuration),intent(in) :: config    !! Input collection of user-defined parameter
      real(DP), intent(in)                    :: t         !! Current time
      real(DP), intent(in)                    :: dt        !! Stepsize
      type(whm_tp)                          :: whm_tp_discard !! Discarded Helio test particles
   end function whm_discard

   !> Get whmctric accelration of massive bodies
   module subroutine whm_getacch_pl(self, config, t, lflag, whm_plA, xh)
      implicit none
      class(whm_pl), intent(inout)                :: self      !! Helio massive body particle data structure. 
      type(swiftest_configuration),intent(in)       :: config    !! Input collection of user-defined parameter
      real(DP), intent(in)                          :: t         !! Current time. This is passed to the user-defined acceleration function.
      logical, intent(in)                           :: lflag     !! Logical flag indicating whether to recompute direct cross term accelrations
      class(whm_pl), optional, intent(in)         :: whm_plA !! Dummy argument used to make this a polymorphic method for both pl and tp classes
      real(DP), dimension(:,:), optional,intent(in) :: xh        !! Dummy argument to make this method polymorphic with the tp class
   end subroutine whm_getacch_pl

   !> Get whmcentric accelration of the test particle
   module subroutine whm_getacch_tp(self, config, t, lflag, whm_plA, xh)
      implicit none
      class(whm_tp), intent(inout)                 :: self      !! Helio test particle data structure
      type(swiftest_configuration),intent(in)        :: config    !! Input collection of user-defined parameter
      real(DP), intent(in)                           :: t         !! Current time. This is passed to the user-defined acceleration function.
      logical, intent(in)                            :: lflag     !! Logical flag indicating whether to recompute direct cross term accelrations
      class(whm_pl), optional, intent(in)          :: whm_plA !! Helio massive body particle data structure. Optional to allow this method to be polymorphic for pl and tp classes
      real(DP), dimension(:,:), optional, intent(in) :: xh        !! Heliocentric positions of massive bodies at time t
   end subroutine whm_getacch_tp

   module subroutine whm_drift_pl(self, dt, mu)
      implicit none
      class(whm_pl), intent(inout) :: self !! Helio massive body particle data structure
      real(DP), intent(in)           :: dt   !! Stepsize
      real(DP), optional, intent(in) :: mu   !! G * (m1 + m2), G = gravitational constant, m1 = mass of central body, m2 = mass of body to drift
   end subroutine whm_drift_pl

   module subroutine whm_drift_tp(self, dt, mu)
      implicit none
      class(whm_tp), intent(inout) :: self !! Helio test particle data structure
      real(DP), intent(in)           :: dt   !! Stepsize
      real(DP), optional, intent(in) :: mu   !! G * m1, G = gravitational constant, m1 = mass of central body
   end subroutine whm_drift_tp

   module subroutine whm_getacch_int_pl(self, whm_plA)
      implicit none
      class(whm_pl), intent(inout)        :: self      !! Helio massive body particle data structure.
      class(whm_pl), optional, intent(in) :: whm_plA !! Dummy argument used to make this a polymorphic method for both pl and tp classes
   end subroutine whm_getacch_int_pl

   module subroutine whm_getacch_int_tp(self, whm_plA)
      implicit none
      class(whm_tp), intent(inout)        :: self      !! Helio test particle data structure
      class(whm_pl), optional, intent(in) :: whm_plA !! Helio massive body particle data structure. Optional to allow this method to be polymorphic for pl and tp classes
   end subroutine whm_getacch_int_tp

   module subroutine whm_kickvb_pl(self, dt)
      implicit none
      class(whm_pl), intent(inout) :: self !! Helio massive body particle data structure.
      real(DP), intent(in)           :: dt   !! Stepsize
   end subroutine whm_kickvb_pl

   module subroutine whm_kickvb_tp(self, dt)
      implicit none
      class(whm_tp), intent(inout) :: self !! Helio test particle data structure
      real(DP), intent(in)           :: dt   !! Stepsize
   end subroutine whm_kickvb_tp

   module subroutine whm_lindrift_pl(self, dt, pt)
      implicit none
      class(whm_pl), intent(inout)           :: self !! Helio massive body particle data structur
      real(DP), intent(in)                     :: dt   !! Stepsize
      real(DP), dimension(NDIM), intent(inout) :: pt   !! Negative barycentric velocity of the Sun
   end subroutine whm_lindrift_pl

   module subroutine whm_lindrift_tp(self, dt, pt)
      implicit none
      class(whm_tp), intent(inout)           :: self !! Helio test particle data structure
      real(DP), intent(in)                     :: dt   !! Stepsize
      real(DP), dimension(NDIM), intent(inout) :: pt   !! Negative barycentric velocity of the Sun
   end subroutine whm_lindrift_tp

   module subroutine whm_step(whm_plA, whm_tpA, config, t, dt, lfirst)
      implicit none
      type(whm_pl), intent(inout)           :: whm_plA !! Helio massive body particle data structure. 
      type(whm_tp), intent(inout)           :: whm_tpA !! Helio test particle data structure
      type(swiftest_configuration),intent(in) :: config    !! Input collection of user-defined parameter
      real(DP), intent(in)                    :: t         !! Current time
      real(DP), intent(in)                    :: dt        !! Stepsize
      logical, intent(inout)                  :: lfirst    !! Logical flag indicating whether current invocation is the first 
                                                           !!     TODO: Replace lfirst with internal flag with save attribute
   end subroutine whm_step

   module subroutine whm_step_pl(self, whm_plA, config, t, dt, lfirst, xbeg, xend, ptb, pte)
      implicit none
      class(whm_pl), intent(inout)           :: self            !! Helio massive body particle data structure.
      class(whm_pl), optional, intent(in)    :: whm_plA       !! Dummy argument used to make this a polymorphic method for both pl and tp classes
      type(swiftest_configuration),intent(in)  :: config          !! Input collection of user-defined parameter
      real(DP), intent(in)                     :: t               !! Current time
      real(DP), intent(in)                     :: dt              !! Stepsize
      logical, intent(inout)                   :: lfirst          !! Logical flag indicating whether current invocation is the first 
                                                                  !!     TODO: Replace lfirst with internal flag with save attribute
      real(DP), dimension(:), intent(inout)    :: ptb             !! Negative barycentric velocity of the Sun at beginning of time step
      real(DP), dimension(:), intent(inout)    :: pte             !! Negative barycentric velocity of the Sun at end of time step
      real(DP), dimension(:, :), intent(inout) :: xbeg            !! Heliocentric massive body positions at beginning of time step
      real(DP), dimension(:, :), intent(inout) :: xend            !! Heliocentric massive body positions at end of time step
   end subroutine whm_step_pl

   module subroutine whm_step_tp(self, whm_plA, config, t, dt, lfirst, xbeg, xend, ptb, pte)
      implicit none
      class(whm_tp), intent(inout)           :: self      !! Helio test particle data structure.
      class(whm_pl), optional, intent(in)    :: whm_plA !! Helio massive body particle data structure.
      type(swiftest_configuration),intent(in)  :: config    !! Input collection of user-defined parameter
      real(DP), intent(in)                     :: t         !! Current time
      real(DP), intent(in)                     :: dt        !! Stepsize
      logical, intent(inout)                   :: lfirst    !! Logical flag indicating whether current invocation is the first 
                                                            !!     TODO: Replace lfirst with internal flag with save attribute
      real(DP), dimension(:), intent(inout)    :: ptb       !! Negative barycentric velocity of the Sun at beginning of time step
      real(DP), dimension(:), intent(inout)    :: pte       !! Negative barycentric velocity of the Sun at end of time step
      real(DP), dimension(:, :), intent(inout) :: xbeg      !! Heliocentric massive body positions at beginning of time step
      real(DP), dimension(:, :), intent(inout) :: xend      !! Heliocentric massive body positions at end of time step
   end subroutine whm_step_tp

   module subroutine whm_user_getacch_pl(whm_plA, t)
      implicit none
      class(whm_pl), intent(inout) :: whm_plA  !! Helio massive body particle data structure
      real(DP), intent(in)           :: t          !! Current time
   end subroutine whm_user_getacch_pl

   module subroutine whm_user_getacch_tp(whm_tpA, t)
      implicit none
      class(whm_tp), intent(inout) :: whm_tpA  !! Helio massive body particle data structure
      real(DP), intent(in)           :: t          !! Current tim
   end subroutine whm_user_getacch_tp
end interface

!> Only the constructor and destructor method implementations are listed here. All other methods are implemented in the whm submodules.
   contains





      

end module whm
