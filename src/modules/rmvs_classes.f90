module rmvs_classes
   !! author: David A. Minton
   !!
   !! Definition of classes and methods specific to the Regularized Mixed Variable Symplectic (RMVS) integrator
   !! Partially adapted from David E. Kaufmann's Swifter module: module_rmvs.f90
   use swiftest_globals
   use whm_classes, only : whm_cb, whm_pl, whm_tp, whm_nbody_system
   implicit none

   public

   integer(I4B), parameter :: NTENC = 10
   integer(I4B), parameter :: NTPHENC = 3
   integer(I4B), parameter :: NTPENC = NTENC * NTPHENC
   real(DP), parameter     :: RHSCALE = 3.5_DP
   real(DP), parameter     :: RHPSCALE = 1.0_DP
   real(DP), parameter     :: FACQDT = 2.0_DP

   !********************************************************************************************************************************
   !  rmvs_nbody_system class definitions and method interfaces
   !********************************************************************************************************************************
   type, public, extends(whm_nbody_system) :: rmvs_nbody_system
      !> In the RMVS integrator, only test particles are discarded
      logical                               :: lplanetocentric = .false.  !! Flag that indicates that the object is a planetocentric set of masive bodies used for close encounter calculations
      real(DP)                              :: rts                        !! fraction of Hill's sphere radius to use as radius of encounter region
      real(DP), dimension(:,:), allocatable :: vbeg                       !! Planet velocities at beginning ot step
   contains
      private
      !> Replace the abstract procedures with concrete ones
      procedure, public :: initialize    => rmvs_setup_system  !! Performs RMVS-specific initilization steps, like calculating the Jacobi masses
      procedure, public :: step          => rmvs_step_system
      procedure, public :: set_beg_end   => rmvs_setup_set_beg_end  !! Sets the beginning and ending values of planet positions. Also adds the end velocity for RMVS
   end type rmvs_nbody_system

   type, private :: rmvs_interp
      real(DP), dimension(:, :), allocatable :: x    !! interpolated heliocentric planet position for outer encounter
      real(DP), dimension(:, :), allocatable :: v    !! interpolated heliocentric planet velocity for outer encounter
      real(DP), dimension(:, :), allocatable :: aobl !! Encountering planet's oblateness acceleration value
   end type rmvs_interp

   !********************************************************************************************************************************
   ! rmvs_cb class definitions and method interfaces
   !*******************************************************************************************************************************
   !> RMVS central body particle class 
   type, public, extends(whm_cb) :: rmvs_cb
      type(rmvs_interp), dimension(:), allocatable :: outer !! interpolated heliocentric central body position for outer encounters
      type(rmvs_interp), dimension(:), allocatable :: inner !! interpolated heliocentric central body position for inner encounters
      logical                                      :: lplanetocentric = .false.  !! Flag that indicates that the object is a planetocentric set of masive bodies used for close encounter calculations
   end type rmvs_cb

   !********************************************************************************************************************************
   !  rmvs_tp class definitions and method interfaces
   !*******************************************************************************************************************************

   !! RMVS test particle class
   type, public, extends(whm_tp) :: rmvs_tp
      !! Note to developers: If you add componenets to this class, be sure to update methods and subroutines that traverse the
      !!    component list, such as rmvs_setup_tp and rmvs_spill_tp
      ! encounter steps)
      logical,      dimension(:),   allocatable :: lperi  !! planetocentric pericenter passage flag (persistent for a full rmvs time step) over a full RMVS time step)
      integer(I4B), dimension(:),   allocatable :: plperP !! index of planet associated with pericenter distance peri (persistent over a full RMVS time step)
      integer(I4B), dimension(:),   allocatable :: plencP !! index of planet that test particle is encountering (not persistent for a full RMVS time step)

      ! The following are used to correctly set the oblateness values of the acceleration during an inner encounter with a planet
      type(rmvs_cb)                             :: cb_heliocentric !! Copy of original central body object passed to close encounter (used for oblateness acceleration during planetocentric encoountters)
      real(DP),     dimension(:,:), allocatable :: xheliocentric   !! original heliocentric position (used for oblateness calculation during close encounters)
      integer(I4B)                              :: index           !!  inner substep number within current set
      integer(I4B)                              :: ipleP           !!  index value of encountering planet
      logical                                   :: lplanetocentric = .false.  !! Flag that indicates that the object is a planetocentric set of masive bodies used for close encounter calculations
   contains
      procedure, public :: discard           => rmvs_discard_tp         !! Check to see if test particles should be discarded based on pericenter passage distances with respect to planets encountered
      procedure, public :: encounter_check   => rmvs_encounter_check_tp !! Checks if any test particles are undergoing a close encounter with a massive body
      procedure, public :: fill              => rmvs_fill_tp            !! "Fills" bodies from one object into another depending on the results of a mask (uses the MERGE intrinsic)
      procedure, public :: getacch           => rmvs_getacch_tp         !!  Calculates either the standard or modified version of the acceleration depending if the
                                                                        !! if the test particle is undergoing a close encounter or not
      procedure, public :: setup             => rmvs_setup_tp           !! Constructor method - Allocates space for number of particles
      procedure, public :: spill             => rmvs_spill_tp           !! "Spills" bodies from one object to another depending on the results of a mask (uses the PACK intrinsic)
   end type rmvs_tp

   !********************************************************************************************************************************
   !                                    rmvs_pl class definitions and method interfaces
   !*******************************************************************************************************************************
 
   !> RMVS massive body particle class
   type, private, extends(whm_pl) :: rmvs_pl
      integer(I4B),            dimension(:), allocatable :: nenc    !! number of test particles encountering planet this full rmvs time step
      integer(I4B),            dimension(:), allocatable :: tpenc1P !! index of first test particle encountering planet
      integer(I4B),            dimension(:), allocatable :: plind ! Connects the planetocentric indices back to the heliocentric planet list
      type(rmvs_interp),       dimension(:), allocatable :: outer !! interpolated heliocentric central body position for outer encounters
      type(rmvs_interp),       dimension(:), allocatable :: inner !! interpolated heliocentric central body position for inner encounters
      class(rmvs_nbody_system), dimension(:), allocatable :: planetocentric
      logical                                            :: lplanetocentric = .false.  !! Flag that indicates that the object is a planetocentric set of masive bodies used for close encounter calculations
   contains
      procedure, public :: fill                => rmvs_fill_pl    !! "Fills" bodies from one object into another depending on the results of a mask (uses the MERGE intrinsic)
      procedure, public :: setup               => rmvs_setup_pl    !! Constructor method - Allocates space for number of particles
      procedure, public :: spill               => rmvs_spill_pl    !! "Spills" bodies from one object to another depending on the results of a mask (uses the PACK intrinsic)
   end type rmvs_pl

   interface
      module subroutine rmvs_discard_tp(self, system, param)
         use swiftest_classes, only : swiftest_nbody_system, swiftest_parameters
         implicit none
         class(rmvs_tp),               intent(inout) :: self   !! RMVS test particle object
         class(swiftest_nbody_system), intent(inout) :: system !! Swiftest nbody system object
         class(swiftest_parameters),   intent(in)    :: param  !! Current run configuration parameters 
      end subroutine rmvs_discard_tp

      module function rmvs_encounter_check_tp(self, system, dt) result(lencounter)
         implicit none
         class(rmvs_tp),           intent(inout) :: self       !! RMVS test particle object  
         class(rmvs_nbody_system), intent(inout) :: system     !! RMVS nbody system object
         real(DP),                 intent(in)    :: dt         !! step size
         logical                                 :: lencounter !! Returns true if there is at least one close encounter      
      end function rmvs_encounter_check_tp

      module subroutine rmvs_fill_pl(self, inserts, lfill_list)
         use swiftest_classes, only : swiftest_body
         implicit none
         class(rmvs_pl),        intent(inout) :: self       !! RMVS massive body object 
         class(swiftest_body),  intent(inout) :: inserts    !! Inserted object 
         logical, dimension(:), intent(in)    :: lfill_list !! Logical array of bodies to merge into the keeps
      end subroutine rmvs_fill_pl

      module subroutine rmvs_fill_tp(self, inserts, lfill_list)
         use swiftest_classes, only : swiftest_body
         implicit none
         class(rmvs_tp),        intent(inout) :: self        !! RMVS massive body object
         class(swiftest_body),  intent(inout) :: inserts     !!  Inserted object 
         logical, dimension(:), intent(in)    :: lfill_list !! Logical array of bodies to merge into the keeps
      end subroutine rmvs_fill_tp

      module subroutine rmvs_getacch_tp(self, system, param, t, xhp)
         use swiftest_classes, only : swiftest_nbody_system, swiftest_parameters
         implicit none
         class(rmvs_tp),               intent(inout) :: self   !! RMVS test particle data structure
         class(swiftest_nbody_system), intent(inout) :: system !! Swiftest central body particle data structuree 
         class(swiftest_parameters),   intent(in)    :: param  !! Current run configuration parameters
         real(DP),                     intent(in)    :: t      !! Current time
         real(DP), dimension(:,:),     intent(in)    :: xhp    !! Heliocentric positions of planets at current substep
      end subroutine rmvs_getacch_tp

      module subroutine rmvs_setup_pl(self,n)
         implicit none
         class(rmvs_pl), intent(inout) :: self !! RMVS test particle object
         integer,             intent(in)    :: n    !! Number of test particles to allocate
      end subroutine rmvs_setup_pl

      module subroutine rmvs_setup_set_beg_end(self, xbeg, xend, vbeg)
         implicit none
         class(rmvs_nbody_system), intent(inout)          :: self !! RMVS nbody system object
         real(DP), dimension(:,:), intent(in),   optional :: xbeg, xend, vbeg
      end subroutine rmvs_setup_set_beg_end

      module subroutine rmvs_setup_system(self, param)
         use swiftest_classes, only : swiftest_parameters
         implicit none
         class(rmvs_nbody_system),   intent(inout) :: self    !! RMVS system object
         class(swiftest_parameters), intent(inout) :: param  !! Current run configuration parameters 
      end subroutine rmvs_setup_system

      module subroutine rmvs_setup_tp(self,n)
         implicit none
         class(rmvs_tp), intent(inout)   :: self !! RMVS test particle object
         integer,        intent(in)      :: n    !! Number of test particles to allocate
      end subroutine rmvs_setup_tp

      module subroutine rmvs_spill_pl(self, discards, lspill_list)
         use swiftest_classes, only : swiftest_body
         implicit none
         class(rmvs_pl),   intent(inout) :: self      !! RMVS massive body object
         class(swiftest_body),  intent(inout) :: discards    !! Discarded object 
         logical, dimension(:), intent(in)    :: lspill_list !! Logical array of bodies to spill into the discards
      end subroutine rmvs_spill_pl

      module subroutine rmvs_spill_tp(self, discards, lspill_list)
         use swiftest_classes, only : swiftest_body
         implicit none
         class(rmvs_tp),        intent(inout) :: self       !! RMVS test particle object
         class(swiftest_body),  intent(inout) :: discards    !! Discarded object 
         logical, dimension(:), intent(in)    :: lspill_list !! Logical array of bodies to spill into the discards
      end subroutine rmvs_spill_tp

      module subroutine rmvs_step_system(self, param, t, dt)
         use swiftest_classes, only : swiftest_parameters
         implicit none
         class(rmvs_nbody_system),   intent(inout) :: self    !! RMVS nbody system object
         class(swiftest_parameters), intent(inout) :: param  !! Current run configuration parameters 
         real(DP),                   intent(in)    :: t      !! Simulation time
         real(DP),                   intent(in)    :: dt     !! Current stepsize
      end subroutine rmvs_step_system

   end interface

end module rmvs_classes
