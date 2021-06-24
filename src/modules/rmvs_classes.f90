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
      logical                                  :: lplanetocentric = .false.   !! Flag that indicates that the object is a planetocentric set of masive bodies used for close encounter calculations
   contains
      private
      !> Replace the abstract procedures with concrete ones
      procedure, public :: initialize    => rmvs_setup_system  !! Performs RMVS-specific initilization steps, like calculating the Jacobi masses
      procedure, public :: step          => rmvs_step_system
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
      logical                                   :: lplanetocentric = .false.  !! Flag that indicates that the object is a planetocentric set of test particles used for close encounter calculations
      type(rmvs_interp),           dimension(:), allocatable :: outer !! interpolated heliocentric central body position for outer encounters
      type(rmvs_interp),           dimension(:), allocatable :: inner !! interpolated heliocentric central body position for inner encounters
   end type rmvs_cb

   !********************************************************************************************************************************
   !  rmvs_tp class definitions and method interfaces
   !*******************************************************************************************************************************

   !! RMVS test particle class
   type, public, extends(whm_tp) :: rmvs_tp
      !! Note to developers: If you add componenets to this class, be sure to update methods and subroutines that traverse the
      !!    component list, such as rmvs_setup_tp and rmvs_spill_tp
      ! encounter steps)
      logical                                   :: lplanetocentric = .false.  !! Flag that indicates that the object is a planetocentric set of test particles used for close encounter calculations
      logical,      dimension(:),   allocatable :: lperi  !! planetocentric pericenter passage flag (persistent for a full rmvs time step) over a full RMVS time step)
      integer(I4B), dimension(:),   allocatable :: plperP !! index of planet associated with pericenter distance peri (persistent over a full RMVS time step)
      integer(I4B), dimension(:),   allocatable :: plencP !! index of planet that test particle is encountering (not persistent for a full RMVS time step)
      real(DP),     dimension(:,:), allocatable :: vbeg   !! Planet velocities at beginning ot step

      ! The following are used to correctly set the oblateness values of the acceleration during an inner encounter with a planet
      type(rmvs_cb)                             :: cb_heliocentric     !! Copy of original central body object passed to close encounter (used for oblateness acceleration during planetocentric encoountters)
      real(DP),     dimension(:,:), allocatable :: xheliocentric !! original heliocentric position (used for oblateness calculation during close encounters)
      integer(I4B)                              :: index   !!  inner substep number within current set
      integer(I4B)                              :: ipleP   !!  index value of encountering planet
   contains
      procedure, public :: discard_pl        => rmvs_discard_pl_tp
      procedure, public :: encounter_check   => rmvs_encounter_check_tp !! Checks if any test particles are undergoing a close encounter with a massive body
      procedure, public :: fill              => rmvs_fill_tp            !! "Fills" bodies from one object into another depending on the results of a mask (uses the MERGE intrinsic)
      procedure, public :: getacch           => rmvs_getacch_tp         !!  Calculates either the standard or modified version of the acceleration depending if the
                                                                        !! if the test particle is undergoing a close encounter or not
      procedure, public :: set_beg_end       => rmvs_setup_set_beg_end  !! Sets the beginning and ending values of planet positions. Also adds the end velocity for RMVS
      procedure, public :: setup             => rmvs_setup_tp           !! Constructor method - Allocates space for number of particles
      procedure, public :: spill             => rmvs_spill_tp           !! "Spills" bodies from one object to another depending on the results of a mask (uses the PACK intrinsic)
   end type rmvs_tp

   !********************************************************************************************************************************
   !                                    rmvs_pl class definitions and method interfaces
   !*******************************************************************************************************************************
 

   !> RMVS massive body particle class
   type, private, extends(whm_pl) :: rmvs_base_pl
      logical                                                :: lplanetocentric = .false.   !! Flag that indicates that the object is a planetocentric set of masive bodies used for close encounter calculations
      integer(I4B),                dimension(:), allocatable :: nenc    !! number of test particles encountering planet this full rmvs time step
      integer(I4B),                dimension(:), allocatable :: tpenc1P !! index of first test particle encountering planet
      integer(I4B),                dimension(:), allocatable :: plind ! Connects the planetocentric indices back to the heliocentric planet list
      type(rmvs_interp),           dimension(:), allocatable :: outer !! interpolated heliocentric central body position for outer encounters
      type(rmvs_interp),           dimension(:), allocatable :: inner !! interpolated heliocentric central body position for inner encounters
   contains
      procedure, public :: fill                => rmvs_fill_pl    !! "Fills" bodies from one object into another depending on the results of a mask (uses the MERGE intrinsic)
      procedure, public :: setup               => rmvs_setup_pl    !! Constructor method - Allocates space for number of particles
      procedure, public :: spill               => rmvs_spill_pl    !! "Spills" bodies from one object to another depending on the results of a mask (uses the PACK intrinsic)
   end type rmvs_base_pl

   type, public :: rmvs_encounter_system
      class(rmvs_cb),      allocatable :: cb
      class(rmvs_tp),      allocatable :: tp
      class(rmvs_base_pl), allocatable :: pl
   end type rmvs_encounter_system

   type, public, extends(rmvs_base_pl) :: rmvs_pl
      type(rmvs_encounter_system), dimension(:), allocatable :: planetocentric
      !! Note to developers: If you add componenets to this class, be sure to update methods and subroutines that traverse the
      !!    component list, such as rmvs_setup_pl and rmvs_spill_pl
   end type rmvs_pl
   
   interface
      module subroutine rmvs_discard_pl_tp(self, pl, t, dt)
         use swiftest_classes, only : swiftest_cb, swiftest_pl, swiftest_parameters
         implicit none
         class(rmvs_tp),                intent(inout) :: self
         class(swiftest_pl),            intent(inout) :: pl     !! WHM massive body particle data structure. 
         real(DP),                      intent(in)    :: t      !! Current simulation time
         real(DP),                      intent(in)    :: dt     !! Stepsize
      end subroutine rmvs_discard_pl_tp

      module subroutine rmvs_getacch_tp(self, cb, pl, param, t, xh)
         use swiftest_classes, only : swiftest_cb, swiftest_parameters
         use whm_classes, only : whm_pl
         implicit none
         class(rmvs_tp),                intent(inout) :: self   !! RMVS test particle data structure
         class(swiftest_cb),            intent(inout) :: cb     !! Swiftest central body particle data structuree 
         class(whm_pl),                 intent(inout) :: pl     !! WHM massive body particle data structure. 
         class(swiftest_parameters), intent(in)    :: param !! Input collection of  parameter
         real(DP),                      intent(in)    :: t      !! Current time
         real(DP), dimension(:,:),      intent(in)    :: xh     !! Heliocentric positions of planets
      end subroutine rmvs_getacch_tp

      module subroutine rmvs_setup_system(self, param)
         use swiftest_classes, only : swiftest_parameters
         implicit none
         class(rmvs_nbody_system),      intent(inout) :: self    !! RMVS system object
         class(swiftest_parameters), intent(inout) :: param  !! Input collection of  parameters parameters 
      end subroutine rmvs_setup_system

      module subroutine rmvs_setup_pl(self,n)
         implicit none
         class(rmvs_base_pl), intent(inout) :: self !! RMVS test particle object
         integer,             intent(in)    :: n    !! Number of test particles to allocate
      end subroutine rmvs_setup_pl

      module subroutine rmvs_setup_set_beg_end(self, xbeg, xend, vbeg)
         implicit none
         class(rmvs_tp),           intent(inout)          :: self !! RMVS test particle object
         real(DP), dimension(:,:), intent(in),   optional :: xbeg, xend, vbeg
      end subroutine rmvs_setup_set_beg_end

      module subroutine rmvs_step_system(self, param)
         use swiftest_classes, only : swiftest_parameters
         implicit none
         class(rmvs_nbody_system),      intent(inout) :: self    !! RMVS nbody system object
         class(swiftest_parameters), intent(in)    :: param  !! Input collection of  parameters parameters 
      end subroutine rmvs_step_system

      module subroutine rmvs_setup_tp(self,n)
         implicit none
         class(rmvs_tp), intent(inout)   :: self !! RMVS test particle object
         integer,        intent(in)      :: n    !! Number of test particles to allocate
      end subroutine rmvs_setup_tp

      module function rmvs_encounter_check_tp(self, cb, pl, dt, rts) result(lencounter)
         implicit none
         class(rmvs_tp),       intent(inout) :: self        !! RMVS test particle object  
         class(rmvs_cb),       intent(inout) :: cb          !! RMVS central body object  
         class(rmvs_pl),       intent(inout) :: pl          !! RMVS massive body object  
         real(DP),             intent(in)    :: dt          !! step size
         real(DP),             intent(in)    :: rts         !! fraction of Hill's sphere radius to use as radius of encounter regio
         logical                             :: lencounter  !! Returns true if there is at least one close encounter      
      end function rmvs_encounter_check_tp

      module subroutine rmvs_spill_pl(self, discards, lspill_list)
         use swiftest_classes, only : swiftest_body
         implicit none
         class(rmvs_base_pl),   intent(inout) :: self      !! RMVS massive body object
         class(swiftest_body),  intent(inout) :: discards    !! Discarded object 
         logical, dimension(:), intent(in)   :: lspill_list !! Logical array of bodies to spill into the discards
      end subroutine rmvs_spill_pl

      module subroutine rmvs_fill_pl(self, inserts, lfill_list)
         use swiftest_classes, only : swiftest_body
         implicit none
         class(rmvs_base_pl),   intent(inout) :: self       !! RMVS massive body object 
         class(swiftest_body),  intent(inout) :: inserts    !! Inserted object 
         logical, dimension(:), intent(in)    :: lfill_list !! Logical array of bodies to merge into the keeps
      end subroutine rmvs_fill_pl

      module subroutine rmvs_spill_tp(self, discards, lspill_list)
         use swiftest_classes, only : swiftest_body
         implicit none
         class(rmvs_tp),        intent(inout) :: self       !! RMVS test particle object
         class(swiftest_body),  intent(inout) :: discards    !! Discarded object 
         logical, dimension(:), intent(in)    :: lspill_list !! Logical array of bodies to spill into the discards
      end subroutine rmvs_spill_tp

      module subroutine rmvs_fill_tp(self, inserts, lfill_list)
         use swiftest_classes, only : swiftest_body
         implicit none
         class(rmvs_tp),        intent(inout) :: self        !! RMVS massive body object
         class(swiftest_body),  intent(inout) :: inserts     !!  Inserted object 
         logical, dimension(:), intent(in)    :: lfill_list !! Logical array of bodies to merge into the keeps
      end subroutine rmvs_fill_tp
   end interface

end module rmvs_classes
