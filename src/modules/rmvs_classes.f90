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
   ! rmvs_cb class definitions and method interfaces
   !*******************************************************************************************************************************
   !> RMVS central body particle class
   type, public, extends(whm_cb) :: rmvs_cb
      real(DP), dimension(NDIM) :: xin = (/0.0_DP, 0.0_DP, 0.0_DP/) !! interpolated heliocentric central body position for inner encounter
      real(DP), dimension(NDIM) :: vin = (/0.0_DP, 0.0_DP, 0.0_DP/) !! interpolated heliocentric central body velocity for inner encounter
   contains
   end type rmvs_cb

   !********************************************************************************************************************************
   !  rmvs_tp class definitions and method interfaces
   !*******************************************************************************************************************************

   !! RMVS test particle class
   type, public, extends(whm_tp) :: rmvs_tp
      !! Note to developers: If you add componenets to this class, be sure to update methods and subroutines that traverse the
      !!    component list, such as rmvs_setup_tp and rmvs_spill_tp
      ! encounter steps)
      logical                                   :: lplanetocentric = .false.   !! Flag that indicates that the object is a planetocentric set of test particles used for close encounter calculations
      logical,      dimension(:),   allocatable :: lperi  !! planetocentric pericenter passage flag (persistent for a full rmvs time step) over a full RMVS time step)
      integer(I4B), dimension(:),   allocatable :: plperP !! index of planet associated with pericenter distance peri (persistent over a full RMVS time step)
      integer(I4B), dimension(:),   allocatable :: plencP !! index of planet that test particle is encountering (not persistent for a full RMVS time step)
      real(DP),     dimension(:,:), allocatable :: vbeg   !! Planet velocities at beginning ot step

      ! The following are used to correctly set the oblateness values of the acceleration during an inner encounter with a planet
      type(rmvs_cb)                                  :: cb     !! Copy of original central body object passed to close encounter (used for oblateness acceleration during planetocentric encoountters)
      real(DP),     dimension(:,:),      allocatable :: xheliocen !! original heliocentric position (used for oblateness calculation during close encounters)
      integer(I4B)                                   :: index   !!  inner substep number within current set
      integer(I4B)                                   :: ipleP   !!  index value of encountering planet
      real(DP),     dimension(:,:),      allocatable :: aoblin_pl !! Encountering planet's oblateness acceleration value
      real(DP),     dimension(:,:),      allocatable :: xh_pl     !! Encountering planet's heliocentric position values
   contains
      procedure, public :: discard_pl        => rmvs_discard_pl_tp
      procedure, public :: encounter_check   => rmvs_encounter_check_tp !! Checks if any test particles are undergoing a close encounter with a massive body
      procedure, public :: fill              => rmvs_fill_tp            !! "Fills" bodies from one object into another depending on the results of a mask (uses the MERGE intrinsic)
      procedure, public :: getacch           => rmvs_getacch_in_tp      !!  Calculates either the standard or modified version of the oblateness acceleration depending if the
                                                                        !! if the test particle is undergoing a close encounter or not
      procedure, public :: peri_pass         => rmvs_peri_tp            !! Determine planetocentric pericenter passages for test particles in close encounters with a planet
      procedure, public :: set_beg_end       => rmvs_setup_set_beg_end  !! Sets the beginning and ending values of planet positions. Also adds the end velocity for RMVS
      procedure, public :: setup             => rmvs_setup_tp           !! Constructor method - Allocates space for number of particles
      procedure, public :: spill             => rmvs_spill_tp           !! "Spills" bodies from one object to another depending on the results of a mask (uses the PACK intrinsic)
   end type rmvs_tp

   !********************************************************************************************************************************
   !                                    rmvs_pl class definitions and method interfaces
   !*******************************************************************************************************************************

   !> RMVS massive body particle class
   type, public, extends(whm_pl) :: rmvs_pl
      integer(I4B),  dimension(:),       allocatable :: nenc    !! number of test particles encountering planet this full rmvs time step
      integer(I4B),  dimension(:),       allocatable :: tpenc1P !! index of first test particle encountering planet
      real(DP),      dimension(:, :, :), allocatable :: xout    !! interpolated heliocentric planet position for outer encounter
      real(DP),      dimension(:, :, :), allocatable :: vout    !! interpolated heliocentric planet velocity for outer encounter
      real(DP),      dimension(:, :, :), allocatable :: xin     !! interpolated heliocentric planet position for inner encounter
      real(DP),      dimension(:, :, :), allocatable :: vin     !! interpolated heliocentric planet velocity for inner encounter
      type(rmvs_tp), dimension(:),       allocatable :: tpenc   !! array of encountering test particles with this planet in planetocentric coordinates
      type(whm_pl) , dimension(:),       allocatable :: plenc   !! array of massive bodies that includes the Sun, but not the encountering planet in planetocentric coordinates
      type(rmvs_cb), dimension(:),       allocatable :: cbenc   !! The planet acting as a central body for close encounters
      real(DP),      dimension(:, :, :), allocatable :: aoblin  !! barycentric acceleration on planets due to central body oblateness during inner encounter
      integer(I4B),  dimension(:,:),     allocatable :: plind
      !! Note to developers: If you add componenets to this class, be sure to update methods and subroutines that traverse the
      !!    component list, such as rmvs_setup_pl and rmvs_spill_pl
   contains
      procedure, public :: fill                => rmvs_fill_pl    !! "Fills" bodies from one object into another depending on the results of a mask (uses the MERGE intrinsic)
      procedure, public :: interp_in           => rmvs_interp_in   !! Interpolate planet positions between two Keplerian orbits in inner encounter region
      procedure, public :: interp_out          => rmvs_interp_out  !! Interpolate planet positions between two Keplerian orbits in outer encounter region
      procedure, public :: obl_acc_in          => rmvs_obl_acc_in  !! Compute the oblateness acceleration in the inner encounter region with planets 
      procedure, public :: setup               => rmvs_setup_pl    !! Constructor method - Allocates space for number of particles
      procedure, public :: step_in             => rmvs_step_in_pl  !! Step active test particles ahead in the inner encounter region with planets   
      procedure, public :: step_out            => rmvs_step_out    !! Step active test particles ahead in the outer encounter region with planets
      procedure, public :: make_planetocentric => rmvs_step_make_planetocentric !! Creates encountering test particle structure for the planets
      procedure, public :: end_planetocentric  => rmvs_step_end_planetocentric  !! Creates encountering test particle structure for the planets
      procedure, public :: spill               => rmvs_spill_pl    !! "Spills" bodies from one object to another depending on the results of a mask (uses the PACK intrinsic)
   end type rmvs_pl

   !********************************************************************************************************************************
   !  rmvs_nbody_system class definitions and method interfaces
   !********************************************************************************************************************************
   type, public, extends(whm_nbody_system) :: rmvs_nbody_system
      !> In the RMVS integrator, only test particles are discarded
   contains
      private
      !> Replace the abstract procedures with concrete ones
      procedure, public :: initialize    => rmvs_setup_system  !! Performs RMVS-specific initilization steps, like calculating the Jacobi masses
   end type rmvs_nbody_system
   
   interface
      module subroutine rmvs_discard_pl_tp(self, pl, t, dt)
         use swiftest_classes, only : swiftest_cb, swiftest_pl, swiftest_configuration
         implicit none
         class(rmvs_tp),                intent(inout) :: self
         class(swiftest_pl),            intent(inout) :: pl     !! WHM massive body particle data structure. 
         real(DP),                      intent(in)    :: t      !! Current simulation time
         real(DP),                      intent(in)    :: dt     !! Stepsize
      end subroutine rmvs_discard_pl_tp

      module subroutine rmvs_getacch_in_tp(self, cb, pl, config, t, xh)
         use swiftest_classes, only : swiftest_cb, swiftest_configuration
         use whm_classes, only : whm_pl
         implicit none
         class(rmvs_tp),                intent(inout) :: self   !! RMVS test particle data structure
         class(swiftest_cb),            intent(inout) :: cb     !! Swiftest central body particle data structuree 
         class(whm_pl),                 intent(inout) :: pl     !! WHM massive body particle data structure. 
         class(swiftest_configuration), intent(in)    :: config !! Input collection of  parameter
         real(DP),                      intent(in)    :: t      !! Current time
         real(DP), dimension(:,:),      intent(in)    :: xh     !! Heliocentric positions of planets
      end subroutine rmvs_getacch_in_tp

      module subroutine rmvs_obl_acc_in(self, cb)
         use swiftest_classes, only : swiftest_cb
         implicit none
         class(rmvs_pl),         intent(inout) :: self !! RMVS massive body object
         class(swiftest_cb),     intent(inout) :: cb   !! Swiftest central body object
      end subroutine rmvs_obl_acc_in

      module subroutine rmvs_setup_system(self, config)
         use swiftest_classes, only : swiftest_configuration
         implicit none
         class(rmvs_nbody_system),      intent(inout) :: self    !! RMVS system object
         class(swiftest_configuration), intent(inout) :: config  !! Input collection of  configuration parameters 
      end subroutine rmvs_setup_system

      module subroutine rmvs_setup_pl(self,n)
         implicit none
         class(rmvs_pl),  intent(inout) :: self !! RMVS test particle object
         integer,         intent(in)    :: n    !! Number of test particles to allocate
      end subroutine rmvs_setup_pl

      module subroutine rmvs_step_make_planetocentric(self, cb, tp, config)
         use swiftest_classes, only : swiftest_configuration
         implicit none
         class(rmvs_pl),   intent(inout)  :: self !! RMVS test particle object
         class(rmvs_cb),   intent(inout)  :: cb   !! RMVS central body particle type
         class(rmvs_tp),   intent(inout)  :: tp   !! RMVS test particle object
         class(swiftest_configuration), intent(in) :: config !! Input collection of configuration parameters 
      end subroutine rmvs_step_make_planetocentric

      module subroutine rmvs_step_end_planetocentric(self, cb, tp)
         implicit none
         class(rmvs_pl),   intent(inout)  :: self !! RMVS test particle object
         class(rmvs_cb),   intent(inout)  :: cb   !!  RMVS central body object
         class(rmvs_tp),   intent(inout)  :: tp   !! RMVS test particle object
      end subroutine rmvs_step_end_planetocentric

      module subroutine rmvs_setup_set_beg_end(self, xbeg, xend, vbeg)
         implicit none
         class(rmvs_tp),   intent(inout)  :: self !! RMVS test particle object
         real(DP), dimension(:,:), optional :: xbeg, xend, vbeg
      end subroutine rmvs_setup_set_beg_end

      module subroutine rmvs_destruct_planetocentric_encounter(self)
         implicit none
         class(rmvs_pl),   intent(inout)  :: self !! RMVS test particle object
      end subroutine rmvs_destruct_planetocentric_encounter

      module subroutine rmvs_interp_in(self, cb, dt)
         implicit none
         class(rmvs_pl), intent(inout)   :: self !! RMVS test particle object
         class(rmvs_cb), intent(in)      :: cb   !! RMVS central body particle type
         real(DP), intent(in)            :: dt   !! Step size
      end subroutine rmvs_interp_in

      module subroutine rmvs_interp_out(self, cb, dt)
         implicit none
         class(rmvs_pl), intent(inout)   :: self !! RMVS test particle object
         class(rmvs_cb), intent(in)      :: cb   !! RMVS central body particle type
         real(DP), intent(in)            :: dt   !! Step size
      end subroutine rmvs_interp_out

      module subroutine rmvs_step_in_pl(self, cb, tp, config, t, dt)
         use swiftest_classes, only : swiftest_configuration
         implicit none 
         class(rmvs_pl),                intent(inout)  :: self !! RMVS massive body object
         class(rmvs_cb),                intent(inout)  :: cb   !! RMVS central body object
         class(rmvs_tp),                intent(inout)  :: tp   !! RMVS test particle object
         class(swiftest_configuration), intent(in)     :: config  !! Input collection of configuration parameters 
         real(DP),                      intent(in)     :: t    !! Current time
         real(DP),                      intent(in)     :: dt   !! Step size
      end subroutine rmvs_step_in_pl

      module subroutine rmvs_step_out(self, cb, tp, dt, config)
         use swiftest_classes, only : swiftest_configuration
         implicit none
         class(rmvs_pl),                intent(inout)  :: self !! RMVS massive body object
         class(rmvs_cb),                intent(inout)  :: cb   !! RMVS central body object
         class(rmvs_tp),                intent(inout)  :: tp   !! RMVS test particle object
         real(DP),                      intent(in)     :: dt   !! Step size
         class(swiftest_configuration), intent(in)     :: config  !! Input collection of  configuration parameters 
      end subroutine rmvs_step_out

      module subroutine rmvs_step_out2(cb, pl, tp, t, dt, index, config)
         use swiftest_classes, only : swiftest_configuration
         implicit none
         class(rmvs_cb),                intent(inout) :: cb   !! RMVS central body object
         class(rmvs_pl),                intent(inout)  :: pl   !! RMVS massive body object
         class(rmvs_tp),                intent(inout) :: tp   !! RMVS test particle object
         real(DP),                      intent(in)    :: t     !! Simulation time
         real(DP),                      intent(in)    :: dt    !! Step size
         integer(I4B),                  intent(in)    :: index !! outer substep number within current set
         class(swiftest_configuration), intent(in)    :: config  !! Input collection of  configuration parameters 
      end subroutine rmvs_step_out2   

      module elemental function rmvs_chk_ind(r2, v2, vdotr, dt, r2crit) result(lflag)
         implicit none
         real(DP), intent(in)                :: r2, v2, vdotr, dt, r2crit
         logical                             :: lflag
      end function rmvs_chk_ind

      module subroutine rmvs_step_system(cb, pl, tp, config)
         use swiftest_classes, only : swiftest_configuration
         implicit none
         class(rmvs_cb),                intent(inout) :: cb      !! RMVS central body object  
         class(rmvs_pl),                intent(inout) :: pl      !! RMVS central body object  
         class(rmvs_tp),                intent(inout) :: tp      !! RMVS central body object  
         class(swiftest_configuration), intent(in)    :: config  !! Input collection of  configuration parameters 
      end subroutine rmvs_step_system

      module subroutine rmvs_setup_tp(self,n)
         implicit none
         class(rmvs_tp), intent(inout)   :: self !! RMVS test particle object
         integer,        intent(in)      :: n    !! Number of test particles to allocate
      end subroutine rmvs_setup_tp

      module function rmvs_encounter_check_tp(self, cb, pl, dt, rts) result(lencounter)
         implicit none
         class(rmvs_tp),            intent(inout) :: self        !! RMVS test particle object  
         class(rmvs_cb),            intent(inout) :: cb          !! RMVS central body object  
         class(rmvs_pl),            intent(inout) :: pl          !! RMVS massive body object  
         real(DP),                  intent(in)    :: dt          !! step size
         real(DP),                  intent(in)    :: rts         !! fraction of Hill's sphere radius to use as radius of encounter regio
         logical                                  :: lencounter  !! Returns true if there is at least one close encounter      

      end function rmvs_encounter_check_tp

      module subroutine rmvs_peri_tp(self, cb, pl, t, dt, lfirst, index, ipleP, config)
         use swiftest_classes, only : swiftest_configuration
         class(rmvs_tp),                intent(inout) :: self   !! RMVS test particle object  
         class(rmvs_cb),                intent(inout) :: cb     !! RMVS central body object  
         class(rmvs_pl),                intent(inout) :: pl     !! RMVS massive body object  
         real(DP),                      intent(in)    :: t      !! current time
         real(DP),                      intent(in)    :: dt     !! step size
         logical,                       intent(in)    :: lfirst !! Logical flag indicating whether current invocation is the first
         integer(I4B),                  intent(in)    :: index !! outer substep number within current set
         integer(I4B),                  intent(in)    :: ipleP !!index of RMVS planet being closely encountered
         class(swiftest_configuration), intent(in)    :: config  !! Input collection of  configuration parameters 
      end subroutine rmvs_peri_tp

      module subroutine rmvs_spill_pl(self, discards, lspill_list)
         use swiftest_classes, only : swiftest_body
         implicit none
         class(rmvs_pl), intent(inout) :: self      !! Swiftest massive body body object
         class(swiftest_body), intent(inout) :: discards    !! Discarded object 
         logical, dimension(:), intent(in)   :: lspill_list !! Logical array of bodies to spill into the discards
      end subroutine rmvs_spill_pl

      module subroutine rmvs_fill_pl(self, inserts, lfill_list)
         use swiftest_classes, only : swiftest_body
         implicit none
         class(rmvs_pl),        intent(inout) :: self       !! RMVS massive body object
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
