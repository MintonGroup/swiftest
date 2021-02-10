module rmvs_classes
   !! author: David A. Minton
   !!
   !! Definition of classes and methods specific to the Democratic Heliocentric Method
   !! Partially adapted from David E. Kaufmann's Swifter module: module_whm.f90
   use swiftest_globals
   use swiftest_classes
   use whm_classes
   implicit none

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
   contains
   end type rmvs_cb

   !********************************************************************************************************************************
   !                                    rmvs_pl class definitions and method interfaces
   !*******************************************************************************************************************************

   !> RMVS massive body particle class
   type, public, extends(whm_pl) :: rmvs_pl
      integer(I4B), dimension(:), allocatable :: nenc    ! number of test particles encountering planet this full rmvs time step
      real(DP), dimension(:, :, :),  allocatable :: xout    ! interpolated heliocentric planet position for outer encounter
      real(DP), dimension(:, :, :),  allocatable :: vout    ! interpolated heliocentric planet velocity for outer encounter
      real(DP), dimension(:, :, :),  allocatable :: xin     ! interpolated heliocentric planet position for inner encounter
      real(DP), dimension(:, :, :),  allocatable :: vin     ! interpolated heliocentric planet velocity for inner encounter
      real(DP), dimension(:, :, :),  allocatable :: xpc     ! interpolated planetocentric planet position for inner encounter
      real(DP), dimension(:, :, :),  allocatable :: aoblin  ! barycentric acceleration due to central body oblateness during
                                                            ! inner encounter
      integer(I4B), dimension(:), allocatable :: tpenc1P ! index of first test particle encountering planet 
      !! Note to developers: If you add componenets to this class, be sure to update methods and subroutines that traverse the
      !!    component list, such as rmvs_setup_pl and rmvs_discard_spill_pl
      contains
      procedure, public :: setup        => rmvs_setup_pl    !! Constructor method - Allocates space for number of particles
      procedure, public :: interp_in    => rmvs_interp_in   !! Interpolate planet positions between two Keplerian orbits in inner encounter region
      procedure, public :: interp_out   => rmvs_interp_out  !! Interpolate planet positions between two Keplerian orbits in outer encounter region
      procedure, public :: step_in      => rmvs_step_in     !! Step active test particles ahead in the inner encounter region
      procedure, public :: step_out     => rmvs_step_out    !! Step active test particles ahead in the outer encounter region
   end type rmvs_pl



   !********************************************************************************************************************************
   !  rmvs_tp class definitions and method interfaces
   !*******************************************************************************************************************************

   !! RMVS test particle class
   type, public, extends(whm_tp) :: rmvs_tp
      !! Note to developers: If you add componenets to this class, be sure to update methods and subroutines that traverse the
      !!    component list, such as rmvs_setup_tp and rmvs_discard_spill_tp
      ! encounter steps)
      logical,      dimension(:),   allocatable :: lperi  ! planetocentric pericenter passage flag (persistent for a full rmvs time step)
         ! over a full RMVS time step)
      real(DP),     dimension(:,:), allocatable :: xpc    ! planetocentric position
      real(DP),     dimension(:,:), allocatable :: vpc    ! planetocentric velocity
      real(DP),     dimension(:,:), allocatable :: apc    ! total planetocentric acceleration
      integer(I4B), dimension(:),   allocatable :: plperP ! index of planet associated with pericenter distance peri (persistent over a
      ! full RMVS time step)
      integer(I4B), dimension(:),   allocatable :: plencP ! index of planet that test particle is encountering (not persistent for a full
      ! RMVS time step)
      integer(I4B), dimension(:),   allocatable :: tpencP ! index of next test particle encountering planet

      real(DP),     dimension(:,:), allocatable :: vbeg   ! Planet velocities at beginning ot step
   contains
      procedure, public :: setup             => rmvs_setup_tp    !! Constructor method - Allocates space for number of particles
      procedure, public :: encounter_check   => rmvs_encounter_check   !! Checks if any test particles are undergoing a close encounter with a massive body
   end type rmvs_tp

   interface
      module subroutine rmvs_setup_tp(self,n)
         implicit none
         class(rmvs_tp), intent(inout)   :: self !! Swiftest test particle object
         integer,        intent(in)      :: n    !! Number of test particles to allocate
      end subroutine rmvs_setup_tp

      module function rmvs_encounter_check(self, cb, pl, dt, rts) result(lencounter)
         implicit none
         class(rmvs_tp),            intent(inout) :: self    !! RMVS test particle object  
         class(rmvs_cb),            intent(inout) :: cb      !! RMVS central body object  
         class(rmvs_pl),            intent(inout) :: pl      !! RMVS massive body object  
         real(DP),                  intent(in)    :: dt  !! time step
         real(DP),                  intent(in)    :: rts !! fraction of Hill's sphere radius to use as radius of encounter regio
         logical                                  :: lencounter        

      end function rmvs_encounter_check
   end interface
   !********************************************************************************************************************************
   !  rmvs_nbody_system class definitions and method interfaces
   !********************************************************************************************************************************
   !> An abstract class for the RMVS integrator nbody system 
   type, public, extends(whm_nbody_system) :: rmvs_nbody_system
      !> In the RMVS integrator, only test particles are discarded
   contains
      private
      !> Replace the abstract procedures with concrete ones
      procedure, public :: initialize    => rmvs_setup_system  !! Performs RMVS-specific initilization steps, like calculating the Jacobi masses
   end type rmvs_nbody_system
   

   interface
      module subroutine rmvs_setup_system(self, config)
         implicit none
         class(rmvs_nbody_system),      intent(inout) :: self    !! Swiftest system object
         class(swiftest_configuration), intent(inout) :: config  !! Input collection of user-defined configuration parameters 
      end subroutine rmvs_setup_system

   end interface

   interface
   module subroutine rmvs_setup_pl(self,n)
      implicit none
      class(rmvs_pl), intent(inout)    :: self !! Swiftest test particle object
      integer, intent(in)              :: n    !! Number of test particles to allocate
   end subroutine rmvs_setup_pl

   module subroutine rmvs_interp_in(self, cb, dt)
      implicit none
      class(rmvs_pl), intent(inout)   :: self !! Swiftest test particle object
      class(rmvs_cb), intent(in)      :: cb   !! RMVS central body particle type
      real(DP), intent(in)            :: dt   !! Step size
   end subroutine rmvs_interp_in

   module subroutine rmvs_interp_out(self, cb, dt)
      implicit none
      class(rmvs_pl), intent(inout)   :: self !! Swiftest test particle object
      class(rmvs_cb), intent(in)      :: cb   !! RMVS central body particle type
      real(DP), intent(in)            :: dt   !! Step size
   end subroutine rmvs_interp_out
   module subroutine rmvs_step_in(self, cb, tp, config, dt)
      implicit none
      class(rmvs_pl),                intent(inout) :: self !! RMVS massive body object
      class(rmvs_cb),                intent(inout) :: cb   !! RMVS central body object
      class(rmvs_tp),                intent(inout)  :: tp   !! RMVS test particle object
      class(swiftest_configuration), intent(in)    :: config  !! Input collection of user-defined configuration parameters 
      real(DP),                      intent(in)    :: dt   !! Step size
   end subroutine rmvs_step_in
   module subroutine rmvs_step_out(self, cb, tp, dt, config)
      implicit none
      class(rmvs_pl),                intent(inout)  :: self !! RMVS massive body object
      class(rmvs_cb),                intent(inout)  :: cb   !! RMVS central body object
      class(rmvs_tp),                intent(inout)  :: tp   !! RMVS test particle object
      real(DP),                      intent(in)     :: dt   !! Step size
      class(swiftest_configuration), intent(in)     :: config  !! Input collection of user-defined configuration parameters 
   end subroutine rmvs_step_out

   module subroutine rmvs_step_out2(cb, pl, tp, t, dt, index, config)
      implicit none
      class(rmvs_cb),                intent(inout) :: cb   !! RMVS central body object
      class(rmvs_pl),                intent(inout)  :: pl   !! RMVS massive body object
      class(rmvs_tp),                intent(inout) :: tp   !! RMVS test particle object
      real(DP),                      intent(in)    :: t     !! Simulation time
      real(DP),                      intent(in)    :: dt    !! Step size
      integer(I4B),                  intent(in)    :: index !! outer substep number within current step
      class(swiftest_configuration), intent(in)    :: config  !! Input collection of user-defined configuration parameters 
   end subroutine rmvs_step_out2   
end interface

   !> Interfaces for all non-type bound whm methods that are implemented in separate submodules 
   interface

      module function rmvs_chk_ind(xr, vr, dt, r2crit) result(lflag)
         implicit none
         real(DP), intent(in)                :: dt, r2crit
         real(DP), dimension(:), intent(in)  :: xr, vr
         logical                             :: lflag
      end function rmvs_chk_ind

      module subroutine rmvs_step_system(cb, pl, tp, config)
         implicit none
         class(rmvs_cb),                intent(inout) :: cb      !! WHM central body object  
         class(rmvs_pl),                intent(inout) :: pl      !! WHM central body object  
         class(rmvs_tp),                intent(inout) :: tp      !! WHM central body object  
         class(swiftest_configuration), intent(in)    :: config  !! Input collection of user-defined configuration parameters 
      end subroutine rmvs_step_system
   end interface

end module rmvs_classes
