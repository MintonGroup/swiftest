module rmvs_classes
   !! author: David A. Minton
   !!
   !! Definition of classes and methods specific to the Democratic Heliocentric Method
   !! Partially adapted from David E. Kaufmann's Swifter module: module_whm.f90
   use swiftest_globals
   use swiftest_classes
   use whm_classes
   implicit none

   !********************************************************************************************************************************
   ! rmvs_configuration class definitions and method interfaces
   !*******************************************************************************************************************************
   type, public, extends(whm_configuration) :: rmvs_configuration
   contains
   end type

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
      integer(I4B)              :: nenc    ! number of test particles encountering planet this full rmvs time step
      real(DP), dimension(:, :),  allocatable :: xout    ! interpolated heliocentric planet position for outer encounter
      real(DP), dimension(:, :),  allocatable :: vout    ! interpolated heliocentric planet velocity for outer encounter
      real(DP), dimension(:, :),  allocatable :: xin     ! interpolated heliocentric planet position for inner encounter
      real(DP), dimension(:, :),  allocatable :: vin     ! interpolated heliocentric planet velocity for inner encounter
      real(DP), dimension(:, :),  allocatable :: xpc     ! interpolated planetocentric planet position for inner encounter
      integer(I4B), dimension(:), allocatable :: tpenc1P ! index of first test particle encountering planet 
      !! Note to developers: If you add componenets to this class, be sure to update methods and subroutines that traverse the
      !!    component list, such as rmvs_setup_pl and rmvs_discard_spill_pl
      contains
      procedure, public :: setup        => rmvs_setup_pl    !! Constructor method - Allocates space for number of particles
   end type rmvs_pl

   interface
      module subroutine rmvs_setup_pl(self,n)
         implicit none
         class(rmvs_pl), intent(inout)    :: self !! Swiftest test particle object
         integer, intent(in)             :: n    !! Number of test particles to allocate
      end subroutine rmvs_setup_pl
   end interface

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

   contains
      procedure, public :: setup        => rmvs_setup_tp    !! Constructor method - Allocates space for number of particles
   end type rmvs_tp

   interface
      module subroutine rmvs_setup_tp(self,n)
         implicit none
         class(rmvs_tp), intent(inout)    :: self !! Swiftest test particle object
         integer, intent(in)             :: n    !! Number of test particles to allocate
      end subroutine rmvs_setup_tp
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
      procedure, public :: initialize    => rmvs_setup_system  !! Performs WHM-specific initilization steps, like calculating the Jacobi masses
   end type rmvs_nbody_system
   

   interface
      module subroutine rmvs_setup_system(self, config)
         implicit none
         class(rmvs_nbody_system),      intent(inout) :: self    !! Swiftest system object
         class(swiftest_configuration), intent(inout) :: config  !! Input collection of user-defined configuration parameters 
      end subroutine rmvs_setup_system

   end interface

   !> Interfaces for all non-type bound whm methods that are implemented in separate submodules 
   interface

      module function rmvs_chk(pl, tp, xh, vh, dt, rts) result(lencounter)
         implicit none
         class(rmvs_pl),            intent(in) :: pl      !! WHM central body object  
         class(rmvs_tp),            intent(in) :: tp      !! WHM central body object  
         real(DP), dimension(:, :), intent(in) :: xh, vh  !! planet positions/velocities at beginning of time step
         real(DP),                  intent(in) :: dt  !! time step
         real(DP),                  intent(in) :: rts !! fraction of Hill's sphere radius to use as radius of encounter regio
         logical                    :: lencounter        

      end function rmvs_chk

      module subroutine rmvs_chk_ind(xr, vr, dt, r2crit, iflag)
         implicit none
         real(DP), intent(in)                :: dt, r2crit
         real(DP), dimension(:), intent(in)  :: xr, vr
         integer(I4B),           intent(out) :: iflag
      end subroutine rmvs_chk_ind

      module subroutine rmvs_step_system(cb, pl, tp, config)
         implicit none
         class(rmvs_cb),            intent(inout) :: cb      !! WHM central body object  
         class(rmvs_pl),            intent(inout) :: pl      !! WHM central body object  
         class(rmvs_tp),            intent(inout) :: tp      !! WHM central body object  
         class(rmvs_configuration), intent(in)    :: config  !! Input collection of user-defined configuration parameters 
      end subroutine rmvs_step_system
   end interface

end module rmvs_classes
