!**********************************************************************************************************************************
!
!  Unit Name   : module_ringmoons
!  Unit Type   : module
!  Project     : SWIFTER
!  Package     : module
!  Language    : Fortran 90/95
!
!  Description : Definition of data and structures specific to the RINGMOONS system
!
!  Input
!    Arguments : N/A
!    Terminal  : N/A
!    File      : N/A
!
!  Output
!    Arguments : N/A
!    Terminal  : N/A
!    File      : N/A
!
!  Invocation  : N/A
!
!  Notes       : 
!
!**********************************************************************************************************************************
!  Author(s)   : David A. Minton
!**********************************************************************************************************************************
module module_ringmoons
   use module_parameters
   use module_symba
   implicit none

   real(DP),public,parameter  :: FEEDING_ZONE_FACTOR = 4.0_DP  ! Size of the feeding zone relative to Hill's sphere

   type ringmoons_ring
      real(DP)     :: r_pdisk             ! disk particle radius
      real(DP)     :: Gm_pdisk            ! disk particle mass
      real(DP)     :: rho_pdisk           ! disk particle density
      integer(I4B) :: N                   ! number of bins in disk
      integer(I4B) :: inside = 1          ! bin id of innermost ring bin (can increase if primary accretes a lot mass through updates)
      real(DP)     :: r_F                 ! outside radius of disk
      real(DP)     :: r_I                 ! inside radius of disk
      real(DP)     :: deltaX              ! variable changed bin width used for viscosity calculations
      real(DP)     :: stability_factor    ! constant factor used to determine the stability of the pde solver
      real(DP)     :: RRL,FRL             ! Rigid and fluid Roche limits
      integer(I4B) :: iRRL,iFRL           ! Indexes of Roche limit bins
      real(DP), dimension(:), allocatable :: r                 ! radial distance of center of bin
      real(DP), dimension(:), allocatable :: rinner            ! inner edge of the bin used to determine if the bin is inside the planet
      real(DP), dimension(:), allocatable :: router            ! outer edge of the bin used to determine if the bin is inside the planet
      real(DP), dimension(:), allocatable :: r_hstar           ! normalized ring Hill's radius
      real(DP), dimension(:), allocatable :: X                 ! distance variable X bin center used for viscosity calculations
      real(DP), dimension(:), allocatable :: X2                ! distance variable X**2 bin center  used for viscosity calculations
      real(DP), dimension(:), allocatable :: deltaA            ! differential surface area of ring
      real(DP), dimension(:), allocatable :: Gm                ! mass of ring particles in bin
      real(DP), dimension(:), allocatable :: Gsigma            ! surface mass density of ring
      real(DP), dimension(:), allocatable :: nu                ! viscocity of the ring
      real(DP), dimension(:), allocatable :: sigma_threshold   ! bounds are set up from lindblad resonance locations at frl and synch
      real(DP), dimension(:), allocatable :: Iz                ! polar moment of inertia of ring bin
      real(DP), dimension(:), allocatable :: w                 ! Keplerian angular velocity of ring bin
      real(DP), dimension(:), allocatable :: torque_to_disk    ! total satellite torque acting on the ring bin
   end type ringmoons_ring

   type ringmoons_seeds ! Satellite "seeds" that eventually turn into SyMBA massive bodies
      integer(I4B)                              :: N            ! Number of satellite seeds
      real(DP)                                  :: Gminit       ! initial mass of seeds
      real(DP), dimension(:), allocatable       :: a            ! Semimajor axis of seed
      real(DP), dimension(:), allocatable       :: Gm           ! Mass of seed
      real(DP), dimension(:), allocatable       :: Rhill        ! Hill's sphere radius of seed
      integer(I4B), dimension(:), allocatable   :: rbin         ! Ring bin location of seed
   end type ringmoons_seeds


end module module_ringmoons
