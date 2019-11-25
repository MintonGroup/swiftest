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
   implicit none

   real(DP),public,parameter  :: FEEDING_ZONE_FACTOR = 10.0_DP  ! Size of the feeding zone relative to Hill's sphere
   integer(I4B),public,parameter  :: N_DISK_FACTOR = 1000 ! Minimum number of particles in a bin to consider it a fluid disk
   logical(LGT),public        :: DESTRUCTION_EVENT ! A destruction event has occurred when satellite/seed crosses the RRL
   integer(I4B),public        :: DESTRUCTION_COUNTER = 0
   integer(I4B),public,parameter  :: DESTRUCTION_SAVE_FRAMES = 1 ! Number steps to save as output frames after a destruction event
   integer(I4B),public,parameter :: M_MAX = 200
   real(DP),dimension(2:M_MAX,-1:1),public,save :: lapm,dlapm,marr
   real(DP),dimension(2:M_MAX),public,save :: mfac
   real(DP),parameter         :: RK_FACTOR = 0.01_DP

   type ringmoons_ring
      integer(I4B) :: N                   ! number of bins in disk
      integer(I4B) :: inside = 1          ! bin id of innermost ring bin (can increase if primary accretes a lot mass through updates)
      real(DP)     :: r_F                 ! outside radius of disk
      real(DP)     :: X_F                 ! outside radius of disk
      real(DP)     :: r_I                 ! inside radius of disk
      real(DP)     :: X_I                 ! inside radius of disk
      real(DP)     :: deltaX              ! variable changed bin width used for viscosity calculations
      real(DP)     :: RRL,FRL             ! Rigid and fluid Roche limits
      integer(I4B) :: iRRL,iFRL           ! Indexes of Roche limit bins
      real(DP)     :: GMPi,dGMP           ! Original planet mass and an accumulator term
      real(DP)     :: LPi,dLP             ! Original planet angular momentum and an accumulator term
      real(DP)     :: RPi,dRP
      real(DP)     :: rotPi,drotP
      real(DP), dimension(:), allocatable :: r                 ! radial distance of center of bin
      real(DP), dimension(:), allocatable :: r_hstar           ! normalized ring Hill's radius
      real(DP), dimension(:), allocatable :: X                 ! distance variable X bin center used for viscosity calculations
      real(DP), dimension(:), allocatable :: X2                ! distance variable X**2 bin center  used for viscosity calculations
      real(DP), dimension(:), allocatable :: deltaA            ! differential surface area of ring
      real(DP), dimension(:), allocatable :: Gm                ! mass of ring particles in bin
      real(DP), dimension(:), allocatable :: Gsigma            ! surface mass density of ring bin
      real(DP), dimension(:), allocatable :: nu                ! viscocity of the ring bin
      real(DP), dimension(:), allocatable :: Q                 ! Toomre parameter of the ring bin
      real(DP), dimension(:), allocatable :: Iz                ! polar moment of inertia of ring bin
      real(DP), dimension(:), allocatable :: w                 ! Keplerian angular velocity of ring bin
      real(DP), dimension(:), allocatable :: Torque            ! total satellite torque density acting on the ring bin
      real(DP), dimension(:), allocatable :: r_pdisk           ! ring particle radius
      real(DP), dimension(:), allocatable :: Gm_pdisk          ! ring particle mass
      real(DP), dimension(:), allocatable :: rho_pdisk         ! ring particle mass density
      real(DP), dimension(:), allocatable :: vrel_pdisk        ! ring particle relative velocity
      real(DP), dimension(:), allocatable :: tau               ! ring optical depth
   end type ringmoons_ring

   type ringmoons_seeds ! Satellite "seeds" that eventually turn into SyMBA massive bodies
      integer(I4B)                              :: N            ! Number of satellite seeds
      real(DP)                                  :: Gminit       ! initial mass of seeds
      logical(LGT), dimension(:), allocatable   :: active       ! Flag to determine whether this body is active or not
      real(DP), dimension(:), allocatable       :: a            ! Semimajor axis of seed
      real(DP), dimension(:), allocatable       :: Gm           ! Mass of seed
      real(DP), dimension(:), allocatable       :: Rhill        ! Hill's sphere radius of seed
      integer(I4B), dimension(:), allocatable   :: rbin         ! Ring bin location of seed
      integer(I4B), dimension(:), allocatable   :: fz_bin_inner ! Ring bin location of inner edge of seed feeding zone
      integer(I4B), dimension(:), allocatable   :: fz_bin_outer ! Ring bin location of inner edge of seed feeding zone
      real(DP), dimension(:), allocatable       :: Torque       ! Total torque acting on the seed
      real(DP), dimension(:), allocatable       :: Ttide        ! Tidal torque acting on the seed
   end type ringmoons_seeds


end module module_ringmoons
