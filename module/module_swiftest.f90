!**********************************************************************************************************************************
!
!  Unit Name   : module_swiftest
!  Unit Type   : module
!  Project     : SWIFTEST
!  Package     : module
!  Language    : Fortran 2003 
!
!  Description : Definition of data and structures generic to all integrators
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
!
!  Author(s)   : David A. Minton
!
!**********************************************************************************************************************************
module module_swiftest

   use module_parameters 
   implicit none


   
   !integer(I4B), parameter :: NENMAX = 32767
   !integer(I4B), parameter :: NTENC = 3
   !real(DP), parameter     :: RHSCALE = 6.5_DP
   !real(DP), parameter     :: RSHELL = 0.48075_DP

   type swiftest_pl
       integer(I4B), dimension(:),     allocatable :: id     ! external identifier
       integer(I4B), dimension(:),     allocatable :: status ! status
       real(DP),     dimension(:),     allocatable :: mass   ! mass
       real(DP),     dimension(:),     allocatable :: radius ! radius
       real(DP),     dimension(:),     allocatable :: rhill  ! hill's sphere radius
       real(DP),     dimension(:,:),   allocatable :: xh     ! heliocentric position
       real(DP),     dimension(:,:),   allocatable :: vh     ! heliocentric velocity
       real(DP),     dimension(:,:),   allocatable :: xb     ! barycentric position
       real(DP),     dimension(:,:),   allocatable :: vb     ! barycentric velocity
   end type swiftest_pl

   type swiftest_tp
       integer(I4B), dimension(:),     allocatable :: id     ! external identifier
       integer(I4B), dimension(:),     allocatable :: status ! status
       integer(I4B), dimension(:),     allocatable :: isperi ! perihelion passage flag
       real(DP),     dimension(:),     allocatable :: peri   ! perihelion distance
       real(DP),     dimension(:),     allocatable :: atp    ! semimajor axis following perihelion passage
       real(DP),     dimension(:,:),   allocatable :: xh     ! heliocentric position
       real(DP),     dimension(:,:),   allocatable :: vh     ! heliocentric velocity
       real(DP),     dimension(:,:),   allocatable :: xb     ! barycentric position
       real(DP),     dimension(:,:),   allocatable :: vb     ! barycentric velocity
   end type swiftest_tp


   type swiftest_helio_pl
       real(DP),     dimension(:,:),   allocatable :: ah     ! total heliocentric acceleration
       real(DP),     dimension(:,:),   allocatable :: ahi    ! heliocentric acceleration due to interactions
       type(swiftest_pl)                           :: swiftest  ! swifter planet structure
   end type swiftest_helio_pl

   type swiftest_helio_tp
       real(DP),     dimension(:,:),   allocatable :: ah       ! total heliocentric acceleration
       real(DP),     dimension(:,:),   allocatable :: ahi      ! heliocentric acceleration due to interactions
       type(swiftest_tp)                           :: swiftest  ! swifter test particle structure
   end type swiftest_helio_tp


  type swiftest_symba_pl
       logical(LGT), dimension(:),     allocatable :: lmerged ! flag indicating whether body has merged with another this time step
       integer(I4B), dimension(:),     allocatable :: nplenc  ! number of encounters with other planets this time step
       integer(I4B), dimension(:),     allocatable :: ntpenc  ! number of encounters with test particles this time step
       integer(I4B), dimension(:),     allocatable :: levelg  ! level at which this body should be moved
       integer(I4B), dimension(:),     allocatable :: levelm  ! deepest encounter level achieved this time step
       integer(I4B), dimension(:),     allocatable :: nchild  ! number of children in merger list
       integer(I4B), dimension(:),     allocatable :: isperi  ! perihelion passage flag
       real(DP),     dimension(:),     allocatable :: peri    ! perihelion distance
       real(DP),     dimension(:),     allocatable :: atp     ! semimajor axis following perihelion passage
       type(swiftest_helio_pl)                     :: helio   ! HELIO planet structure
  end type swiftest_symba_pl

  type swiftest_symba_tp
       integer(I4B), dimension(:),     allocatable :: nplenc  ! number of encounters with planets this time step
       integer(I4B), dimension(:),     allocatable :: levelg  ! level at which this particle should be moved
       integer(I4B), dimension(:),     allocatable :: levelm  ! deepest encounter level achieved this time step
       type(swiftest_helio_tp)                     :: helio   ! HELIO test particle structure
  end type swiftest_symba_tp

  type swiftest_symba_plplenc
       logical(LGT), dimension(:),     allocatable :: lvdotr ! relative vdotr flag
       integer(I4B), dimension(:),     allocatable :: status ! status of the interaction
       integer(I4B), dimension(:),     allocatable :: level  ! encounter recursion level
       !TODO: Pointer or arrays?
       !type(symba_pl), POINTER :: pl1P   ! pointer to first planet in encounter
       !type(symba_pl), POINTER :: pl2P   ! pointer to second planet in encounter
  end type swiftest_symba_plplenc

  type swiftest_symba_pltpenc
       logical(LGT), dimension(:),     allocatable :: lvdotr ! relative vdotr flag
       integer(I4B), dimension(:),     allocatable :: status ! status of the interaction
       integer(I4B), dimension(:),     allocatable :: level  ! encounter recursion level
       !TODO: Pointer or arrays?
       !type(symba_pl), POINTER :: plP    ! pointer to planet in encounter
       !type(symba_tp), POINTER :: tpP    ! pointer to test particle in encounter
  end type swiftest_symba_pltpenc

  type swiftest_symba_merger
       integer(I4B), dimension(:),     allocatable :: id     ! external identifier
       integer(I4B), dimension(:),     allocatable :: status ! status
       integer(I4B), dimension(:),     allocatable :: ncomp  ! number of component bodies in this one during this merger
       real(DP),     dimension(:,:),   allocatable :: xh     ! heliocentric position
       real(DP),     dimension(:,:),   allocatable :: vh     ! heliocentric velocity
  end type swiftest_symba_merger

end module module_swiftest
