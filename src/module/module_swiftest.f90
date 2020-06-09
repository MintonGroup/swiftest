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

   use module_globals
   implicit none



   !integer(I4B), parameter :: NENMAX = 32767
   !integer(I4B), parameter :: NTENC = 3
   !real(DP), parameter     :: RHSCALE = 6.5_DP
   !real(DP), parameter     :: RSHELL = 0.48075_DP

   type swiftest_pl
      integer(I4B)                                :: npl    !! Number of massive bodies
      integer(I4B), dimension(:),     allocatable :: name   !! External identifier (hash)
      integer(I4B), dimension(:),     allocatable :: status !! Status
      real(DP),     dimension(:),     allocatable :: mass   !! Mass
      real(DP),     dimension(:),     allocatable :: radius !! Radius
      real(DP),     dimension(:),     allocatable :: rhill  !! Hill's sphere radius
      real(DP),     dimension(:,:),   allocatable :: xh     !! Heliocentric position
      real(DP),     dimension(:,:),   allocatable :: vh     !! Heliocentric velocity
      real(DP),     dimension(:,:),   allocatable :: xb     !! Barycentric position
      real(DP),     dimension(:,:),   allocatable :: vb     !! Barycentric velocity
   end type swiftest_pl

   type swiftest_tp
      integer(I4B)                                :: ntp    !! Number of testparticles
      integer(I4B), dimension(:),     allocatable :: name   !! External identifier (hash)
      integer(I4B), dimension(:),     allocatable :: status !! Status
      integer(I4B), dimension(:),     allocatable :: isperi !! Perihelion passage flag
      real(DP),     dimension(:),     allocatable :: peri   !! Perihelion distance
      real(DP),     dimension(:),     allocatable :: atp    !! Semimajor axis following perihelion passage
      real(DP),     dimension(:,:),   allocatable :: xh     !! Heliocentric position
      real(DP),     dimension(:,:),   allocatable :: vh     !! Heliocentric velocity
      real(DP),     dimension(:,:),   allocatable :: xb     !! Barycentric position
      real(DP),     dimension(:,:),   allocatable :: vb     !! Barycentric velocity
   end type swiftest_tp

end module module_swiftest
