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

   use swiftest 
   implicit none


   
   !integer(I4B), parameter :: NENMAX = 32767
   !integer(I4B), parameter :: NTENC = 3
   !real(DP), parameter     :: RHSCALE = 6.5_DP
   !real(DP), parameter     :: RSHELL = 0.48075_DP

   type swiftest_pl
       integer(I4B), dimension(:),     allocatable :: name     ! external identifier (hash)
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
       integer(I4B), dimension(:),     allocatable :: name     ! external identifier (hash)
       integer(I4B), dimension(:),     allocatable :: status ! status
       integer(I4B), dimension(:),     allocatable :: isperi ! perihelion passage flag
       real(DP),     dimension(:),     allocatable :: peri   ! perihelion distance
       real(DP),     dimension(:),     allocatable :: atp    ! semimajor axis following perihelion passage
       real(DP),     dimension(:,:),   allocatable :: xh     ! heliocentric position
       real(DP),     dimension(:,:),   allocatable :: vh     ! heliocentric velocity
       real(DP),     dimension(:,:),   allocatable :: xb     ! barycentric position
       real(DP),     dimension(:,:),   allocatable :: vb     ! barycentric velocity
   end type swiftest_tp

end module module_swiftest
