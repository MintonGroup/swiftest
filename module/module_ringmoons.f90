!**********************************************************************************************************************************
!
!  Unit Name   : module_ringmoons
!  Unit Type   : module
!  Project     : SWIFTER
!  Package     : module
!  Language    : Fortran 90/95
!
!  Description : Definition of data and structures specific to the RING-MOONS system
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
MODULE module_ringmoons
     USE module_parameters
     USE module_symba
     IMPLICIT NONE

     TYPE ringmoons_ring
         REAL(DP)     :: r_pdisk     ! disk particle size (radius)
         REAL(DP)     :: m_pdisk     ! disk particle size (mass)
         INTEGER(I4B) :: N           ! number of bins in disk
         INTEGER(I4B) :: inside = 0  ! bin id of innermost ring bin (can increase if primary accretes a lot mass through updates)
         REAL(DP)     :: r_F         ! outside radius of disk
         REAL(DP)     :: r_I         ! inside radius of disk
         REAL(DP)     :: deltar      ! width of a bin
         REAL(DP)     :: deltaX      ! variable changed bin width used for viscosity calculations
         REAL(DP), dimension(:), allocatable :: r        ! Radial distance of center of bin
         REAL(DP), dimension(:), allocatable :: X        ! variable changed bin center used for viscosity calculations
         REAL(DP), dimension(:), allocatable :: R_P      ! Radial distance of center of bin relative to central body radius
         REAL(DP), dimension(:), allocatable :: deltaA   ! Differential surface area of ring
         REAL(DP), dimension(:), allocatable :: m        ! mass of ring particles in bin
         REAL(DP), dimension(:), allocatable :: sigma    ! Surface mass density of ring
         REAL(DP), dimension(:), allocatable :: nu       ! viscocity of the ring
         REAL(DP), dimension(:), allocatable :: sigma_threshold ! bounds are set up from lindblad resonance locations at FRL and synch
         REAL(DP), dimension(:), allocatable :: RR
         REAL(DP), dimension(:), allocatable :: I
         REAL(DP), dimension(:), allocatable :: w
         REAL(DP), dimension(:), allocatable :: Torque_to_disk
     END TYPE ringmoons_ring

END MODULE module_ringmoons
!**********************************************************************************************************************************
!
!  Author(s)   : David E. Kaufmann
!
!  Revision Control System (RCS) Information
!
!  Source File : $RCSfile$
!  Full Path   : $Source$
!  Revision    : $Revision$
!  Date        : $Date$
!  Programmer  : $Author$
!  Locked By   : $Locker$
!  State       : $State$
!
!  Modification History:
!
!  $Log$
!**********************************************************************************************************************************
