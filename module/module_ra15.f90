!**********************************************************************************************************************************
!
!  Unit Name   : module_ra15
!  Unit Type   : module
!  Project     : SWIFTER
!  Package     : module
!  Language    : Fortran 90/95
!
!  Description : Definition of data and structures specific to the 15th-order RADAU integrator
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
MODULE module_ra15

     USE module_parameters
     USE module_swifter
     IMPLICIT NONE

     INTEGER(I4B), DIMENSION(8), PARAMETER :: NW = (/ 0, 0, 1, 3, 6, 10, 15, 21 /)
     REAL(DP), DIMENSION(8), PARAMETER     :: H = (/ 0.00000000000000000_DP,                                                      &
                                                     0.05626256053692215_DP,                                                      &
                                                     0.18024069173689236_DP,                                                      &
                                                     0.35262471711316964_DP,                                                      &
                                                     0.54715362633055538_DP,                                                      &
                                                     0.73421017721541053_DP,                                                      &
                                                     0.88532094683909577_DP,                                                      &
                                                     0.97752061356128750_DP /)
     REAL(DP), PARAMETER                   :: ZERO = 0.0_DP, HALF = 0.5_DP, ONE = 1.0_DP, SR = 1.4_DP, PW = 1.0_DP/9.0_DP
     REAL(DP), PARAMETER                   :: DELTARA15 = 1.0E-7_DP
     REAL(DP), DIMENSION(7)                :: u, w
     REAL(DP), DIMENSION(21)               :: c, d, r
     REAL(DP)                              :: eps, w1

     TYPE ra15_pl
          REAL(DP), DIMENSION(NDIM)    :: xbsav    ! barycentric position saved at the start of each RADAU step
          REAL(DP), DIMENSION(NDIM)    :: vbsav    ! barycentric velocity saved at the start of each RADAU step
          REAL(DP), DIMENSION(NDIM)    :: absav    ! barycentric acceleration saved at the start of each RADAU step
          REAL(DP), DIMENSION(NDIM)    :: ab       ! total barycentric acceleration
          REAL(DP), DIMENSION(7, NDIM) :: b        ! B-values
          REAL(DP), DIMENSION(7, NDIM) :: e        ! B(new)-values saved from previous RADAU step
          REAL(DP), DIMENSION(7, NDIM) :: bd       ! corrections to apply to current B(new)-values
          REAL(DP), DIMENSION(7, NDIM) :: g        ! G-values
          TYPE(swifter_pl)             :: swifter  ! SWIFTER planet structure
          TYPE(ra15_pl), POINTER       :: prevP    ! pointer to previous RA15 planet
          TYPE(ra15_pl), POINTER       :: nextP    ! pointer to next RA15 planet
     END TYPE ra15_pl

     TYPE ra15_tp
          REAL(DP), DIMENSION(NDIM)    :: xbsav    ! barycentric position saved at the start of each RADAU step
          REAL(DP), DIMENSION(NDIM)    :: vbsav    ! barycentric velocity saved at the start of each RADAU step
          REAL(DP), DIMENSION(NDIM)    :: absav    ! barycentric acceleration saved at the start of each RADAU step
          REAL(DP), DIMENSION(NDIM)    :: ab       ! total barycentric acceleration
          REAL(DP), DIMENSION(7, NDIM) :: b        ! B-values
          REAL(DP), DIMENSION(7, NDIM) :: e        ! B(new)-values saved from previous RADAU step
          REAL(DP), DIMENSION(7, NDIM) :: bd       ! corrections to apply to current B(new)-values
          REAL(DP), DIMENSION(7, NDIM) :: g        ! G-values
          TYPE(swifter_tp)             :: swifter  ! SWIFTER test particle structure
          TYPE(ra15_tp), POINTER       :: prevP    ! pointer to previous RA15 test particle
          TYPE(ra15_tp), POINTER       :: nextP    ! pointer to next RA15 test particle
     END TYPE ra15_tp

END MODULE module_ra15
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
