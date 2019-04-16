!**********************************************************************************************************************************
!
!  Unit Name   : module_bs
!  Unit Type   : module
!  Project     : SWIFTER
!  Package     : module
!  Language    : Fortran 90/95
!
!  Description : Definition of data and structures specific to the Bulirsch-Stoer integrator
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
MODULE module_bs

     USE module_parameters
     USE module_swifter
     IMPLICIT NONE

     INTEGER(I4B), PARAMETER :: IEST_MAX = 7
     REAL(DP),     PARAMETER :: DELTABS  = 1.0E-7_DP
     REAL(DP),     PARAMETER :: TINYBS   = 1.0E-30_DP

     TYPE bs_pl
          REAL(DP), DIMENSION(NDIM2)           :: y       ! dependent variable vector
          REAL(DP), DIMENSION(NDIM2)           :: dydx    ! derivative of dependent variable vector
          REAL(DP), DIMENSION(NDIM2)           :: ysav    ! dependent variable vector (saved at bsstep start)
          REAL(DP), DIMENSION(NDIM2)           :: dydxsav ! derivative of dependent variable vector (saved at bsstep start)
          REAL(DP), DIMENSION(NDIM2)           :: yscal   ! error scaling vector
          REAL(DP), DIMENSION(NDIM2)           :: ym      ! intermediate work vector for mmid
          REAL(DP), DIMENSION(NDIM2)           :: yseq    ! current sequence estimate of final dependent variable vector
          REAL(DP), DIMENSION(NDIM2)           :: yerr    ! extrapolation error of final dependent variable vector
          REAL(DP), DIMENSION(NDIM2, IEST_MAX) :: qcol    ! extrapolation tableau
          REAL(DP), DIMENSION(NDIM)            :: ab      ! total barycentric acceleration
          TYPE(swifter_pl)                     :: swifter ! SWIFTER planet structure
          TYPE(bs_pl), POINTER                 :: prevP   ! pointer to previous BS planet
          TYPE(bs_pl), POINTER                 :: nextP   ! pointer to next BS planet
     END TYPE bs_pl

     TYPE bs_tp
          REAL(DP), DIMENSION(NDIM2)           :: y       ! dependent variable vector
          REAL(DP), DIMENSION(NDIM2)           :: dydx    ! derivative of dependent variable vector
          REAL(DP), DIMENSION(NDIM2)           :: ysav    ! dependent variable vector (saved at bsstep start)
          REAL(DP), DIMENSION(NDIM2)           :: dydxsav ! derivative of dependent variable vector (saved at bsstep start)
          REAL(DP), DIMENSION(NDIM2)           :: yscal   ! error scaling vector
          REAL(DP), DIMENSION(NDIM2)           :: ym      ! intermediate work vector for mmid
          REAL(DP), DIMENSION(NDIM2)           :: yseq    ! current sequence estimate of final dependent variable vector
          REAL(DP), DIMENSION(NDIM2)           :: yerr    ! extrapolation error of final dependent variable vector
          REAL(DP), DIMENSION(NDIM2, IEST_MAX) :: qcol    ! extrapolation tableau
          REAL(DP), DIMENSION(NDIM)            :: ab      ! total barycentric acceleration
          TYPE(swifter_tp)                     :: swifter ! SWIFTER test particle structure
          TYPE(bs_tp), POINTER                 :: prevP   ! pointer to previous BS test particle
          TYPE(bs_tp), POINTER                 :: nextP   ! pointer to next BS test particle
     END TYPE bs_tp

END MODULE module_bs
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
