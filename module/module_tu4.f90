!**********************************************************************************************************************************
!
!  Unit Name   : module_tu4
!  Unit Type   : module
!  Project     : SWIFTER
!  Package     : module
!  Language    : Fortran 90/95
!
!  Description : Definition of data and structures specific to the 4th-order T + V integrator
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
MODULE module_tu4

     USE module_parameters
     USE module_swifter
     IMPLICIT NONE

     REAL(DP)               :: w0, w1
     REAL(DP), DIMENSION(4) :: asymp, bsymp

     TYPE tu4_pl
          REAL(DP), DIMENSION(NDIM) :: ab       ! total barycentric acceleration
          TYPE(swifter_pl)          :: swifter  ! SWIFTER planet structure
          TYPE(tu4_pl), POINTER     :: prevP    ! pointer to previous TU4 planet
          TYPE(tu4_pl), POINTER     :: nextP    ! pointer to next TU4 planet
     END TYPE tu4_pl

     TYPE tu4_tp
          REAL(DP), DIMENSION(NDIM) :: ab       ! total barycentric acceleration
          TYPE(swifter_tp)          :: swifter  ! SWIFTER test particle structure
          TYPE(tu4_tp), POINTER     :: prevP    ! pointer to previous TU4 test particle
          TYPE(tu4_tp), POINTER     :: nextP    ! pointer to next TU4 test particle
     END TYPE tu4_tp

END MODULE module_tu4
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
