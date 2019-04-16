!**********************************************************************************************************************************
!
!  Unit Name   : module_whm
!  Unit Type   : module
!  Project     : SWIFTER
!  Package     : module
!  Language    : Fortran 90/95
!
!  Description : Definition of data and structures specific to the Wisdom-Holman Mapping
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
MODULE module_whm

     USE module_parameters
     USE module_swifter
     IMPLICIT NONE

     REAL(DP), DIMENSION(NDIM)      :: ah0     ! zeroth term of heliocentric acceleration

     TYPE whm_pl
          REAL(DP)                  :: eta     ! Jacobi mass
          REAL(DP), DIMENSION(NDIM) :: xj      ! Jacobi position
          REAL(DP), DIMENSION(NDIM) :: vj      ! Jacobi velocity
          REAL(DP), DIMENSION(NDIM) :: ah1     ! first term of heliocentric acceleration
          REAL(DP), DIMENSION(NDIM) :: ah2     ! second term of heliocentric acceleration
          REAL(DP), DIMENSION(NDIM) :: ah3     ! third term of heliocentric acceleration
          REAL(DP), DIMENSION(NDIM) :: ah      ! total heliocentric acceleration
          TYPE(swifter_pl)          :: swifter ! SWIFTER planet structure
          TYPE(whm_pl), POINTER     :: prevP   ! pointer to previous WHM planet
          TYPE(whm_pl), POINTER     :: nextP   ! pointer to next WHM planet
     END TYPE whm_pl

     TYPE whm_tp
          REAL(DP), DIMENSION(NDIM) :: ah      ! total heliocentric acceleration
          TYPE(swifter_tp)          :: swifter ! SWIFTER test particle structure
          TYPE(whm_tp), POINTER     :: prevP   ! pointer to previous WHM test particle
          TYPE(whm_tp), POINTER     :: nextP   ! pointer to next WHM test particle
     END TYPE whm_tp

END MODULE module_whm
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
