!**********************************************************************************************************************************
!
!  Unit Name   : module_rmvs
!  Unit Type   : module
!  Project     : SWIFTER
!  Package     : module
!  Language    : Fortran 90/95
!
!  Description : Definition of data and structures specific to the Regularized Mixed Variable Symplectic algorithm
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
MODULE module_rmvs

     USE module_parameters
     USE module_whm
     IMPLICIT NONE

     INTEGER(I4B), PARAMETER :: NTENC = 10
     INTEGER(I4B), PARAMETER :: NTPHENC = 3
     INTEGER(I4B), PARAMETER :: NTPENC = NTENC*NTPHENC
     REAL(DP), PARAMETER     :: RHSCALE = 3.5_DP
     REAL(DP), PARAMETER     :: RHPSCALE = 1.0_DP
     REAL(DP), PARAMETER     :: FACQDT = 2.0_DP

     TYPE rmvs_pl
          INTEGER(I4B)                         :: nenc    ! number of test particles encountering planet this full RMVS time step
          REAL(DP), DIMENSION(NDIM, 0:NTENC)   :: xout    ! interpolated heliocentric planet position for outer encounter
          REAL(DP), DIMENSION(NDIM, 0:NTENC)   :: vout    ! interpolated heliocentric planet velocity for outer encounter
          REAL(DP), DIMENSION(NDIM, 0:NTPHENC) :: xin     ! interpolated heliocentric planet position for inner encounter
          REAL(DP), DIMENSION(NDIM, 0:NTPHENC) :: vin     ! interpolated heliocentric planet velocity for inner encounter
          REAL(DP), DIMENSION(NDIM, 0:NTPHENC) :: xpc     ! interpolated planetocentric planet position for inner encounter
          REAL(DP), DIMENSION(NDIM, 0:NTPHENC) :: aobl    ! barycentric acceleration due to central body oblateness during
                                                          ! inner encounter
          TYPE(rmvs_tp), POINTER               :: tpenc1P ! pointer to first test particle encountering planet
          TYPE(whm_pl)                         :: whm     ! WHM planet structure
          TYPE(rmvs_pl), POINTER               :: prevP   ! pointer to previous RMVS planet
          TYPE(rmvs_pl), POINTER               :: nextP   ! pointer to next RMVS planet
     END TYPE rmvs_pl

     TYPE rmvs_tp
          INTEGER(I4B)              :: isperi ! planetocentric pericenter passage flag (not persistent for a set of inner
                                              ! encounter steps)
          LOGICAL(LGT)              :: lperi  ! planetocentric pericenter passage flag (persistent for a full RMVS time step)
          REAL(DP)                  :: peri   ! planetocentric pericenter distance (set to smallest pericenter distance achieved
                                              ! over a full RMVS time step)
          REAL(DP), DIMENSION(NDIM) :: xpc    ! planetocentric position
          REAL(DP), DIMENSION(NDIM) :: vpc    ! planetocentric velocity
          REAL(DP), DIMENSION(NDIM) :: apc    ! total planetocentric acceleration
          TYPE(rmvs_pl), POINTER    :: plperP ! pointer to planet associated with pericenter distance peri (persistent over a
                                              ! full RMVS time step)
          TYPE(rmvs_pl), POINTER    :: plencP ! pointer to planet that test particle is encountering (not persistent for a full
                                              ! RMVS time step)
          TYPE(rmvs_tp), POINTER    :: tpencP ! pointer to next test particle encountering planet
          TYPE(whm_tp)              :: whm    ! WHM test particle structure
          TYPE(rmvs_tp), POINTER    :: prevP  ! pointer to previous RMVS test particle
          TYPE(rmvs_tp), POINTER    :: nextP  ! pointer to next RMVS test particle
     END TYPE rmvs_tp

END MODULE module_rmvs
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
