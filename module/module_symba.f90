!**********************************************************************************************************************************
!
!  Unit Name   : module_symba
!  Unit Type   : module
!  Project     : SWIFTER
!  Package     : module
!  Language    : Fortran 90/95
!
!  Description : Definition of data and structures specific to the Symplectic Massive Body Algorithm
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
MODULE module_symba

     USE module_parameters
     USE module_helio
     IMPLICIT NONE

     INTEGER(I4B), PARAMETER :: NENMAX = 32767
     INTEGER(I4B), PARAMETER :: NTENC = 3
     REAL(DP), PARAMETER     :: RHSCALE = 6.5_DP
     REAL(DP), PARAMETER     :: RSHELL = 0.48075_DP

     TYPE symba_pl
          LOGICAL(LGT)            :: lmerged ! flag indicating whether body has merged with another this time step
          INTEGER(I4B)            :: nplenc  ! number of encounters with other planets this time step
          INTEGER(I4B)            :: ntpenc  ! number of encounters with test particles this time step
          INTEGER(I4B)            :: levelg  ! level at which this body should be moved
          INTEGER(I4B)            :: levelm  ! deepest encounter level achieved this time step
          INTEGER(I4B)            :: nchild  ! number of children in merger list
          INTEGER(I4B)            :: isperi  ! perihelion passage flag
          REAL(DP)                :: peri    ! perihelion distance
          REAL(DP)                :: atp     ! semimajor axis following perihelion passage
          TYPE(helio_pl)          :: helio   ! HELIO planet structure
          TYPE(symba_pl), POINTER :: parentP ! pointer to parent of merger list
          TYPE(symba_pl), POINTER :: childP  ! pointer to next child in merger list
          TYPE(symba_pl), POINTER :: prevP   ! pointer to previous SyMBA planet
          TYPE(symba_pl), POINTER :: nextP   ! pointer to next SyMBA planet
     END TYPE symba_pl

     TYPE symba_tp
          INTEGER(I4B)            :: nplenc  ! number of encounters with planets this time step
          INTEGER(I4B)            :: levelg  ! level at which this particle should be moved
          INTEGER(I4B)            :: levelm  ! deepest encounter level achieved this time step
          TYPE(helio_tp)          :: helio   ! HELIO test particle structure
          TYPE(symba_tp), POINTER :: prevP   ! pointer to previous SyMBA test particle
          TYPE(symba_tp), POINTER :: nextP   ! pointer to next SyMBA test particle
     END TYPE symba_tp

     TYPE symba_plplenc
          LOGICAL(LGT)            :: lvdotr ! relative vdotr flag
          INTEGER(I4B)            :: status ! status of the interaction
          INTEGER(I4B)            :: level  ! encounter recursion level
          TYPE(symba_pl), POINTER :: pl1P   ! pointer to first planet in encounter
          TYPE(symba_pl), POINTER :: pl2P   ! pointer to second planet in encounter
     END TYPE symba_plplenc

     TYPE symba_pltpenc
          LOGICAL(LGT)            :: lvdotr ! relative vdotr flag
          INTEGER(I4B)            :: status ! status of the interaction
          INTEGER(I4B)            :: level  ! encounter recursion level
          TYPE(symba_pl), POINTER :: plP    ! pointer to planet in encounter
          TYPE(symba_tp), POINTER :: tpP    ! pointer to test particle in encounter
     END TYPE symba_pltpenc

     TYPE symba_merger
          INTEGER(I4B)              :: id     ! external identifier
          INTEGER(I4B)              :: status ! status
          INTEGER(I4B)              :: ncomp  ! number of component bodies in this one during this merger
          REAL(DP), DIMENSION(NDIM) :: xh     ! heliocentric position
          REAL(DP), DIMENSION(NDIM) :: vh     ! heliocentric velocity
     END TYPE symba_merger

END MODULE module_symba
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
