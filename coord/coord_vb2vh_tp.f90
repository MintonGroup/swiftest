!**********************************************************************************************************************************
!
!  Unit Name   : coord_vb2vh_tp
!  Unit Type   : subroutine
!  Project     : Swiftest
!  Package     : coord
!  Language    : Fortran 90/95
!
!  Description : Convert from barycentric to heliocentric coordinates, active test particle velocities only
!
!  Input
!    Arguments : ntp          : number of active test particles
!                swifter_tp1P : pointer to head of active Swifter test particle structure linked-list
!                vs           : barycentric velocity of the Sun
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : swifter_tp1P : pointer to head of active Swifter test particle structure linked-list
!    Terminal  : none
!    File      : none
!
!  Invocation  : CALL coord_vb2vh_tp(ntp, swifter_tp1P, vs)
!
!  Notes       : Adapted from Hal Levison's Swift routine coord_vb2h_tp.f
!
!**********************************************************************************************************************************
SUBROUTINE coord_vb2vh_tp(ntp, symba_tpA, vs)

! Modules
     USE module_parameters
     USE module_swiftest
     USE module_symba
     USE module_interfaces, EXCEPT_THIS_ONE => coord_vb2vh_tp
     IMPLICIT NONE

! Arguments
     INTEGER(I4B), INTENT(IN)              :: ntp
     REAL(DP), DIMENSION(NDIM), INTENT(IN) :: vs
     TYPE(symba_tp), DIMENSION(:), INTENT(INOUT)             :: symba_tpA

! Internals
     INTEGER(I4B)              :: i

! Executable code
     DO i = 1, ntp
          IF (symba_tpA%helio%swiftest%status(i) == ACTIVE) symba_tpA%helio%swiftest%vh(:,i) = symba_tpA%helio%swiftest%vb(:,i) - vs(:)
     END DO

     RETURN

END SUBROUTINE coord_vb2vh_tp
!**********************************************************************************************************************************
!
!  Author(s)   : David E. Kaufmann (Checked by Jennifer Pouplin & Carlisle Wishard)
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
