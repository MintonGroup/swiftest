!**********************************************************************************************************************************
!
!  Unit Name   : helio_getacch_int
!  Unit Type   : subroutine
!  Project     : Swiftest
!  Package     : helio
!  Language    : Fortran 90/95
!
!  Description : Compute direct cross term heliocentric accelerations of planets
!
!  Input
!    Arguments : npl        : number of planets
!                helio_pl1P : pointer to head of helio planet structure linked-list
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : helio_pl1P : pointer to head of helio planet structure linked-list
!    Terminal  : none
!    File      : none
!
!  Invocation  : CALL helio_getacch_int(npl, helio_pl1P)
!
!  Notes       : Adapted from Hal Levison's Swift routine getacch_ah3.f
!
!**********************************************************************************************************************************
SUBROUTINE helio_getacch_int(npl, helio_plA)

! Modules
     USE module_parameters
     USE module_helio
     USE module_interfaces, EXCEPT_THIS_ONE => helio_getacch_int
     IMPLICIT NONE

! Arguments
     INTEGER(I4B), INTENT(IN)       :: npl
     TYPE(helio_pl), INTENT(INOUT)  :: helio_plA

! Internals
     INTEGER(I4B)              :: i, j
     REAL(DP)                  :: rji2, irij3, faci, facj
     REAL(DP), DIMENSION(NDIM) :: dx


! Executable code
     DO i = 2, npl - 1
          DO j = i + 1, npl
               dx(:) = helio_plA%swiftest%xh(:,j) - helio_plA%swiftest%xh(:,i)
               rji2 = DOT_PRODUCT(dx(:), dx(:))
               irij3 = 1.0_DP/(rji2*SQRT(rji2))
               faci = helio_plA%swiftest%mass(i)*irij3
               facj = helio_plA%swiftest%mass(j)*irij3
               helio_plA%ahi(:,i) = helio_plA%ahi(:,i) + facj*dx(:)
               helio_plA%ahi(:,i) = helio_plA%ahi(:,j) - faci*dx(:)
          END DO
     END DO

     RETURN

END SUBROUTINE helio_getacch_int
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
