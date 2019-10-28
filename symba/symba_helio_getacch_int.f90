!**********************************************************************************************************************************
!
!  Unit Name   : symba_helio_getacch_int
!  Unit Type   : subroutine
!  Project     : Swiftest
!  Package     : symba
!  Language    : Fortran 90/95
!
!  Description : Compute direct cross term heliocentric accelerations of planets
!
!  Input
!    Arguments : npl        : number of planets
!                nplm       : number of planets with mass > mtiny
!                helio_pl1P : pointer to head of helio planet structure linked-list
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : helio_pl1P : pointer to head of helio planet structure linked-list
!    Terminal  : none
!    File      : none
!
!  Invocation  : CALL symba_helio_getacch_int(npl, nplm, helio_pl1P)
!
!  Notes       : Adapted from Hal Levison's Swift routine symba5_helio_getacch.f
!
!**********************************************************************************************************************************
SUBROUTINE symba_helio_getacch_int(npl, nplm, helio_plA)

! Modules
     USE module_parameters
     USE module_helio
     USE module_symba
     USE module_interfaces, EXCEPT_THIS_ONE => symba_helio_getacch_int
     IMPLICIT NONE

! Arguments
     INTEGER(I4B), INTENT(IN)                    :: npl, nplm
     TYPE(helio_pl), INTENT(INOUT)               :: helio_plA

! Internals
     INTEGER(I4B)                                :: i, j
     REAL(DP)                                    :: rji2, irij3, faci, facj
     REAL(DP), DIMENSION(NDIM)                   :: dx

! Executable code
     DO i = 2, nplm
          DO j = i + 1, npl
               dx(:) = helio_plA%swiftest%xh(:,j) - helio_plA%swiftest%xh(:,i)
               rji2 = DOT_PRODUCT(dx(:), dx(:))
               irij3 = 1.0_DP/(rji2*SQRT(rji2))
               faci = helio_plA%swiftest%mass(i)*irij3
               facj = helio_plA%swiftest%mass(j)*irij3
               helio_plA%ahi(:,i) = helio_plA%ahi(:,i) + facj*dx(:)
               helio_plA%ahi(:,j) = helio_plA%ahi(:,j) - faci*dx(:)
          END DO
     END DO

     RETURN

END SUBROUTINE symba_helio_getacch_int
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
