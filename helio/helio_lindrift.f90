!**********************************************************************************************************************************
!
!  Unit Name   : helio_lindrift
!  Unit Type   : subroutine
!  Project     : Swiftest
!  Package     : helio
!  Language    : Fortran 90/95
!
!  Description : Perform linear drift of planets due to barycentric momentum of Sun
!
!  Input
!    Arguments : npl          : number of planets
!                swifter_pl1P : pointer to head of Swifter planet structure linked-list
!                dt           : time step
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : swifter_pl1P : pointer to head of Swifter planet structure linked-list
!                pt           : negative barycentric velocity of the Sun
!    Terminal  : none
!    File      : none
!
!  Invocation  : CALL helio_lindrift(npl, swifter_pl1P, dt, pt)
!
!  Notes       : Adapted from Hal Levison's Swift routine helio_lindrift.f
!
!**********************************************************************************************************************************
SUBROUTINE helio_lindrift(npl, swiftest_plA, dt, pt)

! Modules
     USE module_parameters
     USE module_swiftest
     USE module_symba
     USE module_interfaces, EXCEPT_THIS_ONE => helio_lindrift
     IMPLICIT NONE

! Arguments
     INTEGER(I4B), INTENT(IN)                          :: npl
     REAL(DP), INTENT(IN)                              :: dt
     REAL(DP), DIMENSION(NDIM), INTENT(OUT)            :: pt
     TYPE(swiftest_pl), INTENT(INOUT)    :: swiftest_plA

! Internals
     INTEGER(I4B)              :: i

! Added by D. Minton
     REAL(DP) :: ptx,pty,ptz
     REAL(DP),DIMENSION(NDIM) :: pttmp !INTENT(OUT) variables don't play nicely 
                                       !with OpenMP's reduction for some reason


! EDIT THIS PARALLELIZATION

! Executable code
    
     pttmp(:) = (/ 0.0_DP, 0.0_DP, 0.0_DP /)
     DO i = 2, npl
          pttmp(:) = pttmp(:) + swiftest_plA%mass(i)*swiftest_plA%vb(:,i)
     END DO
     pttmp(:) = pttmp(:)/swiftest_plA%mass(1)
     DO i = 2, npl
          swiftest_plA%xh(:,i) = swiftest_plA%xh(:,i) + pttmp(:)*dt
     END DO
     pt(:)=pttmp(:)

     RETURN

END SUBROUTINE helio_lindrift
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
