!**********************************************************************************************************************************
!
!  Unit Name   : symba_step_helio_pl
!  Unit Type   : subroutine
!  Project     : Swiftest
!  Package     : symba
!  Language    : Fortran 90/95
!
!  Description : Step planets ahead in democratic heliocentric coordinates
!
!  Input
!    Arguments : lfirst       : logical flag indicating whether current invocation is the first
!                lextra_force : logical flag indicating whether to include user-supplied accelerations
!                t            : time
!                npl          : number of planets
!                nplm         : number of planets with mass > mtiny
!                nplmax       : maximum allowed number of planets
!                helio_pl1P   : pointer to head of helio planet structure linked-list
!                j2rp2        : J2 * R**2 for the Sun
!                j4rp4        : J4 * R**4 for the Sun
!                dt           : time step
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : lfirst       : logical flag indicating whether current invocation is the first
!                helio_pl1P   : pointer to head of helio planet structure linked-list
!                xbeg         : heliocentric positions of planets with mass > mtiny prior to Kepler drift
!                xend         : heliocentric positions of planets with mass > mtiny after Kepler drift
!                ptb          : negative barycentric velocity of the Sun prior to the first kick
!                pte          : negative barycentric velocity of the Sun after the second kick
!    Terminal  : none
!    File      : none
!
!  Invocation  : CALL symba_step_helio_pl(lfirst, lextra_force, t, npl, nplm, nplmax, helio_pl1P, j2rp2, j4rp4, dt, xbeg, xend,
!                                         ptb, pte)
!
!  Notes       : Adapted from Hal Levison's Swift routines symba5_step_helio.f and helio_step_pl.f
!
!**********************************************************************************************************************************
SUBROUTINE symba_step_helio_pl(lfirst, lextra_force, t, npl, nplm, nplmax, helio_plA, j2rp2, j4rp4, dt, xbeg, xend, ptb, pte)

! Modules
     USE swiftest
     USE module_swiftest
     USE module_helio
     USE module_symba
     USE module_interfaces, EXCEPT_THIS_ONE => symba_step_helio_pl
     IMPLICIT NONE

! Arguments
     LOGICAL(LGT), INTENT(IN)                     :: lextra_force
     LOGICAL(LGT), INTENT(INOUT)                  :: lfirst
     INTEGER(I4B), INTENT(IN)                     :: npl, nplm, nplmax
     REAL(DP), INTENT(IN)                         :: t, j2rp2, j4rp4, dt
     REAL(DP), DIMENSION(NDIM, nplm), INTENT(OUT) :: xbeg, xend
     REAL(DP), DIMENSION(NDIM), INTENT(OUT)       :: ptb, pte
     TYPE(helio_pl), INTENT(INOUT)                :: helio_plA

! Internals
     LOGICAL(LGT)              :: lflag
     INTEGER(I4B)              :: i
     REAL(DP)                  :: dth, msys

! Executable code
     dth = 0.5_DP*dt
     lflag = lfirst
     IF (lfirst) THEN
          CALL coord_vh2vb(npl, helio_plA%swiftest, msys) 
          lfirst = .FALSE.
     END IF
     
     CALL helio_lindrift(npl, helio_plA%swiftest, dth, ptb)

     CALL symba_helio_getacch(lflag, lextra_force, t, npl, nplm, nplmax, helio_plA, j2rp2, j4rp4) 
     lflag = .TRUE.

     CALL helio_kickvb(npl, helio_plA, dth)

     DO i = 2, nplm
          xbeg(:, i) = helio_plA%swiftest%xh(:,i)
     END DO
     CALL helio_drift(npl, helio_plA%swiftest, dt) 

     DO i = 2, nplm
          xend(:, i) = helio_plA%swiftest%xh(:,i)
     END DO
     CALL symba_helio_getacch(lflag, lextra_force, t+dt, npl, nplm, nplmax, helio_plA, j2rp2, j4rp4) 

     CALL helio_kickvb(npl, helio_plA, dth)

     CALL helio_lindrift(npl, helio_plA%swiftest, dth, pte)

     CALL coord_vb2vh(npl, helio_plA%swiftest)

     RETURN

END SUBROUTINE symba_step_helio_pl
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
