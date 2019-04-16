!**********************************************************************************************************************************
!
!  Unit Name   : bs_mmid
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : bs
!  Language    : Fortran 90/95
!
!  Description : Perform a modified midpoint step to advance planet and test particle positions and velocities one time step
!
!  Input
!    Arguments : npl          : number of planets
!                nplmax       : maximum allowed number of planets
!                ntp          : number of active test particles
!                ntpmax       : maximum allowed number of test particles
!                bs_pl1P      : pointer to head of BS planet structure linked-list
!                bs_tp1P      : pointer to head of active BS test particle structure linked-list
!                xs           : independent variable (time)
!                htot         : total time step to take
!                nstep        : number of substeps to use
!                j2rp2        : J2 * R**2 for the Sun
!                j4rp4        : J4 * R**4 for the Sun
!                lextra_force : logical flag indicating whether to include user-supplied accelerations
!                derivs       : name of subroutine used to compute planet and test particle velocities and accelerations
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : bs_pl1P      : pointer to head of BS planet structure linked-list
!                bs_tp1P      : pointer to head of active BS test particle structure linked-list
!    Terminal  : none
!    File      : none
!
!  Invocation  : CALL bs_mmid(npl, nplmax, ntp, ntpmax, bs_pl1P, bs_tp1P, xs, htot, nstep, j2rp2, j4rp4, lextra_force, derivs)
!
!  Notes       : Adapted from Numerical Recipes in Fortran 90: The Art of Parallel Scientific Computing, by Press, Teukolsky,
!                Vetterling, and Flannery, 2nd ed., pp. 1302-3
!
!**********************************************************************************************************************************
SUBROUTINE bs_mmid(npl, nplmax, ntp, ntpmax, bs_pl1P, bs_tp1P, xs, htot, nstep, j2rp2, j4rp4, lextra_force, derivs)

! Modules
     USE module_parameters
     USE module_bs
     USE module_nrutil
     USE module_interfaces, EXCEPT_THIS_ONE => bs_mmid
     IMPLICIT NONE

! Arguments
     LOGICAL(LGT), INTENT(IN) :: lextra_force
     INTEGER(I4B), INTENT(IN) :: npl, nplmax, ntp, ntpmax, nstep
     REAL(DP), INTENT(IN)     :: xs, htot, j2rp2, j4rp4
     TYPE(bs_pl), POINTER     :: bs_pl1P
     TYPE(bs_tp), POINTER     :: bs_tp1P
     INTERFACE
          SUBROUTINE derivs(lextra_force, t, npl, nplmax, ntp, ntpmax, bs_pl1P, bs_tp1P, j2rp2, j4rp4)
               USE module_parameters
               USE module_bs
               IMPLICIT NONE
               LOGICAL(LGT), INTENT(IN) :: lextra_force
               INTEGER(I4B), INTENT(IN) :: npl, nplmax, ntp, ntpmax
               REAL(DP), INTENT(IN)     :: t, j2rp2, j4rp4
               TYPE(bs_pl), POINTER     :: bs_pl1P
               TYPE(bs_tp), POINTER     :: bs_tp1P
          END SUBROUTINE derivs
     END INTERFACE

! Internals
     INTEGER(I4B)               :: i, n, ndum
     REAL(DP)                   :: h, h2, x
     REAL(DP), DIMENSION(NDIM2) :: dum
     TYPE(bs_pl), POINTER       :: bs_plP
     TYPE(bs_tp), POINTER       :: bs_tpP

! Executable code
     h = htot/nstep
     bs_plP => bs_pl1P
     DO i = 1, npl
          bs_plP%ym(:) = bs_plP%ysav(:)
          bs_plP%y(:) = bs_plP%ym(:) + h*bs_plP%dydxsav(:)
          bs_plP => bs_plP%nextP
     END DO
     bs_tpP => bs_tp1P
     DO i = 1, ntp
          bs_tpP%ym(:) = bs_tpP%ysav(:)
          bs_tpP%y(:) = bs_tpP%ym(:) + h*bs_tpP%dydxsav(:)
          bs_tpP => bs_tpP%nextP
     END DO
     x = xs + h
     CALL derivs(lextra_force, x, npl, nplmax, ntp, ntpmax, bs_pl1P, bs_tp1P, j2rp2, j4rp4)
     h2 = 2.0_DP*h
     DO n = 2, nstep
          bs_plP => bs_pl1P
          DO i = 1, npl
               dum(:) = bs_plP%ym(:)
               bs_plP%ym(:) = bs_plP%y(:)
               bs_plP%y(:) = dum(:)
               bs_plP%y(:) = bs_plP%y(:) + h2*bs_plP%dydx(:)
               bs_plP => bs_plP%nextP
          END DO
          bs_tpP => bs_tp1P
          DO i = 1, ntp
               dum(:) = bs_tpP%ym(:)
               bs_tpP%ym(:) = bs_tpP%y(:)
               bs_tpP%y(:) = dum(:)
               bs_tpP%y(:) = bs_tpP%y(:) + h2*bs_tpP%dydx(:)
               bs_tpP => bs_tpP%nextP
          END DO
          x = x + h
          CALL derivs(lextra_force, x, npl, nplmax, ntp, ntpmax, bs_pl1P, bs_tp1P, j2rp2, j4rp4)
     END DO
     bs_plP => bs_pl1P
     DO i = 1, npl
          bs_plP%y(:) = 0.5_DP*(bs_plP%ym(:) + bs_plP%y(:) + h*bs_plP%dydx(:))
          bs_plP => bs_plP%nextP
     END DO
     bs_tpP => bs_tp1P
     DO i = 1, ntp
          bs_tpP%y(:) = 0.5_DP*(bs_tpP%ym(:) + bs_tpP%y(:) + h*bs_tpP%dydx(:))
          bs_tpP => bs_tpP%nextP
     END DO

     RETURN

END SUBROUTINE bs_mmid
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
