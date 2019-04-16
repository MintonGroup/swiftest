!**********************************************************************************************************************************
!
!  Unit Name   : bs_pzextr
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : bs
!  Language    : Fortran 90/95
!
!  Description : Use polynomial extrapolation to evaluate the positions and velocities of planets and test particles across a
!                full time step using an effectively infinitesimal substep by fitting a polynomial to a sequence of estimates
!                with progressively smaller substep values
!
!  Input
!    Arguments : iest    : number of the current estimate in the sequence
!                xest    : substep size for the current estimate
!                npl     : number of planets
!                ntp     : number of active test particles
!                bs_pl1P : pointer to head of BS planet structure linked-list
!                bs_tp1P : pointer to head of active BS test particle structure linked-list
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : bs_pl1P : pointer to head of BS planet structure linked-list
!                bs_tp1P : pointer to head of active BS test particle structure linked-list
!    Terminal  : error message
!    File      : none
!
!  Invocation  : CALL bs_pzextr(iest, xest, npl, ntp, bs_pl1P, bs_tp1P)
!
!  Notes       : Adapted from Numerical Recipes in Fortran 90: The Art of Parallel Scientific Computing, by Press, Teukolsky,
!                Vetterling, and Flannery, 2nd ed., pp. 1305-6
!
!**********************************************************************************************************************************
SUBROUTINE bs_pzextr(iest, xest, npl, ntp, bs_pl1P, bs_tp1P)

! Modules
     USE module_parameters
     USE module_swifter
     USE module_bs
     USE module_interfaces, EXCEPT_THIS_ONE => bs_pzextr
     IMPLICIT NONE

! Arguments
     INTEGER(I4B), INTENT(IN) :: iest, npl, ntp
     REAL(DP), INTENT(IN)     :: xest
     TYPE(bs_pl), POINTER     :: bs_pl1P
     TYPE(bs_tp), POINTER     :: bs_tp1P

! Internals
     INTEGER(I4B)                        :: i, j
     REAL(DP)                            :: delta, f1, f2
     REAL(DP), DIMENSION(NDIM2)          :: d, tmp, q
     REAL(DP), DIMENSION(IEST_MAX), SAVE :: x
     TYPE(bs_pl), POINTER                :: bs_plP
     TYPE(bs_tp), POINTER                :: bs_tpP

! Executable code
     IF (iest > IEST_MAX) THEN
          WRITE(*, *) "bs_pzextr: probable misuse, too much extrapolation"
          CALL util_exit(FAILURE)
     END IF
     x(iest) = xest
     bs_plP => bs_pl1P
     DO i = 1, npl
          bs_plP%yerr(:) = bs_plP%yseq(:)
          bs_plP%y(:) = bs_plP%yseq(:)
          bs_plP => bs_plP%nextP
     END DO
     bs_tpP => bs_tp1P
     DO i = 1, ntp
          bs_tpP%yerr(:) = bs_tpP%yseq(:)
          bs_tpP%y(:) = bs_tpP%yseq(:)
          bs_tpP => bs_tpP%nextP
     END DO
     IF (iest == 1) THEN
          bs_plP => bs_pl1P
          DO i = 1, npl
               bs_plP%qcol(:, 1) = bs_plP%yseq(:)
               bs_plP => bs_plP%nextP
          END DO
          bs_tpP => bs_tp1P
          DO i = 1, ntp
               bs_tpP%qcol(:, 1) = bs_tpP%yseq(:)
               bs_tpP => bs_tpP%nextP
          END DO
     ELSE
          bs_plP => bs_pl1P
          DO i = 1, npl
               d(:) = bs_plP%yseq(:)
               DO j = 1, iest - 1
                    delta = 1.0_DP/(x(iest-j) - xest)
                    f1 = xest*delta
                    f2 = x(iest-j)*delta
                    q(:) = bs_plP%qcol(:, j)
                    bs_plP%qcol(:, j) = bs_plP%yerr(:)
                    tmp(:) = d(:) - q(:)
                    bs_plP%yerr(:) = f1*tmp(:)
                    d(:) = f2*tmp(:)
                    bs_plP%y(:) = bs_plP%y(:) + bs_plP%yerr(:)
               END DO
               bs_plP%qcol(:, iest) = bs_plP%yerr(:)
               bs_plP => bs_plP%nextP
          END DO
          bs_tpP => bs_tp1P
          DO i = 1, ntp
               d(:) = bs_tpP%yseq(:)
               DO j = 1, iest - 1
                    delta = 1.0_DP/(x(iest-j) - xest)
                    f1 = xest*delta
                    f2 = x(iest-j)*delta
                    q(:) = bs_tpP%qcol(:, j)
                    bs_tpP%qcol(:, j) = bs_tpP%yerr(:)
                    tmp(:) = d(:) - q(:)
                    bs_tpP%yerr(:) = f1*tmp(:)
                    d(:) = f2*tmp(:)
                    bs_tpP%y(:) = bs_tpP%y(:) + bs_tpP%yerr(:)
               END DO
               bs_tpP%qcol(:, iest) = bs_tpP%yerr(:)
               bs_tpP => bs_tpP%nextP
          END DO
     END IF

     RETURN

END SUBROUTINE bs_pzextr
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
