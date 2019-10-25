!**********************************************************************************************************************************
!
!  Unit Name   : symba_step_interp
!  Unit Type   : subroutine
!  Project     : Swiftest
!  Package     : symba
!  Language    : Fortran 90/95
!
!  Description : Step planets and active test particles ahead in democratic heliocentric coordinates, calling the recursive
!                subroutine to descend to the appropriate level to handle close encounters
!
!  Input
!    Arguments : lextra_force   : logical flag indicating whether to include user-supplied accelerations
!                lclose         : logical flag indicating whether to check for mergers
!                t              : time
!                npl            : number of planets
!                nplm           : number of planets with mass > mtiny
!                nplmax         : maximum allowed number of planets
!                ntp            : number of active test particles
!                ntpmax         : maximum allowed number of test particles
!                symba_pl1P     : pointer to head of SyMBA planet structure linked-list
!                symba_tp1P     : pointer to head of active SyMBA test particle structure linked-list
!                j2rp2          : J2 * R**2 for the Sun
!                j4rp4          : J4 * R**4 for the Sun
!                dt             : time step
!                eoffset        : energy offset (net energy lost in mergers)
!                mtiny          : smallest self-gravitating mass
!                nplplenc       : number of planet-planet encounters
!                npltpenc       : number of planet-test particle encounters
!                plplenc_list   : array of planet-planet encounter structures
!                pltpenc_list   : array of planet-test particle encounter structures
!                nmergeadd      : number of merged planets to add
!                nmergesub      : number of merged planets to subtract
!                mergeadd_list  : array of structures of merged planets to add
!                mergesub_list  : array of structures of merged planets to subtract
!                encounter_file : name of output file for encounters
!                out_type       : binary format of output file
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : symba_pl1P     : pointer to head of SyMBA planet structure linked-list
!                symba_tp1P     : pointer to head of active SyMBA test particle structure linked-list
!                eoffset        : energy offset (net energy lost in mergers)
!                plplenc_list   : array of planet-planet encounter structures
!                pltpenc_list   : array of planet-test particle encounter structures
!                nmergeadd      : number of merged planets to add
!                nmergesub      : number of merged planets to subtract
!                mergeadd_list  : array of structures of merged planets to add
!                mergesub_list  : array of structures of merged planets to subtract
!    Terminal  : none
!    File      : none
!
!  Invocation  : CALL symba_step_interp(lextra_force, lclose, t, npl, nplm, nplmax, ntp, ntpmax, symba_pl1P, symba_tp1P, j2rp2,
!                                       j4rp4, dt, eoffset, mtiny, nplplenc, npltpenc, plplenc_list, pltpenc_list, nmergeadd,
!                                       nmergesub, mergeadd_list, mergesub_list, encounter_file, out_type)
!
!  Notes       : Adapted from Hal Levison's Swift routine symba5_step_interp.f
!
!**********************************************************************************************************************************
SUBROUTINE symba_step_interp(lextra_force, lclose, t, npl, nplm, nplmax, ntp, ntpmax, symba_plA, symba_tpA, j2rp2, j4rp4, dt,   &
     eoffset, mtiny, nplplenc, npltpenc, plplenc_list, pltpenc_list, nmergeadd, nmergesub, mergeadd_list, mergesub_list,          &
     encounter_file, out_type)

! Modules
     USE module_parameters
     USE module_swiftest
     USE module_helio
     USE module_symba
     USE module_interfaces, EXCEPT_THIS_ONE => symba_step_interp
     IMPLICIT NONE

! Arguments
     LOGICAL(LGT), INTENT(IN)                         :: lextra_force, lclose
     INTEGER(I4B), INTENT(IN)                         :: npl, nplm, nplmax, ntp, ntpmax, nplplenc, npltpenc
     INTEGER(I4B), INTENT(INOUT)                      :: nmergeadd, nmergesub
     REAL(DP), INTENT(IN)                             :: t, j2rp2, j4rp4, dt, mtiny
     REAL(DP), INTENT(INOUT)                          :: eoffset
     CHARACTER(*), INTENT(IN)                         :: encounter_file, out_type
     TYPE(symba_pl), DIMENSION(:), INTENT(INOUT)      :: symba_plA
     TYPE(symba_tp), DIMENSION(:), INTENT(INOUT)      :: symba_tpA
     TYPE(symba_plplenc), DIMENSION(:), INTENT(INOUT) :: plplenc_list
     TYPE(symba_pltpenc), DIMENSION(:), INTENT(INOUT) :: pltpenc_list
     TYPE(symba_merger), DIMENSION(:), INTENT(INOUT)  :: mergeadd_list, mergesub_list

! Internals
     LOGICAL(LGT), SAVE                           :: lmalloc = .TRUE.
     INTEGER(I4B)                                 :: i, irec
     REAL(DP)                                     :: dth, msys
     REAL(DP), DIMENSION(NDIM)                    :: ptb, pte
     REAL(DP), DIMENSION(:, :), ALLOCATABLE, SAVE :: xbeg, xend

! Executable code
     IF (lmalloc) THEN
          ALLOCATE(xbeg(NDIM, nplmax), xend(NDIM, nplmax))
          lmalloc = .FALSE.
     END IF
     dth = 0.5_DP*dt
     CALL coord_vh2vb(npl, symba_plA%helio%swiftest, msys)
     CALL helio_lindrift(npl, symba_plA%helio%swiftest, dth, ptb)
     IF (ntp > 0) THEN
          CALL coord_vh2vb_tp(ntp, symba_tpA%helio%swiftest, -ptb)
          CALL helio_lindrift_tp(ntp, symba_tpA%helio%swiftest, dth, ptb) 
          DO i = 2, npl
               xbeg(:, i) = symba_plA%helio%swiftest%xh(:,i)
          END DO
     END IF
     CALL symba_getacch(lextra_force, t, npl, nplm, nplmax, symba_plA, j2rp2, j4rp4, nplplenc, plplenc_list)
     IF (ntp > 0) CALL symba_getacch_tp(lextra_force, t, npl, nplm, nplmax, ntp, ntpmax, symba_plA, symba_tpA, xbeg, j2rp2,     &
          j4rp4, npltpenc, pltpenc_list)
     CALL helio_kickvb(npl, symba_plA%helio, dth)
     IF (ntp > 0) CALL helio_kickvb_tp(ntp, symba_tpA%helio, dth)
     irec = -1
     CALL symba_helio_drift(irec, npl, symba_plA, dt)
     IF (ntp > 0) CALL symba_helio_drift_tp(irec, ntp, symba_tpA, symba_plA%helio%swiftest%mass(1), dt)
     irec = 0
     CALL symba_step_recur(lclose, t, irec, npl, nplm, ntp, symba_plA, symba_tpA, dt, eoffset, nplplenc, npltpenc,              &
          plplenc_list, pltpenc_list, nmergeadd, nmergesub, mergeadd_list, mergesub_list, encounter_file, out_type)
     IF (ntp > 0) THEN
          DO i = 2, npl
               xend(:, i) = symba_plA%helio%swiftest%xh(:,i)
          END DO
     END IF
     CALL symba_getacch(lextra_force, t+dt, npl, nplm, nplmax, symba_plA, j2rp2, j4rp4, nplplenc, plplenc_list)
     IF (ntp > 0) CALL symba_getacch_tp(lextra_force, t+dt, npl, nplm, nplmax, ntp, ntpmax, symba_plA, symba_tpA, xend, j2rp2,  &
          j4rp4, npltpenc, pltpenc_list)
     CALL helio_kickvb(npl, symba_plA%helio, dth)
     IF (ntp > 0) CALL helio_kickvb_tp(ntp, symba_tpA%helio, dth)
     CALL coord_vb2vh(npl, symba_plA%helio%swiftest)
     CALL helio_lindrift(npl, symba_plA%helio%swiftest, dth, pte)
     IF (ntp > 0) THEN
          CALL coord_vb2vh_tp(ntp, symba_tpA%helio%swiftest, -pte)
          CALL helio_lindrift_tp(ntp, symba_tpA%helio%swiftest, dth, pte)
     END IF

     RETURN

END SUBROUTINE symba_step_interp
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
