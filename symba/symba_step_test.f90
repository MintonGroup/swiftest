!**********************************************************************************************************************************
!
!  Unit Name   : symba_step
!  Unit Type   : subroutine
!  Project     : Swiftest
!  Package     : symba
!  Language    : Fortran 90/95
!
!  Description : Step planets and active test particles ahead in democratic heliocentric coordinates, descending the recursive
!                branch if necessary to handle possible close encounters
!
!  Input
!    Arguments : lfirst         : logical flag indicating whether current invocation is the first
!                lextra_force   : logical flag indicating whether to include user-supplied accelerations
!                lclose         : logical flag indicating whether to check for mergers
!                t              : time
!                npl            : number of planets
!                nplmax         : maximum allowed number of planets
!                ntp            : number of active test particles
!                ntpmax         : maximum allowed number of test particles
!                symba_pl1P     : pointer to head of SyMBA planet structure linked-list
!                symba_tp1P     : pointer to head of active SyMBA test particle structure linked-list
!                j2rp2          : J2 * R**2 for the Sun
!                j4rp4          : J4 * R**4 for the Sun
!                dt             : time step
!                nplplenc       : number of planet-planet encounters
!                npltpenc       : number of planet-test particle encounters
!                plplenc_list   : array of planet-planet encounter structures
!                pltpenc_list   : array of planet-test particle encounter structures
!                nmergeadd      : number of merged planets to add
!                nmergesub      : number of merged planets to subtract
!                mergeadd_list  : array of structures of merged planets to add
!                mergesub_list  : array of structures of merged planets to subtract
!                eoffset        : energy offset (net energy lost in mergers)
!                mtiny          : smallest self-gravitating mass
!                encounter_file : name of output file for encounters
!                out_type       : binary format of output file
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : lfirst         : logical flag indicating whether current invocation is the first
!                symba_pl1P     : pointer to head of SyMBA planet structure linked-list
!                symba_tp1P     : pointer to head of active SyMBA test particle structure linked-list
!                nplplenc       : number of planet-planet encounters
!                npltpenc       : number of planet-test particle encounters
!                plplenc_list   : array of planet-planet encounter structures
!                pltpenc_list   : array of planet-test particle encounter structures
!                nmergeadd      : number of merged planets to add
!                nmergesub      : number of merged planets to subtract
!                mergeadd_list  : array of structures of merged planets to add
!                mergesub_list  : array of structures of merged planets to subtract
!                eoffset        : energy offset (net energy lost in mergers)
!    Terminal  : error message
!    File      : none
!
!  Invocation  : CALL symba_step(lfirst, lextra_force, lclose, t, npl, nplmax, ntp, ntpmax, symba_pl1P, symba_tp1P, j2rp2, j4rp4,
!                                dt, nplplenc, npltpenc, plplenc_list, pltpenc_list, nmergeadd, nmergesub, mergeadd_list,
!                                mergesub_list, eoffset, mtiny, encounter_file, out_type)
!
!  Notes       : Adapted from Hal Levison's Swift routine symba5_step_pl.f
!
!**********************************************************************************************************************************
SUBROUTINE symba_step_test(lfirst, lextra_force, lclose, t, npl, nplmax, ntp, ntpmax, symba_plA, symba_tpA, j2rp2, j4rp4, dt,        &
     nplplenc, npltpenc, plplenc_list, pltpenc_list, nmergeadd, nmergesub, mergeadd_list, mergesub_list, eoffset, mtiny,          &
     encounter_file, out_type)

! Modules
     USE module_parameters
     USE module_swiftest
     USE module_helio
     USE module_symba
     USE module_interfaces, EXCEPT_THIS_ONE => symba_step
     IMPLICIT NONE

! Arguments
     LOGICAL(LGT), INTENT(IN)                         :: lextra_force, lclose
     LOGICAL(LGT), INTENT(INOUT)                      :: lfirst
     INTEGER(I4B), INTENT(IN)                         :: npl, nplmax, ntp, ntpmax
     INTEGER(I4B), INTENT(INOUT)                      :: nplplenc, npltpenc, nmergeadd, nmergesub
     REAL(DP), INTENT(IN)                             :: t, j2rp2, j4rp4, dt, mtiny
     REAL(DP), INTENT(INOUT)                          :: eoffset
     CHARACTER(*), INTENT(IN)                         :: encounter_file, out_type
     TYPE(symba_pl), DIMENSION(:), INTENT(INOUT)      :: symba_plA
     TYPE(symba_tp), DIMENSION(:), INTENT(INOUT)      :: symba_tpA
     TYPE(symba_plplenc), DIMENSION(:), INTENT(INOUT) :: plplenc_list
     TYPE(symba_pltpenc), DIMENSION(:), INTENT(INOUT) :: pltpenc_list
     TYPE(symba_merger), DIMENSION(:), INTENT(INOUT)  :: mergeadd_list, mergesub_list

! Internals
     LOGICAL(LGT)              :: lencounter, lvdotr
     INTEGER(I4B)              :: i, j, irec, nplm
     REAL(DP), DIMENSION(NDIM) :: xr, vr
     !TYPE(swifter_pl), POINTER :: swifter_pliP, swifter_pljP
     !Added by D. Minton
     !TYPE(swifter_pl), POINTER :: swifter_pl1P
     !TYPE(swifter_tp), POINTER :: swifter_tp1P
     !^^^^^^^^^^^^^^^^^^
     !TYPE(swifter_tp), POINTER :: swifter_tpP
     !TYPE(helio_pl), POINTER   :: helio_pl1P
     !TYPE(helio_tp), POINTER   :: helio_tp1P
     !TYPE(symba_pl), POINTER   :: symba_pliP, symba_pljP
     !TYPE(symba_tp), POINTER   :: symba_tpP

! Executable code
     !symba_pliP => symba_pl1P
     !helio_pl1P => symba_pl1P%helio
     !Added by D. Minton
     !swifter_pl1P => helio_pl1P%swifter
     !IF (ALLOCATED(symba_pl1P%symba_plPA)) DEALLOCATE(symba_pl1P%symba_plPA)
     !ALLOCATE(symba_pl1P%symba_plPA(npl))
     !IF (ALLOCATED(swifter_pl1P%swifter_plPA)) DEALLOCATE(swifter_pl1P%swifter_plPA)
     !ALLOCATE(swifter_pl1P%swifter_plPA(npl))
     !IF (ALLOCATED(helio_pl1P%helio_plPA)) DEALLOCATE(helio_pl1P%helio_plPA)
     !ALLOCATE(helio_pl1P%helio_plPA(npl))
     !IF (ntp>0) THEN
     !   IF (ALLOCATED(symba_tp1P%symba_tpPA)) DEALLOCATE(symba_tp1P%symba_tpPA)
     !   ALLOCATE(symba_tp1P%symba_tpPA(ntp))
     !END IF
     !^^^^^^^^^^^^^^^^^^
    !DO i = 1, npl
    !      symba_plA%nplenc = 0
    !      symba_tpA%ntpenc = 0
    !      symba_plA%levelg = -1
    !      symba_tpA%levelm = -1
    !      ! Added by D. Minton
    !      symba_pl1P%symba_plPA(i)%thisP => symba_pliP
    !      helio_pl1p%helio_plPA(i)%thisP => symba_pliP%helio
    !      swifter_pl1P%swifter_plPA(i)%thisP => symba_pliP%helio%swifter
    !      !^^^^^^^^^^^^
    !      symba_pliP => symba_pliP%nextP
    ! END DO
     !symba_tpP => symba_tp1P

     symba_plA%nplenc(:) = 0
     symba_plA%ntpenc(:) = 0
     symba_plA%levelg(:) = -1
     symba_plA%levelm(:) = -1
     symba_tpA%nplenc(:) = 0 
     symba_tpA%levelg(:) = -1
     symba_tpA%levelm(:) = -1


     !THERE SHOULD BE SOME PARALLEL BITS IN HERE

     !DO i = 1, ntp
     !     symba_tpP%nplenc = 0
     !     symba_tpP%levelg = -1
     !     symba_tpP%levelm = -1
     !     ! Added by D. Minton
     !     symba_tp1P%symba_tpPA(i)%thisP => symba_tpP
     !     !^^^^^^^^^^^^
     !     symba_tpP => symba_tpP%nextP
     !END DO
     nplplenc = 0
     npltpenc = 0
     IF (symba_plA%helio%swiftest%mass(1) < mtiny) THEN
          nplm = 0
     ELSE
          nplm = 1
     END IF
     irec = 0

     lencounter = ((nplplenc > 0) .OR. (npltpenc > 0))
     IF (lencounter) THEN
          CALL symba_step_interp(lextra_force, lclose, t, npl, nplm, nplmax, ntp, ntpmax, symba_plA, symba_tpA, j2rp2, j4rp4,   &
               dt, eoffset, mtiny, nplplenc, npltpenc, plplenc_list, pltpenc_list, nmergeadd, nmergesub, mergeadd_list,           &
               mergesub_list, encounter_file, out_type)
          lfirst = .TRUE.
     ELSE
          CALL symba_step_helio(lfirst, lextra_force, t, npl, nplm, nplmax, ntp, ntpmax, helio_plA, helio_tpA, j2rp2, j4rp4, dt)
     END IF

     RETURN

END SUBROUTINE symba_step_test
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
