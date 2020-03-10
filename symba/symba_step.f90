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
SUBROUTINE symba_step(lfirst, lextra_force, lclose, t, npl, nplmax, ntp, ntpmax, symba_plA, symba_tpA, j2rp2, j4rp4, dt,        &
     nplplenc, npltpenc, plplenc_list, pltpenc_list, nmergeadd, nmergesub, mergeadd_list, mergesub_list, eoffset, mtiny,          &
     encounter_file, out_type, num_plpl_comparisons, ik_plpl, jk_plpl, ik_pltp, jk_pltp)

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
     TYPE(symba_pl), INTENT(INOUT)                    :: symba_plA
     TYPE(symba_tp), INTENT(INOUT)                    :: symba_tpA
     TYPE(symba_plplenc), INTENT(INOUT)               :: plplenc_list
     TYPE(symba_pltpenc), INTENT(INOUT)               :: pltpenc_list
     TYPE(symba_merger), INTENT(INOUT)                :: mergeadd_list, mergesub_list
     INTEGER(I4B), INTENT(IN)                         :: num_plpl_comparisons
     INTEGER(I4B), DIMENSION(num_plpl_comparisons),INTENT(IN)            :: ik_plpl, jk_plpl ! relates the linear index to the i,j matrix indices for planet-planet collisions
     INTEGER(I4B), DIMENSION((npl-1)*ntp),INTENT(IN)  :: ik_pltp, jk_pltp ! linear index to i,j matrix for planet-test particle

! Internals
     LOGICAL(LGT)              :: lencounter, lvdotr
     INTEGER(I4B)              :: i, j, irec, nplm, k
     REAL(DP), DIMENSION(NDIM) :: xr, vr
     REAL(DP), DIMENSION(NDIM,num_plpl_comparisons) :: dist_plpl_array, vel_plpl_array
     REAL(DP), DIMENSION(NDIM,(npl-1)*ntp) :: dist_pltp_array, vel_pltp_array
     LOGICAL(LGT), dimension(num_plpl_comparisons) :: plpl_encounters, plpl_lvdotr ! array for plpl encounters, will record when comparison results in encounter
     LOGICAL(LGT), dimension((npl-1)*ntp) :: pltp_encounters
     INTEGER(I4B), dimension(num_plpl_comparisons) :: plpl_irec
     REAL(DP) :: start, finish
     
! Executable code

     plpl_encounters = .FALSE.
     pltp_encounters = .FALSE.

     ! initialize planets
     symba_plA%nplenc(1:npl) = 0 ! number of planet encounters this particular planet has
     symba_plA%ntpenc(1:npl) = 0 ! number of test particle encounters this particle planet has
     symba_plA%levelg(1:npl) = -1 ! 
     symba_plA%levelm(1:npl) = -1 ! 
     symba_plA%index_parent(1:npl) = (/ (i, i=1,npl)/)
     symba_plA%index_child(:,1:npl) = 0

     ! initialize test particles
     symba_tpA%nplenc(1:ntp) = 0 
     symba_tpA%levelg(1:ntp) = -1
     symba_tpA%levelm(1:ntp) = -1

     nplplenc = 0 ! number of encounters in the entire run 
     npltpenc = 0
     IF (symba_plA%helio%swiftest%mass(1) < mtiny) THEN
          nplm = 0 ! number of planets > mtiny
     ELSE
          nplm = 1
     END IF
     irec = 0 ! recursion counter, 0 since we're in the top loop
     plpl_irec = 0

! ALL THIS NEEDS TO BE CHANGED TO THE TREE SEARCH FUNCTION FOR ENCOUNTERS

     CALL util_dist_eucl_plpl(npl,symba_plA%helio%swiftest%xh, num_plpl_comparisons, ik_plpl, jk_plpl, dist_plpl_array) ! does not care about mtiny
     CALL util_dist_eucl_plpl(npl,symba_plA%helio%swiftest%vh, num_plpl_comparisons, ik_plpl, jk_plpl, vel_plpl_array) ! does not care about mtiny

!$omp parallel do default(none) &
!$omp private(i, lencounter) &
!$omp shared(dist_plpl_array, vel_plpl_array, ik_plpl, jk_plpl, symba_plA, plpl_encounters, plpl_irec, plpl_lvdotr, dt)

     DO i = 1,num_plpl_comparisons
          CALL symba_chk(dist_plpl_array(:,i), vel_plpl_array(:,i), symba_plA%helio%swiftest%rhill(ik_plpl(i)), &
               symba_plA%helio%swiftest%rhill(jk_plpl(i)), dt, plpl_irec(i), lencounter, plpl_lvdotr(i))
          IF (lencounter) THEN
               plpl_encounters(i) = .TRUE. ! record that this comparison resulted in an encounter
               ! for the i particle
               symba_plA%lmerged(ik_plpl(i)) = .FALSE. ! they have not merged YET
               symba_plA%nplenc(ik_plpl(i)) = symba_plA%nplenc(ik_plpl(i)) + 1 ! number of particles that planet "i" has close encountered
               symba_plA%levelg(ik_plpl(i)) = plpl_irec(i) ! recursion level
               symba_plA%levelm(ik_plpl(i)) = plpl_irec(i) ! recursion level
               symba_plA%nchild(ik_plpl(i)) = 0 
               ! for the j particle
               symba_plA%lmerged(jk_plpl(i)) = .FALSE.
               symba_plA%nplenc(jk_plpl(i)) = symba_plA%nplenc(jk_plpl(i)) + 1
               symba_plA%levelg(jk_plpl(i)) = plpl_irec(i)
               symba_plA%levelm(jk_plpl(i)) = plpl_irec(i)
               symba_plA%nchild(jk_plpl(i)) = 0
          END IF
     ENDDO

!$omp end parallel do

     ! here i'll order the encounters
     if(any(plpl_encounters))then ! first check if there were any encounters
          do i = 1,num_plpl_comparisons
               if(plpl_encounters(i))then
                    nplplenc = nplplenc + 1

                    IF (nplplenc > NENMAX) THEN ! there can only be so many recorded planet-planet encounters
                         WRITE(*, *) "SWIFTER Error:"
                         WRITE(*, *) "   PL-PL encounter list is full."
                         WRITE(*, *) "   STOPPING..."
                         CALL util_exit(FAILURE)
                    END IF

                    plplenc_list%status(nplplenc) = ACTIVE ! you are in an encounter
                    plplenc_list%lvdotr(nplplenc) = plpl_lvdotr(i) ! flag of relative accelerations to say if there will be a close encounter in next timestep
                    plplenc_list%level(nplplenc) = plpl_irec(i) ! recursion level
                    plplenc_list%index1(nplplenc) = ik_plpl(i) ! index of first planet in encounter
                    plplenc_list%index2(nplplenc) = jk_plpl(i) ! index of second planet in encounter
               endif
          enddo
     endif

     if(ntp>0)then
          CALL util_dist_eucl_pltp(npl, ntp, symba_plA%helio%swiftest%xh, symba_tpA%helio%swiftest%xh, &
               ik_pltp, jk_pltp, dist_pltp_array)
          CALL util_dist_eucl_pltp(npl, ntp, symba_plA%helio%swiftest%vh, symba_tpA%helio%swiftest%vh, &
               ik_pltp, jk_pltp, vel_pltp_array)

!$omp parallel do

          DO i = 1,(npl-1)*ntp
               CALL symba_chk(dist_pltp_array(:,i), vel_pltp_array(:,i), &
                    symba_plA%helio%swiftest%rhill(ik_pltp(i)), 0.0_DP, dt, irec, lencounter, lvdotr)
               IF (lencounter) THEN
                    ! for the planet
                    symba_plA%ntpenc(ik_pltp(i)) = symba_plA%ntpenc(ik_pltp(i)) + 1
                    symba_plA%levelg(ik_pltp(i)) = irec
                    symba_plA%levelm(ik_pltp(i)) = irec
                    ! for the test particle
                    symba_tpA%nplenc(jk_pltp(i)) = symba_tpA%nplenc(jk_pltp(i)) + 1
                    symba_tpA%levelg(jk_pltp(i)) = irec
                    symba_tpA%levelm(jk_pltp(i)) = irec
                    pltp_encounters(i) = .TRUE.
               END IF
          ENDDO

!$omp end parallel do

          if(any(pltp_encounters))then
               do i = 1,(npl-1)*ntp
                    if(pltp_encounters(i))then

                         npltpenc = npltpenc + 1 ! increment number of planet-test particle interactions
                         IF (npltpenc > NENMAX) THEN
                              WRITE(*, *) "SWIFTER Error:"
                              WRITE(*, *) "   PL-TP encounter list is full."
                              WRITE(*, *) "   STOPPING..."
                              CALL util_exit(FAILURE)
                         END IF
                         pltpenc_list%status(npltpenc) = ACTIVE
                         pltpenc_list%lvdotr(npltpenc) = lvdotr
                         pltpenc_list%level(npltpenc) = irec
                         pltpenc_list%indexpl(npltpenc) = ik_pltp(i)
                         pltpenc_list%indextp(npltpenc) = jk_pltp(i)
                    endif
               enddo
          endif
     endif


     ! ! double loop through all planet-planet and then planet-test particle interactions
     ! DO i = 2, npl ! start at i=2, since we don't need the central body
     !      IF (symba_plA%helio%swiftest%mass(i) < mtiny) EXIT ! don't do this if it's below mtiny
     !      nplm = nplm + 1 ! otherwise, increase the count of planets > mtiny
     !      DO j = i + 1, npl ! through all planets
     !           ! xr = relative distance, vr = relative velocity
     !           ! each of xh and vh arrays is (3,npl), where each row is ~(x, y, z) in heliocentric coordinates
     !           ! xr and vr are thus xr(3)
     !           xr(:) = symba_plA%helio%swiftest%xh(:,j) - symba_plA%helio%swiftest%xh(:,i)
     !           vr(:) = symba_plA%helio%swiftest%vh(:,j) - symba_plA%helio%swiftest%vh(:,i)
     !           ! compare these two particles in symba_chk to check for close encounters
     !           CALL symba_chk(xr(:), vr(:), symba_plA%helio%swiftest%rhill(i), &
     !                symba_plA%helio%swiftest%rhill(j), dt, irec, lencounter, lvdotr)
     !           ! if symba_chk returns a positive lencounter flag, the two planets are in a close encounter
     !           ! OR there may be a close encounter in next timestep (by using vr and the timestep)
     !           IF (lencounter) THEN
     !                nplplenc = nplplenc + 1 ! increment the total number of encounters for this timestep
     !                IF (nplplenc > NENMAX) THEN ! there can only be so many recorded planet-planet encounters
     !                     WRITE(*, *) "SWIFTER Error:"
     !                     WRITE(*, *) "   PL-PL encounter list is full."
     !                     WRITE(*, *) "   STOPPING..."
     !                     CALL util_exit(FAILURE)
     !                END IF
     !                plplenc_list%status(nplplenc) = ACTIVE ! you are in an encounter
     !                plplenc_list%lvdotr(nplplenc) = lvdotr ! flag of relative accelerations to say if there will be a close encounter in next timestep
     !                plplenc_list%level(nplplenc) = irec ! recursion level
     !                plplenc_list%index1(nplplenc) = i ! index of first planet in encounter
     !                plplenc_list%index2(nplplenc) = j ! index of second planet in encounter
     !                ! for the i particle
     !                symba_plA%lmerged(i) = .FALSE. ! they have not merged YET
     !                symba_plA%nplenc(i) = symba_plA%nplenc(i) + 1 ! number of particles that planet "i" has close encountered
     !                symba_plA%levelg(i) = irec ! recursion level
     !                symba_plA%levelm(i) = irec ! recursion level
     !                symba_plA%nchild(i) = 0 
     !                ! for the j particle
     !                symba_plA%lmerged(j) = .FALSE.
     !                symba_plA%nplenc(j) = symba_plA%nplenc(j) + 1
     !                symba_plA%levelg(j) = irec
     !                symba_plA%levelm(j) = irec
     !                symba_plA%nchild(j) = 0
     !           END IF
     !      END DO
     !      DO j = 1, ntp ! through all test particles, same as above
     !           xr(:) = symba_tpA%helio%swiftest%xh(:,j) - symba_plA%helio%swiftest%xh(:,i)
     !           vr(:) = symba_tpA%helio%swiftest%vh(:,j) - symba_plA%helio%swiftest%vh(:,i)
     !           CALL symba_chk(xr(:), vr(:), symba_plA%helio%swiftest%rhill(i), 0.0_DP, dt, irec, lencounter, lvdotr)
     !           IF (lencounter) THEN
     !                ! for the planet
     !                symba_plA%ntpenc(i) = symba_plA%ntpenc(i) + 1
     !                symba_plA%levelg(i) = irec
     !                symba_plA%levelm(i) = irec
     !                ! for the test particle
     !                symba_tpA%nplenc(j) = symba_tpA%nplenc(j) + 1
     !                symba_tpA%levelg(j) = irec
     !                symba_tpA%levelm(j) = irec
     !                npltpenc = npltpenc + 1 ! increment number of planet-test particle interactions
     !                IF (npltpenc > NENMAX) THEN
     !                     WRITE(*, *) "SWIFTER Error:"
     !                     WRITE(*, *) "   PL-TP encounter list is full."
     !                     WRITE(*, *) "   STOPPING..."
     !                     CALL util_exit(FAILURE)
     !                END IF
     !                pltpenc_list%status(npltpenc) = ACTIVE
     !                pltpenc_list%lvdotr(npltpenc) = lvdotr
     !                pltpenc_list%level(npltpenc) = irec
     !                pltpenc_list%indexpl(npltpenc) = i
     !                pltpenc_list%indextp(npltpenc) = j
     !           END IF
     !      END DO
     ! END DO

! END OF THINGS THAT NEED TO BE CHANGED IN THE TREE

     ! temp
     nplm = npl

     ! flag to see if there was an encounter
     lencounter = ((nplplenc > 0) .OR. (npltpenc > 0))
     IF (lencounter) THEN ! if there was an encounter, we need to enter symba_step_interp to see if we need recursion
          CALL symba_step_interp(lextra_force, lclose, t, npl, nplm, nplmax, ntp, ntpmax, symba_plA, symba_tpA, j2rp2, j4rp4,   &
               dt, eoffset, mtiny, nplplenc, npltpenc, plplenc_list, pltpenc_list, nmergeadd, nmergesub, mergeadd_list,           &
               mergesub_list, encounter_file, out_type)
          lfirst = .TRUE.
     ELSE ! otherwise we can just advance the particles
          CALL symba_step_helio(lfirst, lextra_force, t, npl, nplm, nplmax, ntp, ntpmax, symba_plA%helio, symba_tpA%helio, &
               j2rp2, j4rp4, dt)
     END IF

     RETURN

END SUBROUTINE symba_step
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
