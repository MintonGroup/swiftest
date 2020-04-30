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
     encounter_file, out_type, num_plpl_comparisons, k_plpl, num_pltp_comparisons, k_pltp, dist_plpl_array, vel_plpl_array, &
     dist_pltp_array, vel_pltp_array)

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
     INTEGER(I4B), INTENT(IN)                         :: num_plpl_comparisons, num_pltp_comparisons
     INTEGER(I4B), DIMENSION(2,num_plpl_comparisons), INTENT(IN) :: k_plpl
     INTEGER(I4B), DIMENSION(2,num_pltp_comparisons), INTENT(IN)  :: k_pltp
     REAL(DP), DIMENSION(NDIM,num_plpl_comparisons) :: dist_plpl_array, vel_plpl_array
     REAL(DP), DIMENSION(NDIM,num_pltp_comparisons) :: dist_pltp_array, vel_pltp_array

! Internals
     LOGICAL(LGT)              :: lencounter, lvdotr
     INTEGER(I4B)              :: i, j, irec, nplm, k, counter
     INTEGER(I4B), ALLOCATABLE :: plpl_encounters_indices(:), pltp_encounters_indices(:)
     REAL(DP), DIMENSION(NDIM) :: xr, vr
     
     ! INTEGER(I4B), dimension(num_pltp_comparisons) :: pltp_encounters, pltp_lvdotr
     ! INTEGER(I4B), dimension(num_plpl_comparisons) :: plpl_encounters, plpl_lvdotr

     ! REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: dist_plpl_array, vel_plpl_array
     ! REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: dist_pltp_array, vel_pltp_array
     INTEGER(I4B), ALLOCATABLE, DIMENSION(:) :: pltp_encounters, pltp_lvdotr
     INTEGER(I4B), ALLOCATABLE, DIMENSION(:) :: plpl_encounters, plpl_lvdotr
     
! Executable code


     allocate(pltp_encounters(num_pltp_comparisons))
     allocate(pltp_lvdotr(num_pltp_comparisons))
     allocate(plpl_encounters(num_plpl_comparisons))
     allocate(plpl_lvdotr(num_plpl_comparisons))

     plpl_encounters = 0
     plpl_lvdotr = 0

     pltp_encounters = 0
     pltp_lvdotr = 0

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

! ALL THIS NEEDS TO BE CHANGED TO THE TREE SEARCH FUNCTION FOR ENCOUNTERS

     CALL util_dist_eucl_plpl(npl,symba_plA%helio%swiftest%xh, num_plpl_comparisons, k_plpl, dist_plpl_array) 
     CALL util_dist_eucl_plpl(npl,symba_plA%helio%swiftest%vh, num_plpl_comparisons, k_plpl, vel_plpl_array) 
     CALL symba_chk_eucl(num_plpl_comparisons, k_plpl, dist_plpl_array, vel_plpl_array, &
          symba_plA%helio%swiftest%rhill, dt, plpl_encounters, plpl_lvdotr)

     ! here i'll order the encounters
     nplplenc = count(plpl_encounters > 0)
     if(nplplenc>0)then

          allocate(plpl_encounters_indices(nplplenc))

          ! plpl_encounters_indices = pack(plpl_encounters,plpl_encounters > 0)
          ! so it turns out this is significantly faster than the pack command
          counter = 1
          do k = 1,num_plpl_comparisons
               if(plpl_encounters(k).gt.0)then
                    plpl_encounters_indices(counter) = k
                    counter = counter + 1
               endif
          enddo

          symba_plA%lmerged(k_plpl(1,plpl_encounters_indices)) = .FALSE. ! they have not merged YET
          symba_plA%nplenc(k_plpl(1,plpl_encounters_indices)) = symba_plA%nplenc(k_plpl(1,plpl_encounters_indices)) + 1 ! number of particles that planet "i" has close encountered
          symba_plA%levelg(k_plpl(1,plpl_encounters_indices)) = 0 ! recursion level
          symba_plA%levelm(k_plpl(1,plpl_encounters_indices)) = 0 ! recursion level
          symba_plA%nchild(k_plpl(1,plpl_encounters_indices)) = 0 
          ! for the j particle
          symba_plA%lmerged(k_plpl(2,plpl_encounters_indices)) = .FALSE.
          symba_plA%nplenc(k_plpl(2,plpl_encounters_indices)) = symba_plA%nplenc(k_plpl(2,plpl_encounters_indices)) + 1
          symba_plA%levelg(k_plpl(2,plpl_encounters_indices)) = 0
          symba_plA%levelm(k_plpl(2,plpl_encounters_indices)) = 0
          symba_plA%nchild(k_plpl(2,plpl_encounters_indices)) = 0

          plplenc_list%status(1:nplplenc) = ACTIVE ! you are in an encounter
          plplenc_list%lvdotr(1:nplplenc) = plpl_lvdotr(plpl_encounters_indices)! flag of relative accelerations to say if there will be a close encounter in next timestep 
          plplenc_list%level(1:nplplenc)  = 0 ! recursion level
          plplenc_list%index1(1:nplplenc) = k_plpl(1,plpl_encounters_indices) ! index of first planet in encounter
          plplenc_list%index2(1:nplplenc) = k_plpl(2,plpl_encounters_indices) ! index of second planet in encounter
          deallocate(plpl_encounters_indices)
     endif

     if(ntp>0)then

          ! CALL util_dist_eucl_pltp(npl, ntp, symba_plA%helio%swiftest%xh, symba_tpA%helio%swiftest%xh, &
          !      num_pltp_comparisons, k_pltp, dist_pltp_array)
          ! CALL util_dist_eucl_pltp(npl, ntp, symba_plA%helio%swiftest%vh, symba_tpA%helio%swiftest%vh, &
          !      num_pltp_comparisons, k_pltp, vel_pltp_array)
          CALL symba_chk_eucl_pltp(num_pltp_comparisons, k_pltp, symba_plA, symba_tpA, dt, pltp_encounters, pltp_lvdotr)
     
          npltpenc = count(pltp_encounters > 0)
          if(npltpenc>0)then

               allocate(pltp_encounters_indices(npltpenc))

               counter = 1
               do k = 1,num_pltp_comparisons
                    if(pltp_encounters(k).gt.0)then
                         pltp_encounters_indices(counter) = k
                         counter = counter + 1
                    endif
               enddo

               symba_plA%ntpenc(k_pltp(1,pltp_encounters_indices)) = symba_plA%ntpenc(k_pltp(1,pltp_encounters_indices)) + 1
               symba_plA%levelg(k_pltp(1,pltp_encounters_indices)) = 0
               symba_plA%levelm(k_pltp(1,pltp_encounters_indices)) = 0

               symba_tpA%nplenc(k_pltp(2,pltp_encounters_indices)) = symba_tpA%nplenc(k_pltp(2,pltp_encounters_indices)) + 1
               symba_tpA%levelg(k_pltp(2,pltp_encounters_indices)) = 0
               symba_tpA%levelm(k_pltp(2,pltp_encounters_indices)) = 0

               pltpenc_list%status(1:npltpenc) = ACTIVE
               pltpenc_list%lvdotr(1:npltpenc) = pltp_lvdotr(pltp_encounters_indices)
               pltpenc_list%level(1:npltpenc) = 0
               pltpenc_list%indexpl(1:npltpenc) = k_pltp(1,pltp_encounters_indices)
               pltpenc_list%indextp(1:npltpenc) = k_pltp(1,pltp_encounters_indices)

               deallocate(pltp_encounters_indices)
          endif
     endif
     
! END OF THINGS THAT NEED TO BE CHANGED IN THE TREE

     nplm = count(symba_plA%helio%swiftest%mass > mtiny)
     ! flag to see if there was an encounter
     lencounter = ((nplplenc > 0) .OR. (npltpenc > 0))
     IF (lencounter) THEN ! if there was an encounter, we need to enter symba_step_interp to see if we need recursion
          CALL symba_step_interp(lextra_force, lclose, t, npl, nplm, nplmax, ntp, ntpmax, symba_plA, symba_tpA, j2rp2, j4rp4,   &
               dt, eoffset, mtiny, nplplenc, npltpenc, plplenc_list, pltpenc_list, nmergeadd, nmergesub, mergeadd_list,           &
               mergesub_list, encounter_file, out_type, num_plpl_comparisons, k_plpl, num_pltp_comparisons, k_pltp, &
               dist_plpl_array, dist_pltp_array)
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
