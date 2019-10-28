!**********************************************************************************************************************************
!
!  Unit Name   : symba_add
!  Unit Type   : subroutine
!  Project     : Swiftest
!  Package     : symba
!  Language    : Fortran 90/95
!
!  Description : Set up pointers within SyMBA, helio and Swifter planet and test particle structure linked-lists
!
!  Input
!    Arguments : npl          : number of planets
!                ntp          : number of active test particles
!                symba_plA    : SyMBA planet structure array
!                symba_tpA    : SyMBA test particle structure array
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : symba_plA    : SyMBA planet structure array
!                symba_tpA    : SyMBA test particle structure array
!                symba_pl1P   : pointer to head of SyMBA planet structure linked-list
!                symba_tp1P   : pointer to head of active SyMBA test particle structure linked-list
!                swifter_pl1P : pointer to head of Swifter planet structure linked-list
!                swifter_tp1P : pointer to head of active Swifter test particle structure linked-list
!    Terminal  : none
!    File      : none
!
!  Invocation  : CALL symba_add(npl, ntp, symba_plA, symba_tpA, symba_pl1P, symba_tp1P, swifter_pl1P, swifter_tp1P)
!
!  Notes       : 
!
!**********************************************************************************************************************************
SUBROUTINE symba_add(npl, mergeadd_list, nmergeadd, symba_pl1P, swifter_pl1P, mtiny)

! Modules
     USE module_parameters
     USE module_swifter
     USE module_helio
     USE module_symba
     USE module_interfaces, EXCEPT_THIS_ONE => symba_add
     IMPLICIT NONE

! Arguments
     INTEGER(I4B), INTENT(INOUT)                         :: npl
     INTEGER(I4B), INTENT(IN)                            :: nmergeadd
     REAL(DP), INTENT(IN)                                :: mtiny
     TYPE(swifter_pl), POINTER                           :: swifter_pl1P
     TYPE(symba_merger), INTENT(IN)                      :: mergeadd_list
     TYPE(symba_pl), POINTER                             :: symba_pl1P
! Internals
     INTEGER(I4B)              :: i, nplm, add_nplm
     TYPE(swifter_pl), POINTER :: swifter_pliP, swifter_pljP
     TYPE(swifter_pl), POINTER :: add_swifter_pl1P, swifter_plP, end_swifter_pl1P
     TYPE(helio_pl), POINTER   :: add_helio_pl1P, helio_plP, end_helio_pl1P
     TYPE(symba_pl), POINTER   :: add_symba_pl1P, symba_plP, end_symba_pl1P, end_symba_plP
     TYPE(symba_pl), POINTER   :: symba_pliP, symba_pljP ! copied from symba_step
     TYPE(swifter_pl), POINTER   :: swifter_plinit
     TYPE(symba_pl), DIMENSION(:), ALLOCATABLE, TARGET :: fragadd_list
! Executable code
! X 0 . Allocate new structure array of type (symba_pl)
! X 1. pull in the array of bodies to be added to the system in this timestep
! X 2. sort the add_array by decreasing order of mass
!   2.bis call util_hills(nmergeadd,add_swifter_pl1P)
! X 3. get the index of mtiny from nplm
! X 4. see how many particles in the new add_array fall above and below mtiny
! ? 5. break the existing linked list at mtiny+1
! X 6. insert add_array into existing linked list at mtiny+1
! X 7. restructure the linked list
! X 8. find new index of mtiny and update index

!! Allocation of new structure array of type (symba_pl)
     ALLOCATE(fragadd_list(nmergeadd))

! **************************** THIS IS COPIED FROM SYMBA_SETUP AND USED TO SET UP ALL THE POINTERS
! This create a linked list of the particles we are adding at this timestep from the array of particles we are adding at this timestep
     add_symba_pl1P => fragadd_list(1)
     add_helio_pl1P => fragadd_list(1)%helio
     add_swifter_pl1P => fragadd_list(1)%helio%swifter
     NULLIFY(add_symba_pl1P%prevP)
     NULLIFY(add_helio_pl1P%prevP)
     NULLIFY(add_swifter_pl1P%prevP)
     IF (nmergeadd == 1) THEN
          NULLIFY(add_symba_pl1P%nextP)
          NULLIFY(add_helio_pl1P%nextP)
          NULLIFY(add_swifter_pl1P%nextP)
     ELSE
          add_symba_pl1P%nextP => fragadd_list(2)
          add_helio_pl1P%nextP => fragadd_list(2)%helio
          add_swifter_pl1P%nextP => fragadd_list(2)%helio%swifter
          DO i = 2, nmergeadd - 1
               fragadd_list(i)%prevP => fragadd_list(i-1)
               fragadd_list(i)%nextP => fragadd_list(i+1)
               helio_plP => fragadd_list(i)%helio
               helio_plP%prevP => fragadd_list(i-1)%helio
               helio_plP%nextP => fragadd_list(i+1)%helio
               swifter_plP => fragadd_list(i)%helio%swifter
               swifter_plP%prevP => fragadd_list(i-1)%helio%swifter
               swifter_plP%nextP => fragadd_list(i+1)%helio%swifter
          END DO
          fragadd_list(nmergeadd)%prevP => fragadd_list(nmergeadd-1)
          symba_plP => fragadd_list(nmergeadd)
          NULLIFY(symba_plP%nextP)
          helio_plP => fragadd_list(nmergeadd)%helio
          helio_plP%prevP => fragadd_list(nmergeadd-1)%helio
          NULLIFY(helio_plP%nextP)
          swifter_plP => fragadd_list(nmergeadd)%helio%swifter
          swifter_plP%prevP => fragadd_list(nmergeadd-1)%helio%swifter
          NULLIFY(swifter_plP%nextP)
     END IF
! initialize fragadd_list pointers 
     swifter_plinit => add_swifter_pl1P
     swifter_plinit%xh = mergeadd_list(i)%xh
     swifter_plinit%vh = mergeadd_list(i)%vh
     swifter_plinit%mass = mergeadd_list(i)%mass
     swifter_plinit%radius = mergeadd_list(i)%radius
     swifter_plinit%id = mergeadd_list(i)%id
     swifter_plinit%status = mergeadd_list(i)%status
     DO i=2,nmergeadd-1
          swifter_plinit%xh = mergeadd_list(i)%xh
          swifter_plinit%vh = mergeadd_list(i)%vh
          swifter_plinit%mass = mergeadd_list(i)%mass
          swifter_plinit%radius = mergeadd_list(i)%radius
          swifter_plinit%id = mergeadd_list(i)%id
          swifter_plinit%status = mergeadd_list(i)%status
          swifter_plinit => swifter_plinit%nextP
          !PRINT *, 'mergeadd_list mass', fragadd_list(i)%helio%swifter%mass
     END DO



! calculate hill radius inside add_swifter_pl1P particles
     CALL util_hills(nmergeadd,add_swifter_pl1P)
! This reorders the linked list of particles to add at this timestep
     CALL symba_reorder_pl(nmergeadd,add_symba_pl1P)

! **************************** THIS IS COPIED FROM SYMBA_STEP AND USED TO UPDATE NPLM
! This searched the total particle list from the previous timestep (WITHOUT the new additions) to find the location of mtiny (nplm)
     IF (symba_pl1P%helio%swifter%mass < mtiny) THEN
          nplm = 0
     ELSE
          nplm = 1
     END IF
     symba_pliP => symba_pl1P
     DO i = 2, npl
          symba_pliP => symba_pliP%nextP
          swifter_pliP => symba_pliP%helio%swifter
          IF (swifter_pliP%mass < mtiny) EXIT
          nplm = nplm + 1
     END DO

! This searched the new list of particles to add to find the location of mtiny (add_nplm)
     IF (add_symba_pl1P%helio%swifter%mass < mtiny) THEN
          add_nplm = 0
     ELSE
          add_nplm = 1
     END IF
     symba_pliP => add_symba_pl1P
     DO i = 2, nmergeadd
          symba_pliP => symba_pliP%nextP
          swifter_pliP => symba_pliP%helio%swifter
          IF (swifter_pliP%mass < mtiny) EXIT
          add_nplm = add_nplm + 1
     END DO

! **************************** THIS IS COPIED FROM SYMBA_DISCARD_SPILL_PL AND USED TO INSERT OUR ADDED LINKED LIST TO THE EXISTING LINKED LIST
! We might have to actually break the old list before we stick stuff in it, not sure. CHECK ON THIS.
! This section inserts the head of the new list into the old list at mtiny 
     IF ((nmergeadd > 0) .AND. (nplm > 1) .AND. (nplm < npl)) THEN
          symba_plP => symba_pl1P
          end_symba_plP => symba_pl1P
          DO i = 1, nplm-1
               symba_plP => symba_plP%nextP
               end_symba_plP => end_symba_plP%nextP
          END DO
          end_symba_plP => end_symba_plP%nextP
     ! This takes the head of our new list and points the previous to the index of mtiny (nplm) in the old list
          add_symba_pl1P%prevP => symba_plP
          add_helio_pl1P%prevP => helio_plP
          add_swifter_pl1P%prevP => swifter_plP
     ! This takes the index of mtiny (nplm) in the old list and points it to the head of the new list
          symba_plP%nextP => add_symba_pl1P
          symba_plP%helio%nextP => add_helio_pl1P
          symba_plP%helio%swifter%nextP => add_swifter_pl1P
     ! This section takes the add list and goes to the end of it
          DO i = 1, nmergeadd-1
               symba_plP => symba_plP%nextP
          END DO 
     ! This takes the first of the particles below mtiny (nplm) of our old list and points to its previous to the last of the add list 
          end_symba_pl1P%prevP => symba_plP
          end_helio_pl1P%prevP => helio_plP
          end_swifter_pl1P%prevP => swifter_plP
     ! This takes the last of the add list and points its next to the first of the particles below mtiny (nplm) in the old list 
          symba_plP%nextP => end_symba_pl1P
          symba_plP%helio%nextP => end_helio_pl1P
          symba_plP%helio%swifter%nextP => end_swifter_pl1P
     ! This section goes through the rest of the old list up until the second to last particle
          DO i = 1, npl-nplm-1
               symba_plP => symba_plP%nextP
          END DO 
     ! This section removed the next of the very last particle
          NULLIFY(symba_plP%nextP)
          NULLIFY(helio_plP%nextP)
          NULLIFY(swifter_plP%nextP)
     END IF
! This does the same if all the bodies in the system are massive
     IF (nplm == npl) THEN 
          symba_plP => symba_pl1P
          DO i = 1, nplm - 1
               symba_plP => symba_plP%nextP
          END DO
     ! This takes the head of our new list and points the previous to the index of mtiny (nplm) in the old list
          add_symba_pl1P%prevP => symba_plP
          add_helio_pl1P%prevP => helio_plP
          add_swifter_pl1P%prevP => swifter_plP
     ! This takes the index of mtiny (nplm) in the old list and points it to the head of the new list
          symba_plP%nextP => add_symba_pl1P
          symba_plP%helio%nextP => add_helio_pl1P
          symba_plP%helio%swifter%nextP => add_swifter_pl1P
     ! This section takes the add list and goes to the end of it
          DO i = 1, (nmergeadd-1)
               symba_plP => symba_plP%nextP
          END DO 
     ! This section removed the next of the very last particle
          NULLIFY(symba_plP%nextP)
          NULLIFY(helio_plP%nextP)
          NULLIFY(swifter_plP%nextP)
     END IF
! This does the same if the sun is the only massive body in the system
     IF (nplm == 1) THEN 
          symba_plP => symba_pl1P
          end_symba_plP => symba_pl1P
          end_symba_plP => end_symba_plP%nextP
     ! This takes the head of our new list and points the previous to the index of mtiny (nplm) in the old list
          add_symba_pl1P%prevP => symba_plP
          add_helio_pl1P%prevP => helio_plP
          add_swifter_pl1P%prevP => swifter_plP
     ! This takes the index of mtiny (nplm) in the old list and points it to the head of the new list
          symba_plP%nextP => add_symba_pl1P
          symba_plP%helio%nextP => add_helio_pl1P
          symba_plP%helio%swifter%nextP => add_swifter_pl1P
     ! This section takes the add list and goes to the end of it
          DO i = 1, (nmergeadd-1)
               symba_plP => symba_plP%nextP
          END DO 
     ! This takes the first of the particles below mtiny (nplm) of our old list and points to its previous to the last of the add list 
          end_symba_pl1P%prevP => symba_plP
          end_helio_pl1P%prevP => helio_plP
          end_swifter_pl1P%prevP => swifter_plP
     ! This takes the last of the add list and points its next to the first of the particles below mtiny (nplm) in the old list 
          symba_plP%nextP => end_symba_pl1P
          symba_plP%helio%nextP => end_helio_pl1P
          symba_plP%helio%swifter%nextP => end_swifter_pl1P
     ! This section goes through the rest of the old list up until the second to last particle
          DO i = 1, (npl-nplm-1)
               symba_plP => symba_plP%nextP
          END DO 
     ! This section removed the next of the very last particle
          NULLIFY(symba_plP%nextP)
          NULLIFY(helio_plP%nextP)
          NULLIFY(swifter_plP%nextP)
     END IF
     
     npl = npl + nmergeadd

     RETURN

END SUBROUTINE symba_add
!**********************************************************************************************************************************
!
!  Author(s)   : Carlisle Wishard and Jennifer Pouplin (Based on symba_setup by David E. Kaufmann)
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
