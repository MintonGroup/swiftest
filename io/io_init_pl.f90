!**********************************************************************************************************************************
!
!  Unit Name   : io_init_pl
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : io
!  Language    : Fortran 90/95
!
!  Description : Read in planet data
!
!  Input
!    Arguments : inplfile       : name of input file for planets
!                in_type        : format of input data file
!                lclose         : logical flag indicating whether planets' physical radii are present in input file
!                lrhill_present : logical flag indicating whether planets' Hill's sphere radii are present in input file
!                lrotation      : logical flag indicating whether planets' rotation parameters are present in the file
!                npl            : number of planets
!                swifter_pl1P   : pointer to head of Swifter planet structure linked-list
!    Terminal  : none
!    File      : id             : planet identifier               (all planets)
!                mass           : mass                            (all planets)
!                rhill          : Hill's sphere radius (optional) (if present, all planets except the central body)
!                radius         : physical radius (optional)      (if present, all planets including the central body)
!                xh             : heliocentric position           (all planets)
!                vh             : heliocentric velocity           (all planets)
!                Ip             : principal moments of inertia    (if present, all planets including the central body)
!                rot            : rotation vector                 (if present, all planets including the central body)
!
!  Output
!    Arguments : swifter_pl1P   : pointer to head of Swifter planet structure linked-list
!    Terminal  : error message
!    File      : none
!
!  Invocation  : CALL io_init_pl(inplfile, in_type, lclose, lrhill_present, npl, swifter_pl1P)
!
!  Notes       : Adapted from Martin Duncan's Swift routine io_init_pl.f
!                xdr file support removed (D. Minton)
!
!**********************************************************************************************************************************
SUBROUTINE io_init_pl(inplfile, in_type, lclose, lrhill_present, lrotation, npl, swifter_pl1P)

! Modules
     USE module_parameters
     USE module_swifter
     USE module_interfaces, EXCEPT_THIS_ONE => io_init_pl
     IMPLICIT NONE

! Arguments
     LOGICAL(LGT), INTENT(IN)  :: lclose, lrhill_present, lrotation
     INTEGER(I4B), INTENT(IN)  :: npl
     CHARACTER(*), INTENT(IN)  :: inplfile, in_type
     TYPE(swifter_pl), POINTER :: swifter_pl1P

! Internals
     INTEGER(I4B), PARAMETER   :: LUN = 7
     INTEGER(I4B)              :: i, iu, ierr, inpl
     TYPE(swifter_pl), POINTER :: swifter_plP

! Executable code
     swifter_plP => swifter_pl1P
     IF (in_type == "ASCII") THEN
          CALL io_open(LUN, inplfile, "OLD", "FORMATTED", ierr)
          READ(LUN, *) inpl
          READ(LUN, *) swifter_plP%id, swifter_plP%mass
          swifter_plP%rhill = 0.0_DP  
          IF (lclose) THEN
              READ(LUN, *) swifter_plP%radius,swifter_plP%k2,swifter_plP%Q
          ELSE
              swifter_plP%radius = 0.0_DP
          END IF
          IF (lrotation) THEN
              READ(LUN, *) swifter_plP%Ip(:)
              READ(LUN, *) swifter_plP%rot(:)
          END IF
          READ(LUN, *) swifter_plP%xh(:)
          READ(LUN, *) swifter_plP%vh(:)
          DO i = 1, NDIM
               IF ((swifter_plP%xh(i) /= 0.0_DP) .OR. (swifter_plP%vh(i) /= 0.0_DP)) THEN
                    WRITE(*, *) "SWIFTER Error:"
                    WRITE(*, *) " Input MUST be in heliocentric coordinates."
                    WRITE(*, *) " Position/velocity components of Body 1 are"
                    WRITE(*, *) swifter_plP%xh(:)
                    WRITE(*, *) swifter_plP%vh(:)
                    CALL util_exit(FAILURE)
               END IF
          END DO
          swifter_plP%status = ACTIVE
          DO i = 2, npl
               swifter_plP => swifter_plP%nextP
               IF (lrhill_present) THEN
                    READ(LUN, *) swifter_plP%id, swifter_plP%mass, swifter_plP%rhill
               ELSE
                    READ(LUN, *) swifter_plP%id, swifter_plP%mass
                    swifter_plP%rhill = 0.0_DP
               END IF
               IF (lclose) THEN
                    READ(LUN, *) swifter_plP%radius,swifter_plP%k2,swifter_plP%Q
               ELSE
                    swifter_plP%radius = 0.0_DP
               END IF
               IF (lrotation) THEN
                   READ(LUN, *) swifter_plP%Ip(:)
                   READ(LUN, *) swifter_plP%rot(:)
               END IF
               READ(LUN, *) swifter_plP%xh(:)
               READ(LUN, *) swifter_plP%vh(:)
               swifter_plP%status = ACTIVE
          END DO
          CLOSE(UNIT = LUN)
     ELSE
          write(*,*) trim(adjustl(in_type))," file type not supported"
          call util_exit(FAILURE)
     END IF

     RETURN

END SUBROUTINE io_init_pl
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
