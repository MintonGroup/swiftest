!******************************************************************************
!
!  Unit Name   : 
!  Unit Type   : 
!  Project     : SWIFTER
!  Package     :
!  Language    : Fortran 90/95
!
!  Description : 
!
!  Input
!    Arguments :
!    Terminal  :
!    File      :
!
!  Output
!    Arguments :
!    Terminal  :
!    File      :
!
!  Invocation  : 
!
!  Notes       : 
!
!******************************************************************************
PROGRAM tool_encounter_read

! Modules
     USE module_parameters
     USE module_interfaces
     IMPLICIT NONE

! Arguments
     INTEGER(I4B)      :: nplmax         ! Maximum number of planets
     INTEGER(I4B)      :: ntpmax         ! Maximum number of test particles
     INTEGER(I4B)      :: istep_out      ! Time steps between binary outputs
     INTEGER(I4B)      :: istep_dump     ! Time steps between dumps
     REAL(DP)          :: t0             ! Integration start time
     REAL(DP)          :: tstop          ! Integration stop time
     REAL(DP)          :: dt             ! Time step
     REAL(DP)          :: j2rp2          ! J2*R^2 term for central body
     REAL(DP)          :: j4rp4          ! J4*R^4 term for central body
     REAL(DP)          :: rmin           ! Minimum heliocentric radius for test particle
     REAL(DP)          :: rmax           ! Maximum heliocentric radius for test particle
     REAL(DP)          :: rmaxu          ! Maximum unbound heliocentric radius for test particle
     REAL(DP)          :: qmin           ! Minimum pericenter distance for test particle
     REAL(DP)          :: qmin_alo       ! Minimum semimajor axis for qmin
     REAL(DP)          :: qmin_ahi       ! Maximum semimajor axis for qmin
     CHARACTER(STRMAX) :: qmin_coord     ! Coordinate frame to use for qmin
     CHARACTER(STRMAX) :: encounter_file ! Name of output file for encounters
     CHARACTER(STRMAX) :: inplfile       ! Name of input file for planets
     CHARACTER(STRMAX) :: intpfile       ! Name of input file for test particles
     CHARACTER(STRMAX) :: in_type        ! Format of input data files
     CHARACTER(STRMAX) :: outfile        ! Name of output binary file
     CHARACTER(STRMAX) :: out_type       ! Binary format of output file
     CHARACTER(STRMAX) :: out_form       ! Data to write to output file
     CHARACTER(STRMAX) :: out_stat       ! Open status for output binary file
     LOGICAL(LGT)      :: lclose         ! Check for planet-test particle encounters
     LOGICAL(LGT)      :: lextra_force   ! Use user-supplied force routines
     LOGICAL(LGT)      :: lbig_discard   ! Dump planet data with discards
     LOGICAL(LGT)      :: lrhill_present ! Hill's sphere radius present

! Internals
     INTEGER(I4B)              :: i,ierr,id1,id2
     REAL(DP)                  :: t,mass1,mass2
     REAL(DP), DIMENSION(NDIM) :: xh1,xh2,vh1,vh2
     CHARACTER(STRMAX)         :: inparfile

! Executable code
     WRITE(*,100,ADVANCE="NO")"Enter name of parameter data file: "
     READ(*,100)inparfile
 100 FORMAT(A)
     inparfile=TRIM(ADJUSTL(inparfile))
     CALL io_init_param(inparfile,nplmax,ntpmax,t0,tstop,dt,inplfile,intpfile,in_type,istep_out,outfile,out_type,out_form,        &
          out_stat,istep_dump,j2rp2,j4rp4,lclose,rmin,rmax,rmaxu,qmin,qmin_coord,qmin_alo,qmin_ahi,encounter_file,lextra_force,   &
          lbig_discard,lrhill_present)
     ierr=0
     i=0
     DO
          ierr=io_read_encounter(t,id1,id2,mass1,mass2,xh1,xh2,vh1,vh2,encounter_file,out_type)
          IF (ierr /= 0) EXIT
          i=i+1
          WRITE(*,*)"Encounter #",i
          WRITE(*,*)"Time = ",t
          WRITE(*,*)" Body1: ",id1,mass1,xh1,vh1
          WRITE(*,*)" Body2: ",id2,mass2,xh2,vh2
     END DO
     CALL util_exit(SUCCESS)

     STOP

END PROGRAM tool_encounter_read
!******************************************************************************
!
!  Author(s)   : 
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
!******************************************************************************
