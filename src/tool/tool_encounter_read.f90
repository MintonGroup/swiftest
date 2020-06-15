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
     use swiftest
   use io
     IMPLICIT NONE

! Arguments
   type(swiftest_configuration)  :: param    ! derived type containing user-defined parameters
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
     REAL(DP)          :: mtiny          ! Mass cutoff
     CHARACTER(STRMAX) :: qmin_coord     ! Coordinate frame to use for qmin
     CHARACTER(STRMAX) :: encounter_file ! Name of output file for encounters
     CHARACTER(STRMAX) :: inplfile       ! Name of input file for planets
     CHARACTER(STRMAX) :: intpfile       ! Name of input file for test particles
     CHARACTER(STRMAX) :: in_type        ! Format of input data files
     CHARACTER(STRMAX) :: outfile        ! Name of output binary file
     CHARACTER(STRMAX) :: out_type       ! Binary format of output file
     CHARACTER(STRMAX) :: out_form       ! Data to write to output file
     CHARACTER(STRMAX) :: out_stat       ! Open status for output binary file

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
   call config%read_from_file(inparfile)

   ! temporary until the conversion to the derived type argument list is complete
   nplmax = config%nplmax
   ntpmax = config%ntpmax
   t0 = config%t0
   tstop = config%tstop
   dt = config%dt
   inplfile = config%inplfile
   intpfile = config%intpfile
   in_type = config%in_type
   istep_out = config%istep_out
   outfile = config%outfile
   out_type = config%out_type
   out_form = config%out_form
   out_stat = config%out_stat
   istep_dump = config%istep_dump
   j2rp2 = config%j2rp2
   j4rp4 = config%j4rp4
   rmin = config%rmin
   rmax = config%rmax
   rmaxu = config%rmaxu
   qmin = config%qmin
   qmin_coord = config%qmin_coord
   qmin_alo = config%qmin_alo
   qmin_ahi = config%qmin_ahi
   encounter_file = config%encounter_file
   mtiny = config%mtiny
   !^^^^^^^^^^^^^^^^^^^^^^^^^
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
