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
     USE swiftest
     USE module_interfaces
   use io
     IMPLICIT NONE

! Arguments
   type(input_parameters)  :: param    ! derived type containing user-defined parameters
   type(feature_list) :: feature       ! temporary until the parameter derived type conversion is complete
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
   call param%read_from_file(inparfile)

   ! temporary until the conversion to the derived type argument list is complete
   nplmax = param%nplmax
   ntpmax = param%ntpmax
   t0 = param%t0
   tstop = param%tstop
   dt = param%dt
   inplfile = param%inplfile
   intpfile = param%intpfile
   in_type = param%in_type
   istep_out = param%istep_out
   outfile = param%outfile
   out_type = param%out_type
   out_form = param%out_form
   out_stat = param%out_stat
   istep_dump = param%istep_dump
   j2rp2 = param%j2rp2
   j4rp4 = param%j4rp4
   rmin = param%rmin
   rmax = param%rmax
   rmaxu = param%rmaxu
   qmin = param%qmin
   qmin_coord = param%qmin_coord
   qmin_alo = param%qmin_alo
   qmin_ahi = param%qmin_ahi
   encounter_file = param%encounter_file
   mtiny = param%mtiny
   feature = param%feature
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
