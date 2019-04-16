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
PROGRAM left

! Modules
     USE module_parameters
     USE module_swifter
     USE module_whm
     USE module_random_access
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
     INTEGER(I4B)                                    :: npl,ntp,ntp0,nsp,iout,idump,iloop,i
     REAL(DP)                                        :: t,tfrac,tbase
     CHARACTER(STRMAX)                               :: inparfile
     REAL(DP)                                        :: gmsum,msun,a,e,inc,capom,omega,capm,peri,apo
     LOGICAL(LGT)                                    :: lfirst
     TYPE(whm_pl), DIMENSION(:), ALLOCATABLE, TARGET :: whm_plA
     TYPE(whm_tp), DIMENSION(:), ALLOCATABLE, TARGET :: whm_tpA
     TYPE(whm_pl), POINTER                           :: whm_pl1P
     TYPE(whm_tp), POINTER                           :: whm_tp1P,whm_tpd1P
     TYPE(swifter_pl), POINTER                       :: swifter_pl1P,swifter_plP
     TYPE(swifter_tp), POINTER                       :: swifter_tp1P,swifter_tpd1P,swifter_tpP

! Executable code
     CALL util_version
     WRITE(*,100,ADVANCE="NO")"Enter name of parameter dump file: "
     READ (*,100)inparfile
 100 FORMAT(A)
     inparfile=TRIM(ADJUSTL(inparfile))
     CALL io_init_param(inparfile,nplmax,ntpmax,t0,tstop,dt,inplfile,intpfile,in_type,istep_out,outfile,out_type,out_form,        &
          out_stat,istep_dump,j2rp2,j4rp4,lclose,rmin,rmax,rmaxu,qmin,qmin_coord,qmin_alo,qmin_ahi,encounter_file,lextra_force,   &
          lbig_discard,lrhill_present)
     CALL io_getn(inplfile,intpfile,in_type,npl,nplmax,ntp,ntpmax)
     ALLOCATE(whm_plA(nplmax))
     CALL set_point(whm_plA)
     IF (ntp > 0) THEN
          ALLOCATE(whm_tpA(ntpmax))
          CALL set_point(whm_tpA)
     END IF
     CALL whm_setup(npl,ntp,whm_plA,whm_tpA,whm_pl1P,whm_tp1P,swifter_pl1P,swifter_tp1P)
     CALL io_init_pl(inplfile,in_type,lclose,lrhill_present,npl,swifter_pl1P)
     CALL io_init_tp(intpfile,in_type,ntp,swifter_tp1P)
     CALL util_valid(npl,ntp,swifter_pl1P,swifter_tp1P)

     tfrac = t0/tstop

     OPEN(7,FILE="left_pl.out")
     WRITE(7,*)' '
     IF (t0 < 10000.0_DP) THEN
          WRITE(*,*) 'Time = ', t0
          WRITE(7,*) 'Time = ', t0
     ELSE
          WRITE(*,10) t0
          WRITE(7,10) t0
     END IF
 10  FORMAT(' Time = ', ES13.5)
     WRITE(*,*) 'Fraction done = ', tfrac
     WRITE(7,*) 'Fraction done = ', tfrac

     WRITE(7,*)' '
     WRITE(7,*)'Number of planets (including the Sun): ',npl
     WRITE(7,*)' '
     WRITE(7,*)'Planet:'
     WRITE(7,*)'      id      mass               a               e              i'
     msun=swifter_pl1P%mass
     swifter_plP=>swifter_pl1P
     DO i=2,npl
          swifter_plP=>swifter_plP%nextP
          gmsum=swifter_plP%mass+msun
          CALL orbel_xv2el(swifter_plP%xh(:),swifter_plP%vh(:),gmsum,a,e,inc,capom,omega,capm)
          WRITE(7,101)swifter_plP%id,swifter_plP%mass,a,e,inc*DEGRAD
 101      FORMAT(5x,i4,1x,1p1e13.5,3(5x,0p1f10.4))
     END DO
     CLOSE(7)

     OPEN(7,FILE="left_tp.out")
     WRITE(7,*)' '
     IF (t0 < 10000.0_DP) THEN
          WRITE(7,*) 'Time = ', t0
     ELSE
          WRITE(7,10) t0
     END IF
     WRITE(7,*)'Fraction done = ', tfrac

     WRITE(7,*)' '
     WRITE(*,*)'Number of active test particles: ',ntp
     WRITE(7,*)'Number of active test particles: ',ntp
     WRITE(7,*)' '
     WRITE(7,*)'Test particle:'
     WRITE(7,*)'      id           a              q              Q              i'
     swifter_tpP=>swifter_tp1P
     DO i=1,ntp
          CALL orbel_xv2el(swifter_tpP%xh(:),swifter_tpP%vh(:),msun,a,e,inc,capom,omega,capm)
          apo = a*(1.0_DP + e)
          peri = a*(1.0_DP - e)
          WRITE(7,102)swifter_tpP%id,a,peri,apo,inc*DEGRAD
 102      FORMAT(5x,i4,4(5x,f10.4))
          swifter_tpP=>swifter_tpP%nextP
     ENDDO
     CLOSE(7)

     STOP

END PROGRAM left
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
