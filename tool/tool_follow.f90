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
PROGRAM tool_follow

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
     LOGICAL(LGT)      :: lrotation      ! Rotation parameters present

! Internals
     INTEGER(I4B)                                    :: npl,ntp,ntp0,nsp,iout,iloop,i,iu,ifol,nskp,ierr,iout_form,ic,istep,id
     REAL(DP)                                        :: t,tfrac,tbase,tmax
     CHARACTER(STRMAX)                               :: inparfile
     REAL(DP)                                        :: gmsum,msun,a,e,inc,capom,omega,capm,mass,peri,apo,radius
     REAL(DP)                                        :: xh,yh,zh,vxh,vyh,vzh
     LOGICAL(LGT)                                    :: lfirst
     TYPE(whm_pl), DIMENSION(:), ALLOCATABLE, TARGET :: whm_plA
     TYPE(whm_tp), DIMENSION(:), ALLOCATABLE, TARGET :: whm_tpA
     TYPE(whm_pl), POINTER                           :: whm_pl1P
     TYPE(whm_tp), POINTER                           :: whm_tp1P,whm_tpd1P
     TYPE(swifter_pl), POINTER                       :: swifter_pl1P,swifter_plP
     TYPE(swifter_tp), POINTER                       :: swifter_tp1P,swifter_tpd1P,swifter_tpP

! Executable code
     CALL util_version
     WRITE(*,100,ADVANCE="NO")"Enter name of parameter data file: "
     READ (*,100)inparfile
 100 FORMAT(A)
     inparfile=TRIM(ADJUSTL(inparfile))
     CALL io_init_param(inparfile,nplmax,ntpmax,t0,tstop,dt,inplfile,intpfile,in_type,istep_out,outfile,out_type,out_form,        &
          out_stat,j2rp2,j4rp4,lclose,rmin,rmax,rmaxu,qmin,qmin_coord,qmin_alo,qmin_ahi,encounter_file,lextra_force,   &
          lbig_discard,lrhill_present,lrotation)
     CALL io_getn(inplfile,intpfile,in_type,npl,nplmax,ntp,ntpmax)
     ALLOCATE(whm_plA(nplmax))
     CALL set_point(whm_plA)
     IF (ntp > 0) THEN
          ALLOCATE(whm_tpA(ntpmax))
          CALL set_point(whm_tpA)
     END IF
     CALL whm_setup(npl,ntp,whm_plA,whm_tpA,whm_pl1P,whm_tp1P,swifter_pl1P,swifter_tp1P)
     CALL io_init_pl(inplfile,in_type,lclose,lrhill_present,lrotation,npl,swifter_pl1P)
     CALL io_init_tp(intpfile,in_type,ntp,swifter_tp1P)
     CALL util_valid(npl,ntp,swifter_pl1P,swifter_tp1P)
     iu=20
     WRITE(*,100,ADVANCE="NO")"Enter the id of the particle to follow: "
     READ(*,*)ifol
     WRITE(*,*)"Following particle ",ifol
     WRITE(*,100,ADVANCE="NO")"Enter the frequency: "
     READ(*,*)nskp
     CALL io_open(iu,outfile,"OLD","UNFORMATTED",ierr)
     IF (ierr /= 0) THEN
          WRITE(*,*)"SWIFTER Error:"
          WRITE(*,*)"   Unable to open output binary file"
          CALL util_exit(FAILURE)
     END IF
     SELECT CASE (out_form)
          CASE ("EL")
               iout_form=EL
          CASE ("XV")
               iout_form=XV
          CASE ("FILT")
               iout_form=FILT
     END SELECT
     OPEN(UNIT=7,FILE="follow.out")

     tmax=t0
     ic=0
     DO
          ierr=io_read_hdr(iu,t,npl,ntp,iout_form,out_type)
          IF (ierr /= 0) EXIT
          istep=0
          SELECT CASE (iout_form)
               CASE (EL)
                    DO i=2,npl
                         ierr=io_read_line(iu,id,a,e,inc,capom,omega,capm,out_type,MASS=mass,RADIUS=radius)
                         IF (ierr /= 0) THEN
                              WRITE(*,*)"SWIFTER Error:"
                              WRITE(*,*)"   Stop while reading planets ierr = ",ierr
                              CALL util_exit(FAILURE)
                         END IF
                         IF (id == ifol) THEN
                              istep=1
                              ic=ic+1
                              IF (MOD(ic,nskp) == 0) THEN
                                   inc=inc*DEGRAD
                                   capom=capom*DEGRAD
                                   omega=omega*DEGRAD
                                   capm=capm*DEGRAD
                                   peri=a*(1.0_DP-e)
                                   apo=a*(1.0_DP+e)
                                   WRITE(7,1000)t,ifol,a,e,inc,capom,omega,capm,peri,apo
 1000                              FORMAT(1x,e15.7,1x,i3,1x,f10.4,1x,f7.5,4(1x,f9.4),2(1x,f10.4))
                              END IF
                              tmax=t
                         END IF
                    END DO
                    DO i=1,ntp
                         ierr=io_read_line(iu,id,a,e,inc,capom,omega,capm,out_type)
                         IF (ierr /= 0) THEN
                              WRITE(*,*)"SWIFTER Error:"
                              WRITE(*,*)"   Stop while reading test particles ierr = ",ierr
                              CALL util_exit(FAILURE)
                         END IF
                         IF (id == ifol) THEN
                              istep=1
                              ic=ic+1
                              IF (MOD(ic,nskp) == 0) THEN
                                   inc=inc*DEGRAD
                                   capom=capom*DEGRAD
                                   omega=omega*DEGRAD
                                   capm=capm*DEGRAD
                                   peri=a*(1.0_DP-e)
                                   apo=a*(1.0_DP+e)
                                   WRITE(7,1000)t,ifol,a,e,inc,capom,omega,capm,peri,apo
                              END IF
                              tmax=t
                         END IF
                    END DO
               CASE (XV)
                    DO i=2,npl
                         ierr=io_read_line(iu,id,xh,yh,zh,vxh,vyh,vzh,out_type,MASS=mass,RADIUS=radius)
                         IF (ierr /= 0) THEN
                              WRITE(*,*)"SWIFTER Error:"
                              WRITE(*,*)"   Stop while reading planets ierr = ",ierr
                              CALL util_exit(FAILURE)
                         END IF
                         IF (id == ifol) THEN
                              istep=1
                              ic=ic+1
                              IF (MOD(ic,nskp) == 0) THEN
                                   WRITE(7,1001)t,ifol,xh,yh,zh,vxh,vyh,vzh
 1001                              FORMAT(1x,e15.7,1x,i3,6(1x,es15.7))
                              END IF
                              tmax=t
                         END IF
                    END DO
                    DO i=1,ntp
                         ierr=io_read_line(iu,id,xh,yh,zh,vxh,vyh,vzh,out_type)
                         IF (ierr /= 0) THEN
                              WRITE(*,*)"SWIFTER Error:"
                              WRITE(*,*)"   Stop while reading test particles ierr = ",ierr
                              CALL util_exit(FAILURE)
                         END IF
                         IF (id == ifol) THEN
                              istep=1
                              ic=ic+1
                              IF (MOD(ic,nskp) == 0) THEN
                                   WRITE(7,1001)t,ifol,xh,yh,zh,vxh,vyh,vzh
                              END IF
                              tmax=t
                         END IF
                    END DO
               CASE (FILT)
                    ! DEK - to be implemented
          END SELECT
          IF (istep == 0) EXIT
     END DO

     WRITE(*,*)"Tmax = ",tmax

     STOP

END PROGRAM tool_follow
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
