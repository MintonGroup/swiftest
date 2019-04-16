!**********************************************************************************************************************************
!
!  Unit Name   : io_write_line
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : io
!  Language    : Fortran 90/95
!
!  Description : Write a line (record) to output binary file
!
!  Input
!    Arguments : iu       : unit number associated with output binary file
!                id       : planet or test particle identifier
!                d1       : first quantity to write  (semimajor axis (pericentric distance for a parabola) or heliocentric x )
!                d2       : second quantity to write (eccentricity                                         or heliocentric y )
!                d3       : third quantity to write  (inclination                                          or heliocentric z )
!                d4       : fourth quantity to write (longitude of the ascending node                      or heliocentric vx)
!                d5       : fifth quantity to write  (argument of pericenter                               or heliocentric vy)
!                d6       : sixth quantity to write  (mean anomaly                                         or heliocentric vz)
!                out_type : format of output binary file
!                MASS     : optional mass (omitted for massless test particle)
!                RADIUS   : optional radius (omitted for massless test particle; must be present if optional MASS is present)
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : none
!    Terminal  : error message
!    File      : id                   : planet or test particle identifier
!                smass (or dmass)     : optional mass (omitted for massless test particle)
!                sradius (or dradius) : optional radius (omitted for massless test particle)
!                svec (or dvec)       : 6-vector of orbital elements or position and velocity components
!
!  Invocation  : CALL io_write_line(iu, id, d1, d2, d3, d4, d5, d6, out_type, MASS, RADIUS)
!
!  Notes       : Adapted from Hal Levison's Swift routine io_write_line.F
!
!**********************************************************************************************************************************
SUBROUTINE io_write_line(iu, id, d1, d2, d3, d4, d5, d6, out_type, MASS, RADIUS)

! Modules
     USE module_parameters
     USE module_fxdr
     USE module_interfaces, EXCEPT_THIS_ONE => io_write_line
     IMPLICIT NONE

! Arguments
     INTEGER(I4B), INTENT(IN)       :: iu, id
     REAL(DP), INTENT(IN)           :: d1, d2, d3, d4, d5, d6
     REAL(DP), OPTIONAL, INTENT(IN) :: MASS, RADIUS
     CHARACTER(*), INTENT(IN)       :: out_type

! Internals
     INTEGER(I4B)           :: ierr
     REAL(DP), DIMENSION(6) :: dvec
     REAL(SP), DIMENSION(6) :: svec
     REAL(SP)               :: smass, sradius
     LOGICAL(LGT)           :: lmass, lradius

! Executable code
     dvec(1) = d1; dvec(2) = d2; dvec(3) = d3; dvec(4) = d4; dvec(5) = d5; dvec(6) = d6
     svec(1) = d1; svec(2) = d2; svec(3) = d3; svec(4) = d4; svec(5) = d5; svec(6) = d6
     lmass = PRESENT(MASS)
     IF (lmass) THEN
          lradius = PRESENT(RADIUS)
          IF (.NOT. lradius) THEN
               WRITE(*, *) "SWIFTER Error:"
               WRITE(*, *) "   Subroutine io_write_line called with optional MASS but without optional RADIUS"
               CALL util_exit(FAILURE)
          END IF
          smass = MASS
          sradius = RADIUS
     END IF
     SELECT CASE (out_type)
          CASE (REAL4_TYPE)
               IF (lmass) THEN
                    WRITE(iu, IOSTAT = ierr) id, smass, sradius, svec
               ELSE
                    WRITE(iu, IOSTAT = ierr) id, svec
               END IF
               IF (ierr < 0) THEN
                    WRITE(*, *) "SWIFTER Error:"
                    WRITE(*, *) "   Unable to write binary file record"
                    CALL util_exit(FAILURE)
               END IF
          CASE (REAL8_TYPE)
               IF (lmass) THEN
                    WRITE(iu, IOSTAT = ierr) id, MASS, RADIUS, dvec
               ELSE
                    WRITE(iu, IOSTAT = ierr) id, dvec
               END IF
               IF (ierr < 0) THEN
                    WRITE(*, *) "SWIFTER Error:"
                    WRITE(*, *) "   Unable to write binary file record"
                    CALL util_exit(FAILURE)
               END IF
          CASE (XDR4_TYPE)
               ierr = ixdrint(iu, id)
               IF (ierr < 0) THEN
                    WRITE(*, *) "SWIFTER Error:"
                    WRITE(*, *) "   Unable to write binary file record"
                    CALL util_exit(FAILURE)
               END IF
               IF (lmass) THEN
                    ierr = ixdrreal(iu, smass)
                    IF (ierr < 0) THEN
                         WRITE(*, *) "SWIFTER Error:"
                         WRITE(*, *) "   Unable to write binary file record"
                         CALL util_exit(FAILURE)
                    END IF
                    ierr = ixdrreal(iu, sradius)
                    IF (ierr < 0) THEN
                         WRITE(*, *) "SWIFTER Error:"
                         WRITE(*, *) "   Unable to write binary file record"
                         CALL util_exit(FAILURE)
                    END IF
               END IF
               ierr = ixdrrmat(iu, 6, svec)
               IF (ierr < 0) THEN
                    WRITE(*, *) "SWIFTER Error:"
                    WRITE(*, *) "   Unable to write binary file record"
                    CALL util_exit(FAILURE)
               END IF
          CASE (XDR8_TYPE)
               ierr = ixdrint(iu, id)
               IF (ierr < 0) THEN
                    WRITE(*, *) "SWIFTER Error:"
                    WRITE(*, *) "   Unable to write binary file record"
                    CALL util_exit(FAILURE)
               END IF
               IF (lmass) THEN
                    ierr = ixdrdouble(iu, MASS)
                    IF (ierr < 0) THEN
                         WRITE(*, *) "SWIFTER Error:"
                         WRITE(*, *) "   Unable to write binary file record"
                         CALL util_exit(FAILURE)
                    END IF
                    ierr = ixdrdouble(iu, RADIUS)
                    IF (ierr < 0) THEN
                         WRITE(*, *) "SWIFTER Error:"
                         WRITE(*, *) "   Unable to write binary file record"
                         CALL util_exit(FAILURE)
                    END IF
               END IF
               ierr = ixdrdmat(iu, 6, dvec)
               IF (ierr < 0) THEN
                    WRITE(*, *) "SWIFTER Error:"
                    WRITE(*, *) "   Unable to write binary file record"
                    CALL util_exit(FAILURE)
               END IF
     END SELECT

     RETURN

END SUBROUTINE io_write_line
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
