!**********************************************************************************************************************************
!
!  Unit Name   : io_read_line
!  Unit Type   : function
!  Project     : Swifter
!  Package     : io
!  Language    : Fortran 90/95
!
!  Description : Read a line (record) from input binary file
!
!  Input
!    Arguments : iu       : unit number associated with input binary file
!                out_type : format of input binary file
!    Terminal  : none
!    File      : id                  : planet or test particle identifier
!                smass (or MASS)     : optional mass (omitted for massless test particle)
!                sradius (or RADIUS) : optional radius (omitted for massless test particle)
!                svec (or dvec)      : 6-vector of orbital elements or position and velocity components
!
!  Output
!    Arguments : id       : planet or test particle identifier
!                d1       : first quantity  (semimajor axis (pericentric distance for a parabola) or heliocentric x )
!                d2       : second quantity (eccentricity                                         or heliocentric y )
!                d3       : third quantity  (inclination                                          or heliocentric z )
!                d4       : fourth quantity (longitude of the ascending node                      or heliocentric vx)
!                d5       : fifth quantity  (argument of pericenter                               or heliocentric vy)
!                d6       : sixth quantity  (mean anomaly                                         or heliocentric vz)
!                MASS     : optional mass (omitted for massless test particle)
!                RADIUS   : optional radius (omitted for massless test particle)
!    Terminal  : none
!    File      : none
!
!  Invocation  : istat = io_read_line(iu, id, d1, d2, d3, d4, d5, d6, out_type, MASS, RADIUS)
!
!  Notes       : Adapted from Hal Levison's Swift function io_read_line
!
!                Function returns read error status (0 = OK, nonzero = ERROR)
!
!**********************************************************************************************************************************
FUNCTION io_read_line(iu, id, d1, d2, d3, d4, d5, d6, out_type, MASS, RADIUS)

! Modules
     USE module_parameters
     USE module_fxdr
     USE module_interfaces, EXCEPT_THIS_ONE => io_read_line
     IMPLICIT NONE

! Arguments
     INTEGER(I4B)                    :: io_read_line
     INTEGER(I4B), INTENT(IN)        :: iu
     INTEGER(I4B), INTENT(OUT)       :: id
     REAL(DP), INTENT(OUT)           :: d1, d2, d3, d4, d5, d6
     REAL(DP), OPTIONAL, INTENT(OUT) :: MASS, RADIUS
     CHARACTER(*), INTENT(IN)        :: out_type

! Internals
     LOGICAL(LGT)           :: lmass, lradius
     INTEGER(I4B)           :: ierr
     REAL(SP)               :: smass, sradius
     REAL(SP), DIMENSION(6) :: svec
     REAL(DP), DIMENSION(6) :: dvec

! Executable code
     lmass = PRESENT(MASS)
     IF (lmass) THEN
          lradius = PRESENT(RADIUS)
          IF (.NOT. lradius) THEN
               WRITE(*, *) "SWIFTER Error:"
               WRITE(*, *) "   Subroutine io_read_line called with optional MASS but without optional RADIUS"
               CALL util_exit(FAILURE)
          END IF
     END IF
     SELECT CASE (out_type)
          CASE (REAL4_TYPE)
               IF (lmass) THEN
                    READ(iu, IOSTAT = ierr) id, smass, sradius, svec
               ELSE
                    READ(iu, IOSTAT = ierr) id, svec
               END IF
               io_read_line = ierr
               IF (ierr /= 0) RETURN
               IF (lmass) MASS = smass
               d1 = svec(1); d2 = svec(2); d3 = svec(3); d4 = svec(4); d5 = svec(5); d6 = svec(6)
          CASE (REAL8_TYPE)
               IF (lmass) THEN
                    READ(iu, IOSTAT = ierr) id, MASS, RADIUS, dvec
               ELSE
                    READ(iu, IOSTAT = ierr) id, dvec
               END IF
               io_read_line = ierr
               IF (ierr /= 0) RETURN
               d1 = dvec(1); d2 = dvec(2); d3 = dvec(3); d4 = dvec(4); d5 = dvec(5); d6 = dvec(6)
          CASE (XDR4_TYPE)
               ierr = ixdrint(iu, id)
               io_read_line = ierr
               IF (ierr /= 0) RETURN
               IF (lmass) THEN
                    ierr = ixdrreal(iu, smass)
                    io_read_line = ierr
                    IF (ierr /= 0) RETURN
                    MASS = smass
                    ierr = ixdrreal(iu, sradius)
                    io_read_line = ierr
                    IF (ierr /= 0) RETURN
                    RADIUS = sradius
               END IF
               ierr = ixdrrmat(iu, 6, svec)
               io_read_line = ierr
               IF (ierr /= 0) RETURN
               d1 = svec(1); d2 = svec(2); d3 = svec(3); d4 = svec(4); d5 = svec(5); d6 = svec(6)
          CASE (XDR8_TYPE)
               ierr = ixdrint(iu, id)
               io_read_line = ierr
               IF (ierr /= 0) RETURN
               IF (lmass) THEN
                    ierr = ixdrdouble(iu, MASS)
                    io_read_line = ierr
                    IF (ierr /= 0) RETURN
                    ierr = ixdrdouble(iu, RADIUS)
                    io_read_line = ierr
                    IF (ierr /= 0) RETURN
               END IF
               ierr = ixdrdmat(iu, 6, dvec)
               io_read_line = ierr
               IF (ierr /= 0) RETURN
               d1 = dvec(1); d2 = dvec(2); d3 = dvec(3); d4 = dvec(4); d5 = dvec(5); d6 = dvec(6)
     END SELECT

     RETURN

END FUNCTION io_read_line
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
