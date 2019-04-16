!**********************************************************************************************************************************
!
!  Unit Name   : io_write_hdr
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : io
!  Language    : Fortran 90/95
!
!  Description : Write frame header to output binary file
!
!  Input
!    Arguments : iu        : unit number associated with output binary file
!                t         : time
!                npl       : number of planets
!                ntp       : number of active test particles
!                iout_form : specifier of data to write to output file (elements / heliocentric coordinates / filtered elements)
!                out_type  : format of output binary file
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : none
!    Terminal  : error message
!    File      : t (or ttmp) : time
!                npl         : number of planets
!                ntp         : number of active test particles
!                iout_form   : specifier of data to write to output file (elements / heliocentric coordinates / filtered elements)
!
!  Invocation  : CALL io_write_hdr(iu, t, npl, ntp, iout_form, out_type)
!
!  Notes       : Adapted from Hal Levison's Swift routine io_write_hdr.F
!
!**********************************************************************************************************************************
SUBROUTINE io_write_hdr(iu, t, npl, ntp, iout_form, out_type)

! Modules
     USE module_parameters
     USE module_fxdr
     USE module_interfaces, EXCEPT_THIS_ONE => io_write_hdr
     IMPLICIT NONE

! Arguments
     INTEGER(I4B), INTENT(IN) :: iu, npl, ntp, iout_form
     REAL(DP), INTENT(IN)     :: t
     CHARACTER(*), INTENT(IN) :: out_type

! Internals
     INTEGER(I4B)               :: ierr
     INTEGER(I4B), DIMENSION(3) :: nn
     REAL(SP)                   :: ttmp

! Executable code
     ttmp = t
     nn(1) = npl
     nn(2) = ntp
     nn(3) = iout_form
     SELECT CASE (out_type)
          CASE (REAL4_TYPE)
               WRITE(iu, IOSTAT = ierr) ttmp, npl, ntp, iout_form
               IF (ierr < 0) THEN
                    WRITE(*, *) "SWIFTER Error:"
                    WRITE(*, *) "   Unable to write binary file header"
                    CALL util_exit(FAILURE)
               END IF
          CASE (REAL8_TYPE)
               WRITE(iu, IOSTAT = ierr) t, npl, ntp, iout_form
               IF (ierr < 0) THEN
                    WRITE(*, *) "SWIFTER Error:"
                    WRITE(*, *) "   Unable to write binary file header"
                    CALL util_exit(FAILURE)
               END IF
          CASE (XDR4_TYPE)
               ierr = ixdrreal(iu, ttmp)
               IF (ierr < 0) THEN
                    WRITE(*, *) "SWIFTER Error:"
                    WRITE(*, *) "   Unable to write binary file header"
                    CALL util_exit(FAILURE)
               END IF
               ierr = ixdrimat(iu, 3, nn)
               IF (ierr < 0) THEN
                    WRITE(*, *) "SWIFTER Error:"
                    WRITE(*, *) "   Unable to write binary file header"
                    CALL util_exit(FAILURE)
               END IF
          CASE (XDR8_TYPE)
               ierr = ixdrdouble(iu, t)
               IF (ierr < 0) THEN
                    WRITE(*, *) "SWIFTER Error:"
                    WRITE(*, *) "   Unable to write binary file header"
                    CALL util_exit(FAILURE)
               END IF
               ierr = ixdrimat(iu, 3, nn)
               IF (ierr < 0) THEN
                    WRITE(*, *) "SWIFTER Error:"
                    WRITE(*, *) "   Unable to write binary file header"
                    CALL util_exit(FAILURE)
               END IF
     END SELECT

     RETURN

END SUBROUTINE io_write_hdr
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
