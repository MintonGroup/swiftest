!**********************************************************************************************************************************
!
!  Unit Name   : io_read_hdr
!  Unit Type   : function
!  Project     : Swifter
!  Package     : io
!  Language    : Fortran 90/95
!
!  Description : Read frame header from input binary file
!
!  Input
!    Arguments : iu        : unit number associated with input binary file
!                out_type  : format of input binary file
!    Terminal  : none
!    File      : t (or ttmp) : time
!                npl         : number of planets
!                ntp         : number of active test particles
!                iout_form   : specifier of data contained in input file (elements / heliocentric coordinates / filtered elements)
!
!  Output
!    Arguments : t         : time
!                npl       : number of planets
!                ntp       : number of active test particles
!                iout_form : specifier of data contained in input file (elements / heliocentric coordinates / filtered elements)
!    Terminal  : none
!    File      : none
!
!  Invocation  : istat = io_read_hdr(iu, t, npl, ntp, iout_form, out_type)
!
!  Notes       : Adapted from Hal Levison's Swift routine io_read_hdr.F
!
!                Function returns read error status (0 = OK, nonzero = ERROR)
!
!**********************************************************************************************************************************
FUNCTION io_read_hdr(iu, t, npl, ntp, iout_form, out_type)

! Modules
     USE module_parameters
     USE module_fxdr
     USE module_interfaces, EXCEPT_THIS_ONE => io_read_hdr
     IMPLICIT NONE

! Arguments
     INTEGER(I4B)              :: io_read_hdr
     INTEGER(I4B), INTENT(IN)  :: iu
     INTEGER(I4B), INTENT(OUT) :: npl, ntp, iout_form
     REAL(DP), INTENT(OUT)     :: t
     CHARACTER(*), INTENT(IN)  :: out_type

! Internals
     INTEGER(I4B)               :: ierr
     INTEGER(I4B), DIMENSION(3) :: nn
     REAL(SP)                   :: ttmp

! Executable code
     SELECT CASE (out_type)
          CASE (REAL4_TYPE)
               READ(iu, IOSTAT = ierr) ttmp, npl, ntp, iout_form
               io_read_hdr = ierr
               IF (ierr /= 0) RETURN
               t = ttmp
          CASE (REAL8_TYPE)
               READ(iu, IOSTAT = ierr) t, npl, ntp, iout_form
               io_read_hdr = ierr
          CASE (XDR4_TYPE)
               ierr = ixdrreal(iu, ttmp)
               io_read_hdr = ierr
               IF (ierr /= 0) RETURN
               t = ttmp
               ierr = ixdrimat(iu, 3, nn)
               io_read_hdr = ierr
               IF (ierr /= 0) RETURN
               npl = nn(1)
               ntp = nn(2)
               iout_form = nn(3)
          CASE (XDR8_TYPE)
               ierr = ixdrdouble(iu, t)
               io_read_hdr = ierr
               IF (ierr /= 0) RETURN
               ierr = ixdrimat(iu, 3, nn)
               io_read_hdr = ierr
               IF (ierr /= 0) RETURN
               npl = nn(1)
               ntp = nn(2)
               iout_form = nn(3)
     END SELECT

     RETURN

END FUNCTION io_read_hdr
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
