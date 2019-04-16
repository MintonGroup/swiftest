!**********************************************************************************************************************************
!
!  Unit Name   : module_fxdr
!  Unit Type   : module
!  Project     : SWIFTER
!  Package     : module
!  Language    : Fortran 90/95
!
!  Description : Definition of interfaces of functions in the FXDR (Fortran eXternal Data Representation) library
!
!  Input
!    Arguments : N/A
!    Terminal  : N/A
!    File      : N/A
!
!  Output
!    Arguments : N/A
!    Terminal  : N/A
!    File      : N/A
!
!  Invocation  : N/A
!
!  Notes       : FXDR is a library, written and maintained by David W. Pierce, that enables calls to the XDR (eXternal Data
!                Representation) routines from Fortran.
!
!                Reference : http://meteora.ucsd.edu/~pierce/fxdr_home_page.html
!
!**********************************************************************************************************************************
MODULE module_fxdr

     IMPLICIT NONE

     INTERFACE
          FUNCTION initxdr(filename, mode, returnonerror)
               CHARACTER(*), INTENT(IN) :: filename
               CHARACTER(1), INTENT(IN) :: mode
               LOGICAL, INTENT(IN)      :: returnonerror
               INTEGER                  :: initxdr
          END FUNCTION initxdr
     END INTERFACE

     INTERFACE
          FUNCTION ixdrclose(ixdr)
               INTEGER, INTENT(IN) :: ixdr
               INTEGER             :: ixdrclose
          END FUNCTION ixdrclose
     END INTERFACE

     INTERFACE
          FUNCTION ixdrdmat(ixdrs, nels, dval)
               INTEGER, INTENT(IN)                           :: ixdrs, nels
               DOUBLE PRECISION, DIMENSION(nels), INTENT(IN) :: dval
               INTEGER                                       :: ixdrdmat
          END FUNCTION ixdrdmat
     END INTERFACE

     INTERFACE
          FUNCTION ixdrdouble(ixdrs, dval)
               INTEGER, INTENT(IN)          :: ixdrs
               DOUBLE PRECISION, INTENT(IN) :: dval
               INTEGER                      :: ixdrdouble
          END FUNCTION ixdrdouble
     END INTERFACE

     INTERFACE
          FUNCTION ixdrimat(ixdrs, nels, ival)
               INTEGER, INTENT(IN)                  :: ixdrs, nels
               INTEGER, DIMENSION(nels), INTENT(IN) :: ival
               INTEGER                              :: ixdrimat
          END FUNCTION ixdrimat
     END INTERFACE

     INTERFACE
          FUNCTION ixdrint(ixdrs, ival)
               INTEGER, INTENT(IN) :: ixdrs, ival
               INTEGER             :: ixdrint
          END FUNCTION ixdrint
     END INTERFACE

     INTERFACE
          FUNCTION ixdrreal(ixdrs, rval)
               INTEGER, INTENT(IN) :: ixdrs
               REAL, INTENT(IN)    :: rval
               INTEGER             :: ixdrreal
          END FUNCTION ixdrreal
     END INTERFACE

     INTERFACE
          FUNCTION ixdrrmat(ixdrs, nels, rval)
               INTEGER, INTENT(IN)               :: ixdrs, nels
               REAL, DIMENSION(nels), INTENT(IN) :: rval
               INTEGER                           :: ixdrrmat
          END FUNCTION ixdrrmat
     END INTERFACE

END MODULE module_fxdr
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
