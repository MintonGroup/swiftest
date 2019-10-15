!**********************************************************************************************************************************
!
!  Unit Name   : module_interfaces
!  Unit Type   : module
!  Project     : Swifter
!  Package     : module
!  Language    : Fortran 90/95
!
!  Description : Definition of interfaces of subroutines and functions used in Swifter package
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
!  Notes       : 
!
!**********************************************************************************************************************************
MODULE module_ringmoons_interfaces

      IMPLICIT NONE


      INTERFACE
         SUBROUTINE ringmoons_io_init_param()
            USE module_parameters
            IMPLICIT NONE
         END SUBROUTINE ringmoons_io_init_param
      END INTERFACE

      INTERFACE
         SUBROUTINE ringmoons_io_init_ring()
            USE module_parameters
            IMPLICIT NONE
         END SUBROUTINE
      END INTERFACE

      INTERFACE
         SUBROUTINE ringmoons_step(lfirst, t, npl, nplmax, symba_pl1P, j2rp2, j4rp4, eoffset, dt)
            USE module_parameters
            USE module_symba
            IMPLICIT NONE
            LOGICAL(LGT), INTENT(INOUT)                      :: lfirst
            INTEGER(I4B), INTENT(IN)                         :: npl, nplmax
            REAL(DP), INTENT(IN)                             :: t, j2rp2, j4rp4, dt
            REAL(DP), INTENT(INOUT)                          :: eoffset
            TYPE(symba_pl), POINTER                          :: symba_pl1P
         END SUBROUTINE ringmoons_step 
      END INTERFACE



END MODULE module_ringmoons_interfaces
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
