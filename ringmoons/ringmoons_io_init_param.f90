!**********************************************************************************************************************************
!
!  Unit Name   : ringmoons_io_init_param
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : io
!  Language    : Fortran 90/95
!
!  Description : Read in parameters for the RINGMOONS system
!
!  Input
!    Arguments : 
!    Terminal  : none
!    File      : 
!
!  Output
!    Arguments : 
!    Terminal  : status, error messages
!    File      : 
!
!  Invocation  : CALL ringmoons_io_init_param()
!
!  Notes       : Adapted from Andy Hesselbrock's RINGMOONS Python scripts
!
!**********************************************************************************************************************************
SUBROUTINE ringmoons_io_init_param(rm)

! Modules
      USE module_parameters
      USE module_ringmoons
      USE module_ringmoons_interfaces, EXCEPT_THIS_ONE => ringmoons_io_init_param
      IMPLICIT NONE

! Arguments
      type(ringmoons_parameter),INTENT(INOUT) :: rm

! Internals

! Executable code

      RETURN

END SUBROUTINE ringmoons_io_init_param
!**********************************************************************************************************************************
!
!  Author(s)   : David A. Minton  
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
