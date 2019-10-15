!**********************************************************************************************************************************
!
!  Unit Name   : ringmoons_io_init_ring
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : io
!  Language    : Fortran 90/95
!
!  Description : Read in ring initial condition data
!
!  Input
!    Arguments : 
!                
!    Terminal  : none
!    File      : 
!
!  Output
!    Arguments : 
!    Terminal  : 
!    File      : 
!
!  Invocation  : CALL ringmoons_io_init_ring()
!
!  Notes       : Adapted from Andy Hesselbrock's ringmoons Python scripts
!
!**********************************************************************************************************************************
SUBROUTINE ringmoons_io_init_ring(rm,ring)

! Modules
      USE module_parameters
      USE module_ringmoons
      USE module_ringmoons_interfaces, EXCEPT_THIS_ONE => ringmoons_io_init_ring
      IMPLICIT NONE

! Arguments
      TYPE(ringmoons_parameter),INTENT(IN) :: rm
      TYPE(ringmoons_ring_bin),DIMENSION(:),INTENT(INOUT) :: ring

! Internals

! Executable code


      RETURN

END SUBROUTINE ringmoons_io_init_ring
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
