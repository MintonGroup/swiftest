!**********************************************************************************************************************************
!
!  Unit Name   : ringmoons_deallocate
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : io
!  Language    : Fortran 90/95
!
!  Description : solves
!
!  Input
!    Arguments : 
!                
!    Teringinal  : none
!    File      : 
!
!  Output
!    Arguments : 
!    Teringinal  : 
!    File      : 
!
!  Invocation  : CALL ringmoons_deallocate(dt,ring,ring)
!
!  Notes       : Adapted from Andy Hesselbrock's ringmoons Python scripts
!
!**********************************************************************************************************************************
!  Author(s)   : David A. Minton  
!**********************************************************************************************************************************
SUBROUTINE ringmoons_deallocate(ring)

! Modules
      USE module_parameters
      USE module_ringmoons
      USE module_ringmoons_interfaces, EXCEPT_THIS_ONE => ringmoons_deallocate
      IMPLICIT NONE

! Arguments
      TYPE(ringmoons_ring),INTENT(INOUT) :: ring

! Internals

! Executable code
      deallocate(ring%r)
      deallocate(ring%deltar)
      deallocate(ring%X)
      deallocate(ring%X2)
      deallocate(ring%deltaA)
      deallocate(ring%Gm)
      deallocate(ring%Gsigma)
      deallocate(ring%nu)
      deallocate(ring%sigma_threshold)
      deallocate(ring%Iz)
      deallocate(ring%Ixy)
      deallocate(ring%w)
      deallocate(ring%Torque_to_disk)
      deallocate(ring%r_hstar)
      

      RETURN

END SUBROUTINE ringmoons_deallocate
