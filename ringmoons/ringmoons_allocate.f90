!**********************************************************************************************************************************
!
!  Unit Name   : ringmoons_allocate
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
!  Invocation  : CALL ringmoons_allocate(dt,ring,ring)
!
!  Notes       : Adapted from Andy Hesselbrock's ringmoons Python scripts
!
!**********************************************************************************************************************************
!  Author(s)   : David A. Minton  
!**********************************************************************************************************************************
SUBROUTINE ringmoons_allocate(ring)

! Modules
      USE module_parameters
      USE module_ringmoons
      USE module_ringmoons_interfaces, EXCEPT_THIS_ONE => ringmoons_allocate
      IMPLICIT NONE

! Arguments
      TYPE(ringmoons_ring),INTENT(INOUT) :: ring

! Internals

! Executable code
      allocate(ring%r(ring%N))
      allocate(ring%X(ring%N))
      allocate(ring%R_P(ring%N))
      allocate(ring%deltaA(ring%N))
      allocate(ring%m(ring%N))
      allocate(ring%sigma(ring%N))
      allocate(ring%nu(ring%N))
      allocate(ring%sigma_threshold(ring%N))
      allocate(ring%RR(ring%N))
      allocate(ring%I(ring%N))
      allocate(ring%w(ring%N))
      allocate(ring%Torque_to_disk(ring%N))
      

      RETURN

END SUBROUTINE ringmoons_allocate
