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
subroutine ringmoons_allocate(ring,seeds)

! Modules
      use module_parameters
      use module_ringmoons
      use module_ringmoons_interfaces, EXCEPT_THIS_ONE => ringmoons_allocate
      implicit none

! Arguments
      type(ringmoons_ring),intent(inout) :: ring
      type(ringmoons_seeds),intent(inout) :: seeds

! Internals

! Executable code
      allocate(ring%r(ring%N))
      allocate(ring%rinner(ring%N))
      allocate(ring%router(ring%N))
      allocate(ring%X(ring%N))
      allocate(ring%X2(ring%N))
      allocate(ring%deltaA(ring%N))
      allocate(ring%Gm(ring%N))
      allocate(ring%Gsigma(ring%N))
      allocate(ring%nu(0:ring%N+1))
      allocate(ring%sigma_threshold(ring%N))
      allocate(ring%Iz(ring%N))
      allocate(ring%w(ring%N))
      allocate(ring%Torque_to_disk(ring%N))
      allocate(ring%r_hstar(ring%N))

      allocate(seeds%a(seeds%N))
      allocate(seeds%Gm(seeds%N))
      allocate(seeds%Rhill(seeds%N))
      allocate(seeds%rbin(seeds%N))
      

      return

end subroutine ringmoons_allocate
