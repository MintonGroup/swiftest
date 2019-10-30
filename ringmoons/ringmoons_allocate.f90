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
      allocate(ring%r(0:ring%N+1))
      allocate(ring%rinner(0:ring%N+1))
      allocate(ring%router(0:ring%N+1))
      allocate(ring%X(0:ring%N + 1))
      allocate(ring%X2(0:ring%N + 1))
      allocate(ring%deltaA(0:ring%N + 1))
      allocate(ring%Gm(0:ring%N + 1))
      allocate(ring%Gsigma(0:ring%N + 1))
      allocate(ring%nu(0:ring%N + 1))
      allocate(ring%sigma_threshold(0:ring%N + 1))
      allocate(ring%Iz(0:ring%N + 1))
      allocate(ring%w(0:ring%N + 1))
      allocate(ring%Torque(0:ring%N + 1))
      allocate(ring%r_hstar(0:ring%N + 1))

      allocate(seeds%active(seeds%N))
      allocate(seeds%a(seeds%N))
      allocate(seeds%Gm(seeds%N))
      allocate(seeds%Rhill(seeds%N))
      allocate(seeds%rbin(seeds%N))
      allocate(seeds%fz_bin_inner(seeds%N))
      allocate(seeds%fz_bin_outer(seeds%N))
      

      return

end subroutine ringmoons_allocate
