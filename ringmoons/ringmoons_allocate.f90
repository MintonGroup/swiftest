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
      ring%r = 0.0_DP
      allocate(ring%X(0:ring%N + 1))
      ring%X = 0.0_DP
      allocate(ring%X2(0:ring%N + 1))
      ring%X2 = 0.0_DP
      allocate(ring%deltaA(0:ring%N + 1))
      ring%deltaA = 0.0_DP
      allocate(ring%Gm(0:ring%N + 1))
      ring%Gm = 0.0_DP
      allocate(ring%Gsigma(0:ring%N + 1))
      ring%Gsigma = 0.0_DP
      allocate(ring%nu(0:ring%N + 1))
      ring%nu = 0.0_DP
      allocate(ring%Iz(0:ring%N + 1))
      ring%Iz = 0.0_DP
      allocate(ring%w(0:ring%N + 1))
      ring%w = 0.0_DP
      allocate(ring%Torque(0:ring%N + 1))
      ring%Torque = 0.0_DP
      allocate(ring%r_hstar(0:ring%N + 1))
      ring%r_hstar = 0.0_DP

      allocate(seeds%active(seeds%N))
      seeds%active = .false.
      allocate(seeds%a(seeds%N))
      seeds%a = 0.0_DP
      allocate(seeds%Gm(seeds%N))
      seeds%Gm = 0.0_DP
      allocate(seeds%Rhill(seeds%N))
      seeds%Rhill = 0.0_DP
      allocate(seeds%rbin(seeds%N))
      seeds%rbin = 0
      allocate(seeds%fz_bin_inner(seeds%N))
      seeds%fz_bin_inner = 0
      allocate(seeds%fz_bin_outer(seeds%N))
      seeds%fz_bin_outer = 0
      allocate(seeds%Torque(seeds%N))
      seeds%Torque = 0.0_DP
      allocate(seeds%Ttide(seeds%N))
      seeds%Ttide = 0.0_DP
      

      return

end subroutine ringmoons_allocate
