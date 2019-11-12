!**********************************************************************************************************************************
!
!  Unit Name   : ringmoons_seed_spawn
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : io
!  Language    : Fortran 90/95
!
!  Description : Spawns a new satellite at a given location
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
!  Invocation  : CALL ringmoons_seed_spawn(dt,ring,ring)
!
!  Notes       : Adapted from Andy Hesselbrock's ringmoons Python scripts
!
!**********************************************************************************************************************************
!  Author(s)   : David A. Minton  
!**********************************************************************************************************************************
subroutine ringmoons_seed_spawn(swifter_pl1P,ring,seeds,a,Gm)

! Modules
      use module_parameters
      use module_swifter
      use module_ringmoons
      use module_ringmoons_interfaces, EXCEPT_THIS_ONE => ringmoons_seed_spawn
      implicit none

! Arguments
      type(swifter_pl),pointer               :: swifter_pl1P
      type(ringmoons_ring), intent(inout)    :: ring
      type(ringmoons_seeds), intent(inout)   :: seeds
      real(DP), intent(in)                   :: a, Gm

! Internals
      integer(I4B)                        :: i,j,seed_bin,nfz
      real(DP)                            :: fz_width
      type(ringmoons_seeds)               :: new_seeds
      type(ringmoons_ring)                :: tmpring

! Executable code
 
      seed_bin = seeds%N + 1  
      !if (seed_bin > seeds%N) then ! If no previously generated inactive seeds, we'll take advantage of Fortran 2003 automatic allocation and tack it on to the end  
         tmpring%N = 0
         new_seeds%N = seed_bin
         call ringmoons_allocate(tmpring,new_seeds)
         new_seeds%active(1:seeds%N) = seeds%active(1:seeds%N)
         new_seeds%a(1:seeds%N) = seeds%a(1:seeds%N)
         new_seeds%Gm(1:seeds%N) = seeds%Gm(1:seeds%N)
         new_seeds%Rhill(1:seeds%N) = seeds%Rhill(1:seeds%N)
         new_seeds%rbin(1:seeds%N) = seeds%rbin(1:seeds%N)
         new_seeds%fz_bin_inner(1:seeds%N) = seeds%fz_bin_inner(1:seeds%N)
         new_seeds%fz_bin_outer(1:seeds%N) = seeds%fz_bin_outer(1:seeds%N)
         new_seeds%Torque(1:seeds%N) = seeds%Torque(1:seeds%N)
         new_seeds%TTide(1:seeds%N) = seeds%Ttide(1:seeds%N)
         seeds%active = new_seeds%active
         seeds%a = new_seeds%a
         seeds%Gm = new_seeds%Gm
         seeds%rbin = new_seeds%rbin
         seeds%Rhill = new_seeds%Rhill
         seeds%fz_bin_inner = new_seeds%fz_bin_inner
         seeds%fz_bin_outer = new_seeds%fz_bin_outer
         seeds%Torque = new_seeds%Torque
         seeds%Ttide = new_seeds%Ttide
         seeds%N = new_seeds%N
         call ringmoons_deallocate(tmpring,new_seeds)
      !end if

      i = seed_bin 
      seeds%active(i) = .true. 
      j = ringmoons_ring_bin_finder(ring,a)
      seeds%rbin(i) = j 
      seeds%Gm(i) = min(Gm,ring%Gm(j))
      ! Adjust the semimajor axis in order to conserve angular momentum 
      seeds%a(i) = (ring%Iz(j) * ring%w(j))**2 / (swifter_pl1P%mass + seeds%Gm(i))
      seeds%Rhill(i) = seeds%a(i) * (seeds%Gm(i) / (3 * swifter_pl1P%mass))**(1.0_DP / 3.0_DP)

      ! Take away the mass from the ring
      ring%Gm(j) = ring%Gm(j) - seeds%Gm(i)
      ring%Gsigma(j) = ring%Gm(j) / ring%deltaA(j)
      fz_width = FEEDING_ZONE_FACTOR * seeds%Rhill(i)
      seeds%fz_bin_inner(i) = ringmoons_ring_bin_finder(ring,seeds%a(i) - fz_width)
      seeds%fz_bin_outer(i) = ringmoons_ring_bin_finder(ring,seeds%a(i) + fz_width)

      return

end subroutine ringmoons_seed_spawn
