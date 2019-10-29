!**********************************************************************************************************************************
!
!  Unit Name   : ringmoons_seed_construct
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : io
!  Language    : Fortran 90/95
!
!  Description : Constructs the satellite seeds. Given the number of seeds, this will find the mass of each seed needed to generate
!                seeds spaced by spacing_factor times the Hill's radius between the Fluid Roche Limit (FRL) and the outer bin of the
!                ring
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
!  Invocation  : CALL ringmoons_seed_construct(dt,ring,ring)
!
!  Notes       : Adapted from Andy Hesselbrock's ringmoons Python scripts
!
!**********************************************************************************************************************************
!  Author(s)   : David A. Minton  
!**********************************************************************************************************************************
subroutine ringmoons_seed_construct(swifter_pl1P,ring,seeds)

! Modules
      use module_parameters
      use module_swifter
      use module_ringmoons
      use module_ringmoons_interfaces, EXCEPT_THIS_ONE => ringmoons_seed_construct
      implicit none

! Arguments
      type(swifter_pl),pointer :: swifter_pl1P
      type(ringmoons_ring), intent(inout) :: ring
      type(ringmoons_seeds), intent(inout) :: seeds

! Internals
      integer(I4B)                        :: i,j,seed_bin
      real(DP)                            :: dfac, dGm, Gmleft
      real(DP), dimension(seeds%N)        :: fz_width 
      integer(I4B),dimension(seeds%N)     :: nfz

! Executable code
   
       
      seeds%active = .true. 

      ! Find the initial seed spacing
      dfac = exp(1._DP / seeds%N * log(ring%r_F / ring%FRL)) 
      do i = 1,seeds%N
         seeds%a(i) = ring%FRL * (dfac)**(i)
      end do
      ! Set the initial mass of the seed so that their feeding zones touch
      seeds%Gminit = 3 * swifter_pl1P%mass * ((1.0_DP / FEEDING_ZONE_FACTOR) * (dfac - 1.0_DP))**3
      seeds%Gm = seeds%Gminit
      seeds%Rhill = seeds%a * (seeds%Gm / (3 * swifter_pl1P%mass))**(1.0_DP / 3.0_DP)

      ! Find the corresponding ring bins that the seeds are embedded (returns the final ring bin if it is outside, which always has 0 mass)
      seeds%rbin = ringmoons_ring_bin_finder(ring,seeds%a)
      fz_width = FEEDING_ZONE_FACTOR * seeds%Rhill
      seeds%fz_bin_inner = ringmoons_ring_bin_finder(ring,seeds%a - FEEDING_ZONE_FACTOR * seeds%rhill)
      seeds%fz_bin_outer = ringmoons_ring_bin_finder(ring,seeds%a + FEEDING_ZONE_FACTOR * seeds%rhill)
      nfz = seeds%fz_bin_outer - seeds%fz_bin_inner + 1

      ! Take mass out of local ring material if it is available
      do i = 1,seeds%N
         seed_bin = seeds%rbin(i)
         Gmleft = seeds%Gm(i)
         do j = seeds%fz_bin_inner(i),seeds%fz_bin_outer(i) ! loop over bins of the feeding zone and grab mass from them
            dGm = min(Gmleft / nfz(i),ring%Gm(j))
            ring%Gm(j) = ring%Gm(j) - dGm
            ring%Gsigma(j) = ring%Gm(j) / ring%deltaA(j)
            Gmleft = Gmleft - dGm
         end do
         ! If there is still mass left, take it from the innermost feeding zone
         dGm = min(Gmleft,ring%Gm(seed_bin)) 
         ring%Gm(seed_bin) = ring%Gm(seed_bin) - dGm
         ring%Gsigma(seed_bin) = ring%Gm(seed_bin) / ring%deltaA(seed_bin)
      end do
      write(*,*) "Number of seeds created: ",seeds%N
      write(*,*) "Mass of seeds / mass of disk particles: ",seeds%Gminit / ring%Gm_pdisk

end subroutine ringmoons_seed_construct
