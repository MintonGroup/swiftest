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
      integer(I4B)                        :: i,j,iFRL,rbin
      real(DP)                            :: rhill,mseed,dfac,fz_width

! Executable code
   
        
      dfac = exp(1._DP / seeds%N * log(ring%r_F / ring%FRL)) 
      seeds%Gminit = 3 * swifter_pl1P%mass * ((1.0_DP / FEEDING_ZONE_FACTOR) * (dfac - 1.0_DP))**3
      rbin = 1
      do i = 1,seeds%N
         seeds%a(i) = ring%FRL * (dfac)**(i - 1)
         do
            if ((seeds%a(i) >= ring%rinner(rbin)) .and. (seeds%a(i) < ring%router(rbin)) .or. (rbin == ring%N))  exit
            rbin = rbin + 1
         end do
         seeds%rbin(i) = rbin
         seeds%Gm(i) = seeds%Gminit
         seeds%Rhill(i) = seeds%a(i) * (seeds%Gm(i) / (3 * swifter_pl1P%mass))**(1.0_DP / 3.0_DP)
         ring%Gm(rbin) = max(ring%Gm(rbin) - seeds%Gminit, 0.0_DP)
         ring%Gsigma(rbin) = ring%Gm(rbin) / ring%deltaA(rbin)
         fz_width = FEEDING_ZONE_FACTOR * seeds%rhill(i)
         seeds%fz_bin_inner(i) = rbin
         seeds%fz_bin_outer(i) = rbin
         do j = seeds%fz_bin_inner(i),1,-1
            if (ring%rinner(j) < seeds%a(i) - fz_width) exit
            seeds%fz_bin_inner(i) = j
         end do
         do j = seeds%fz_bin_outer(i),ring%N
            if (ring%router(j) > seeds%a(i) + fz_width) exit
            seeds%fz_bin_outer(i) = j
         end do
         seeds%active(i) = .true.
      end do
      write(*,*) "Number of seeds created: ",seeds%N
      write(*,*) "Mass of seeds / mass of disk particles: ",seeds%Gminit / ring%Gm_pdisk

end subroutine ringmoons_seed_construct
