!**********************************************************************************************************************************
!
!  Unit Name   : ringmoons_ring_construct
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
!  Invocation  : CALL ringmoons_ring_construct(dt,ring,ring)
!
!  Notes       : Adapted from Andy Hesselbrock's ringmoons Python scripts
!
!**********************************************************************************************************************************
!  Author(s)   : David A. Minton  
!**********************************************************************************************************************************
subroutine ringmoons_ring_construct(swifter_pl1P,ring,seeds)

! Modules
      use module_parameters
      use module_ringmoons
      use module_swifter
      use module_ringmoons_interfaces, EXCEPT_THIS_ONE => ringmoons_ring_construct
      implicit none

! Arguments
      type(swifter_pl),pointer            :: swifter_pl1P
      type(ringmoons_ring), intent(inout) :: ring
      type(ringmoons_seeds), intent(inout) :: seeds

! Internals
      integer(I4B)                        :: i
      real(DP)                            :: Xlo
      real(DP)                            :: GMP, RP, rhoP


! Executable code
      GMP = swifter_pl1P%mass
      RP  = swifter_pl1P%radius
      rhoP = GMP / ((4.0_DP / 3.0_DP) * PI * RP**3)
      ring%nu = 0.0_DP
      ring%X_I = 2 * sqrt(ring%r_I)
      ring%X_F = 2 * sqrt(ring%r_F)
      ring%deltaX = (ring%X_F - ring%X_I) / ring%N
      ring%rho_pdisk(:) = ring%Gm_pdisk(:) / ((4.0_DP / 3.0_DP) * PI * ring%r_pdisk(:)**3)
      ring%FRL = 2.456_DP * RP * (rhoP / ring%rho_pdisk(ring%N))**(1._DP / 3._DP)
      ring%RRL = 1.44_DP  * RP * (rhoP / ring%rho_pdisk(ring%N))**(1._DP / 3._DP)
      ring%iFRL = ringmoons_ring_bin_finder(ring,ring%FRL)
      ring%iRRL = ringmoons_ring_bin_finder(ring,ring%RRL)

      do i = 0,ring%N + 1
         ! Set up X coordinate system (see Bath & Pringle 1981)
         Xlo = ring%X_I + ring%deltaX * (i - 1)
   
         ring%X(i) = Xlo + 0.5_DP * ring%deltaX
      end do
      ring%X2(:) = ring%X(:)**2

      ! Convert X to r
      ring%r(:) = 0.25_DP * (ring%X(:))**2
        
      ! Factors to convert surface mass density into mass 
      ring%deltaA(:) = 0.25_DP * PI * ring%X(:)**3 * ring%deltaX !2 * PI * deltar * ring%r(i)
      ring%Gm(:) = ring%Gsigma(:) * ring%deltaA(:)
      
      ! Specific moment of inertia of the ring bin
      ring%Iz(:) = (ring%r(:))**2
      ring%w(:) = sqrt(GMP / ring%r(:)**3)

      ring%Torque(:) = 0.0_DP

      call ringmoons_update_ring(swifter_pl1P,ring)

      seeds%Rhill(1:seeds%N)  = seeds%a(1:seeds%N) * (seeds%Gm(1:seeds%N) / (3 * GMP))**(1.0_DP / 3.0_DP)
      seeds%rbin(1:seeds%N)   = ringmoons_ring_bin_finder(ring,seeds%a(1:seeds%N))
      seeds%fz_bin_inner(1:seeds%N) = seeds%rbin(1:seeds%N) !ringmoons_ring_bin_finder(ring,seeds%a(i) - fz_width)
      seeds%fz_bin_outer(1:seeds%N) = seeds%rbin(1:seeds%N) !ringmoons_ring_bin_finder(ring,seeds%a(i) + fz_width)

      return

end subroutine ringmoons_ring_construct
