!**********************************************************************************************************************************
!
!  Unit Name   : ringmoons_ring_velocity_dispersion
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : io
!  Language    : Fortran 90/95
!
!  Description : Evolves the ring aggregate mass and velocity dispersion according to the predator/prey model of 
!                 Esposito et al. (2012)
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
!  Invocation  : CALL ringmoons_ring_velocity_dispersion(swifter_pl1P,ring)
!
!  Notes       : 
!
!**********************************************************************************************************************************
!  Author(s)   : David A. Minton  
!**********************************************************************************************************************************
subroutine ringmoons_ring_velocity_dispersion(swifter_pl1P,ring,lpredprey)

! Modules
   use module_parameters
   use module_swifter
   use module_ringmoons
   use module_ringmoons_interfaces, EXCEPT_THIS_ONE => ringmoons_ring_velocity_dispersion
   implicit none

! Arguments
   type(swifter_pl),pointer             :: swifter_pl1P
   type(ringmoons_ring), intent(inout)  :: ring
   logical(lgt), intent(in)             :: lpredprey

! Internals
   integer(I4B)                         :: i
   real(DP),dimension(0:ring%N+1)       :: kappa_rhstar,eta_rhstar


! Executable code


   if (lpredprey) then
      call ringmoons_ring_predprey(swifter_pl1P,ring) ! Evolve the size and velocity dispersion distribution of the ring 
                                                      ! following the predator/prey model of Esposito et al. (2012)
      ring%r_hstar(:) = ring%r(:) * (2 * ring%Gm_pdisk(:) /(3._DP * swifter_pl1P%mass))**(1._DP/3._DP) / (2 * ring%r_pdisk(:)) 
   else ! Constant particle mass and radius
      ring%r_hstar(:) = ring%r(:) * (2 * ring%Gm_pdisk(:) /(3._DP * swifter_pl1P%mass))**(1._DP/3._DP) / (2 * ring%r_pdisk(:)) 
      ! See Salmon et al. 2010 for this
      kappa_rhstar(:) = ringmoons_transition_function(ring%r_hstar(:))
      eta_rhstar(:) = 1._DP - kappa_rhstar(:)
      ring%vrel_pdisk(:) = kappa_rhstar(:) * sqrt(ring%Gm_pdisk(:) / ring%r_pdisk(:)) + eta_rhstar(:) * (2 * ring%r_pdisk(:) * ring%w(:))
   end if
   

return
end subroutine ringmoons_ring_velocity_dispersion
