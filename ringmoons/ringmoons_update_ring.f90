!**********************************************************************************************************************************
!
!  Unit Name   : ringmoons_update_ring
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : io
!  Language    : Fortran 90/95
!
!  Description : Updates ring bin values for Hill ratio (r_hstar), viscosity (nu), Toomre parameter (Q), and optical depth (tau)
!                given values for particle mass (Gm_pdisk), radius (r_pdisk), and surface mass density (Gsigma)
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
!  Invocation  : CALL ringmoons_update_ring(dt,ring,ring)
!
!  Notes       : Adapted from Andy Hesselbrock's ringmoons Python scripts
!
!**********************************************************************************************************************************
!  Author(s)   : David A. Minton  
!**********************************************************************************************************************************
subroutine ringmoons_update_ring(swifter_pl1P,ring)

! Modules
   use module_parameters
   use module_swifter
   use module_ringmoons
   use module_ringmoons_interfaces, EXCEPT_THIS_ONE => ringmoons_update_ring
   implicit none

! Arguments
   type(swifter_pl),pointer                  :: swifter_pl1P
   type(ringmoons_ring),intent(inout) :: ring
   real(DP) :: rad_limit

! Internals

! Executable code
   rad_limit = 1.1_DP * RAD_LIMIT_CM / DU2CM
   ring%r_hstar(:) = ring%r(:) * (2 * ring%Gm_pdisk(:) /(3._DP * swifter_pl1P%mass))**(1._DP/3._DP) / (2 * ring%r_pdisk(:)) 
   where ((ring%r_pdisk(:) > rad_limit).and.(ring%Gm(:) > N_DISK_FACTOR * ring%Gm_pdisk))
      ring%Q(:) = ring%w(:) * ring%vrel_pdisk(:) / (3.36_DP * ring%Gsigma(:))
      ring%tau(:) = PI * ring%r_pdisk(:)**2 * ring%Gsigma(:) / ring%Gm_pdisk(:)
      ring%nu(:) = ringmoons_viscosity(ring%Gsigma(:), ring%Gm_pdisk(:), (ring%vrel_pdisk(:))**2, &
                                       ring%r_pdisk(:), ring%r_hstar(:), ring%Q(:), ring%tau(:), ring%w(:))
   elsewhere
      ring%Q(:) = huge(1._DP)
      ring%tau(:) = 0.0_DP
      ring%nu(:) = 0.0_DP
   end where

end subroutine ringmoons_update_ring
