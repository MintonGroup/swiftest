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
subroutine ringmoons_ring_construct(swifter_pl1P,ring)

! Modules
      use module_parameters
      use module_ringmoons
      use module_swifter
      use module_ringmoons_interfaces, EXCEPT_THIS_ONE => ringmoons_ring_construct
      implicit none

! Arguments
      type(swifter_pl),pointer            :: swifter_pl1P
      type(ringmoons_ring), intent(inout) :: ring

! Internals
      integer(I4B)                        :: i
      real(DP)                            :: Xlo,Xhi,rlo,rhi,rhill,deltar
      real(DP)                            :: GMP, RP, rhoP

! Executable code
      GMP = swifter_pl1P%mass
      RP  = swifter_pl1P%radius
      rhoP = GMP / ((4.0_DP / 3.0_DP) * PI * RP**3)
      ring%nu = 0.0_DP
      ring%deltaX = (2 * sqrt(ring%r_F) - 2 * sqrt(ring%r_I)) / ring%N
      ring%rho_pdisk = ring%Gm_pdisk / ((4.0_DP / 3.0_DP) * PI * ring%r_pdisk**3)
      ring%FRL = 2.456_DP * RP * (rhoP / ring%rho_pdisk)**(1._DP / 3._DP)
      ring%RRL = 1.44_DP  * RP * (rhoP / ring%rho_pdisk)**(1._DP / 3._DP)
      ring%iFRL = ringmoons_ring_bin_finder(ring,ring%FRL)
      ring%iRRL = ringmoons_ring_bin_finder(ring,ring%RRL)
      do i = 0,ring%N + 1
         ! Set up X coordinate system (see Bath & Pringle 1981)
         Xlo = 2 * sqrt(ring%r_I) + ring%deltaX * (1._DP * i - 1._DP)
         Xhi = Xlo + ring%deltaX
   
         ring%X(i) = Xlo + 0.5_DP * ring%deltaX
         ring%X2(i) = ring%X(i)**2

         ! Convert X to r
         rlo = (0.5_DP * Xlo)**2
         rhi = (0.5_DP * Xhi)**2
         deltar = rhi - rlo
         ring%r(i) = (0.5_DP * ring%X(i))**2
         ring%rinner(i) = rlo
         ring%router(i) = rhi
        
         ! Factors to convert surface mass density into mass 
         ring%deltaA(i) = 2 * PI * deltar * ring%r(i)
         ring%Gm(i) = ring%Gsigma(i) * ring%deltaA(i)
      
         ! Specific moment of inertia of the ring bin
         ring%Iz(i) = 0.5_DP * (rlo**2 + rhi**2)
         ring%w(i) = sqrt(GMP / ring%r(i)**3)

         ring%Torque(i) = 0.0_DP
         rhill = ring%r(i) * (2 * ring%Gm_pdisk /(3._DP * GMP))**(1._DP/3._DP) ! See Salmon et al. 2010 for this
         ring%r_hstar(i) = rhill / (2 * ring%r_pdisk)  
      end do

      
      
      

      return

end subroutine ringmoons_ring_construct
