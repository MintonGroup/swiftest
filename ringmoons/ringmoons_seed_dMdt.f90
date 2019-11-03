!**********************************************************************************************************************************
!
!  Unit Name   : ringmoons_seed_dMdt
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : io
!  Language    : Fortran 90/95
!
!  Description : Models the growth rate of a seed embedded within the ring using a runaway growth model with a factor that corrects
!                for distance from the FRL
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
!  Invocation  : CALL ringmoons_seed_dMdt(dt,ring,ring)
!
!  Notes       : Adapted from Andy Hesselbrock's ringmoons Python scripts, but uses a different growth model
!
!**********************************************************************************************************************************
!  Author(s)   : David A. Minton  
!**********************************************************************************************************************************
elemental function ringmoons_seed_dMdt(ring,GMP,Gsigma,Gmseed,a) result(Gmdot)
!function ringmoons_seed_dMdt(ring,GMP,Gsigma,Gmseed,a) result(Gmdot)

! Modules
      use module_parameters
      use module_ringmoons
      use module_ringmoons_interfaces, EXCEPT_THIS_ONE => ringmoons_seed_dMdt
      implicit none

! Arguments
      type(ringmoons_ring), intent(in)       :: ring
      real(DP), intent(in)                   :: GMP,Gsigma,Gmseed,a
      real(DP)                               :: Gmdot

! Internals
      integer(I4B)                           :: i
      real(DP)                               :: Gmeff,C
      real(DP),parameter                     :: eff2   = 1e-7_DP ! This term gets the growth rate to match up closely to Andy's
      real(DP),parameter                     :: growth_exponent = 4._DP / 3._DP
      real(DP), parameter                    :: SIGLIMIT = 1e-100_DP

! Executable code
     
      if (a < ring%FRL) then
         Gmdot = 0.0_DP
      else
         C = 12 * PI**(2._DP / 3._DP) * (3._DP / (4 * ring%rho_pdisk))**(1._DP / 3._DP) / sqrt(GMP) 
         if (Gsigma > SIGLIMIT) then
            Gmeff = GMseed * (1._DP - (ring%FRL / a)**3) ! Reduces the growth rate near the FRL
            Gmdot = C * Gsigma / (eff2 * sqrt(a)) * Gmeff**(growth_exponent)
         else
            Gmdot = 0.0_DP
         end if
      end if
    
      return
end function ringmoons_seed_dMdt
