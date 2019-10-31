!**********************************************************************************************************************************
!
!  Unit Name   : ringmoons_timestep
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : io
!  Language    : Fortran 90/95
!
!  Description : Grows the seeds by accreting mass from within their local feeding zones from either the disk or other seeds
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
!  Invocation  : CALL ringmoons_timestep(dt,ring,ring)
!
!  Notes       : Adapted from Andy Hesselbrock's ringmoons Python scripts
!
!**********************************************************************************************************************************
!  Author(s)   : David A. Minton  
!**********************************************************************************************************************************
function ringmoons_timestep(swifter_pl1P,ring,seeds,dtin) result(dtout)

! Modules
      use module_parameters
      use module_swifter
      use module_ringmoons
      use module_ringmoons_interfaces, EXCEPT_THIS_ONE => ringmoons_timestep
      implicit none

! Arguments
      type(swifter_pl),pointer               :: swifter_pl1P
      type(ringmoons_ring), intent(in)       :: ring
      type(ringmoons_seeds), intent(in)      :: seeds
      real(DP), intent(in)                   :: dtin
      real(DP)                               :: dtout

! Internals
      integer(I4B)                           :: i
      real(DP),parameter                     :: SEED_GROWTH_FACTOR = 0.01_DP ! smallest increase in fractional mass allowable in a single time step
      real(DP)                               :: dGm_max,nu_max

! Executable code

      ! Start with viscous stability
      nu_max = maxval(abs(ring%nu))
      dtout = dtin
      if (nu_max > 0.0_DP) then
         dtout = min(dtout,ring%stability_factor / nu_max)  ! smallest timestep for the viscous evolution equation
      end if

      ! Now aim for seed growth accuracy
      dGm_max = maxval(ringmoons_seed_dMdt(ring,swifter_pl1P%mass,ring%Gsigma(seeds%rbin(:)), &
                       seeds%Gm(:),seeds%a(:)) / seeds%Gm(:),seeds%active)
      if (dGm_max > 0.0_DP) then
         dtout = min(dtout,SEED_GROWTH_FACTOR / dGm_max)  ! smallest timestep for the seed growth equation 
      end if

      return
end function ringmoons_timestep
