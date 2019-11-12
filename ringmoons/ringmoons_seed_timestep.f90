!**********************************************************************************************************************************
!
!  Unit Name   : ringmoons_seed_timestep
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : io
!  Language    : Fortran 90/95
!
!  Description : Calculates the maximum accurate seed_timestep for the seed growth and migration Runge-Kutta method
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
!  Invocation  : CALL ringmoons_seed_timestep(dt,ring,ring)
!
!  Notes       : Adapted from Andy Hesselbrock's ringmoons Python scripts
!
!**********************************************************************************************************************************
!  Author(s)   : David A. Minton  
!**********************************************************************************************************************************
function ringmoons_seed_timestep(swifter_pl1P,ring,seeds,dtin) result(dtout)

! Modules
      use module_parameters
      use module_swifter
      use module_ringmoons
      use module_ringmoons_interfaces, EXCEPT_THIS_ONE => ringmoons_seed_timestep
      implicit none

! Arguments
      type(swifter_pl),pointer               :: swifter_pl1P
      type(ringmoons_ring), intent(in)       :: ring
      type(ringmoons_seeds), intent(in)      :: seeds
      real(DP), intent(in)                   :: dtin
      real(DP)                               :: dtout

! Internals
      integer(I4B)                           :: i,nfz
      real(DP)                               :: dGm_max,da_max,sigavg,sig_max,nu_max
      real(DP)                               :: torque_term
      real(DP),dimension(seeds%N)            :: mdot
      

! Executable code

      dtout = dtin

      dGm_max = -1._DP
      call ringmoons_update_seeds(swifter_pl1P,ring,seeds)
      do i = 1,seeds%N
         nfz = seeds%fz_bin_outer(i) - seeds%fz_bin_inner(i) + 1
         sigavg = sum(ring%Gsigma(seeds%fz_bin_inner(i):seeds%fz_bin_outer(i))) / real(nfz, kind = DP)
         mdot(i) = ringmoons_seed_dMdt(ring,swifter_pl1P%mass,sigavg,seeds%Gm(i),seeds%a(i))
         dGm_max = max(dGm_max,mdot(i) / seeds%Gm(i))
      end do
      if (dGm_max > 0.0_DP) then
         dtout = min(dtout,RK_FACTOR / dGm_max)  ! smallest seed_timestep for the seed growth equation 
      end if

      ! Now aim for seed migration accuracy
         !write(*,*) 'migration'
      da_max = maxval(abs(ringmoons_seed_dadt(swifter_pl1P%mass,seeds%Gm(:),seeds%a(:),seeds%Torque(:),mdot(:))) / seeds%a(:),seeds%active) 
                     
      if (da_max > 0.0_DP) then
         dtout = min(dtout, RK_FACTOR / da_max) ! smallest step that keeps the body within a approximately single bin size 
      end if

      return
end function ringmoons_seed_timestep
