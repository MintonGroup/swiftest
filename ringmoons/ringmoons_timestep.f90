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
      integer(I4B)                           :: i,nfz
      real(DP),parameter                     :: RK_FACTOR = 0.01_DP ! smallest increase in fractional mass allowable in a single time step
      real(DP)                               :: dGm_max,nu_max,da_max,dadot_max,sigavg,Torque_max
      

! Executable code

      dtout = dtin

      ! Start with viscous stability
      nu_max = maxval(abs(ring%nu))
      if (nu_max > 0.0_DP) then
         dtout = min(dtout,0.5_DP * ring%stability_factor / nu_max)  ! smallest timestep for the viscous evolution equation
         !write(*,*) 'Viscous dt/dtin: ', ring%stability_factor / nu_max / dtin
      end if

      ! Satellite Torque stability
      !Torque_max = maxval(abs(ring%Torque))
      !if (Torque_max > 0.0_DP) then
      !   !write(*,*) dtout
      !   dtout = min(dtout, 1e-3 * ring%deltaX * sqrt(swifter_pl1P%mass) * minval(ring%deltaA(:)) / Torque_max)
      !   !write(*,*) dtout
      !end if 
       

      ! Now aim for seed growth accuracy
      dGm_max = -1._DP
      do concurrent (i = 1:seeds%N, seeds%active(i))
         nfz = seeds%fz_bin_outer(i) - seeds%fz_bin_inner(i) + 1
         sigavg = sum(ring%Gsigma(seeds%fz_bin_inner(i):seeds%fz_bin_outer(i))) / real(nfz, kind = DP)
         dGm_max = max(dGm_max,ringmoons_seed_dMdt(ring,swifter_pl1P%mass,sigavg,seeds%Gm(i),seeds%a(i)) / seeds%Gm(i))
      end do
      if (dGm_max > 0.0_DP) then
         dtout = min(dtout,RK_FACTOR / dGm_max)  ! smallest timestep for the seed growth equation 
      end if

      ! Now aim for seed migration accuracy
      dadot_max = maxval(abs(ringmoons_seed_dadt(swifter_pl1P%mass,seeds%Gm(:),seeds%a(:),seeds%Torque(:))),seeds%active) 
                     
      if (da_max > 0.0_DP) then
         dtout = min(dtout, (ring%deltaX / 2._DP)**2 / dadot_max) ! smallest step that keeps the body within a approximately single bin size 
      end if

      return
end function ringmoons_timestep
