!**********************************************************************************************************************************
!
!  Unit Name   : ringmoons_seed_grow
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
!  Invocation  : CALL ringmoons_seed_grow(dt,ring,ring)
!
!  Notes       : Adapted from Andy Hesselbrock's ringmoons Python scripts
!
!**********************************************************************************************************************************
!  Author(s)   : David A. Minton  
!**********************************************************************************************************************************
subroutine ringmoons_seed_grow(swifter_pl1P,ring,seeds,dt)

! Modules
      use module_parameters
      use module_swifter
      use module_ringmoons
      use module_ringmoons_interfaces, EXCEPT_THIS_ONE => ringmoons_seed_grow
      implicit none

! Arguments
      type(swifter_pl),pointer               :: swifter_pl1P
      type(ringmoons_ring), intent(inout)    :: ring
      type(ringmoons_seeds), intent(inout)   :: seeds
      real(DP), intent(in)                   :: dt

! Internals
      integer(I4B)                           :: i,j,k,seed_bin,fz_bin_outer,fz_bin_inner,nfz
      real(DP)                               :: P,Q,R,S
      real(DP)                               :: Gmleft,fz_width,dGm,Gmdisk,Lfromring,Lseed_original,Ldiff
      real(DP),dimension(seeds%N)            :: dGmseeds

! Executable code

      ! Estimate the growth using Runge Kutta Calculate instantaneous growth ra
      do concurrent(i=1:seeds%N,seeds%active(i)) !local(P,Q,R,S)
         P = ringmoons_seed_dMdt(ring,swifter_pl1P%mass,ring%Gsigma(seeds%rbin(i)),seeds%Gm(i)                ,seeds%a(i))
         Q = ringmoons_seed_dMdt(ring,swifter_pl1P%mass,ring%Gsigma(seeds%rbin(i)),seeds%Gm(i) + 0.5_DP * P,seeds%a(i))
         R = ringmoons_seed_dMdt(ring,swifter_pl1P%mass,ring%Gsigma(seeds%rbin(i)),seeds%Gm(i) + 0.5_DP * Q,seeds%a(i))
         S = ringmoons_seed_dMdt(ring,swifter_pl1P%mass,ring%Gsigma(seeds%rbin(i)),seeds%Gm(i) +          R,seeds%a(i))
         dGmseeds(i) = (dt / 6._DP) * (P + Q + R + S)
      end do

      ! Get the mass out of the feeding zone (if it's there), starting from the bin the seed is in and working its way outward
      do i = 1, seeds%N
         if (.not.seeds%active(i)) cycle
         Lseed_original = seeds%Gm(i) * sqrt(swifter_pl1P%mass * seeds%a(i))

         seed_bin = seeds%rbin(i)
         fz_width = FEEDING_ZONE_FACTOR * seeds%rhill(i)
         do j = seeds%fz_bin_inner(i),1,-1
            if (ring%rinner(j) < seeds%a(i) - fz_width) exit
            seeds%fz_bin_inner(i) = j
         end do
         do j = seeds%fz_bin_outer(i),ring%N
            if (ring%router(j) > seeds%a(i) + fz_width) exit
            seeds%fz_bin_outer(i) = j
         end do
         nfz = seeds%fz_bin_outer(i) - seeds%fz_bin_inner(i) + 1 

         Gmleft = dGmseeds(i)
         Lfromring = 0.0_DP
         do j = seeds%fz_bin_inner(i),seeds%fz_bin_outer(i) ! loop over bins of the feeding zone and grab mass from them
            dGm = min(Gmleft / nfz,ring%Gm(j))
            ring%Gm(j) = ring%Gm(j) - dGm
            ring%Gsigma(j) = ring%Gm(j) / ring%deltaA(j)
            Gmleft = Gmleft - dGm
            ! Because the angular momentum of the bin is computed from the bin center, but the seed can be anywhere, 
            ! we have to make sure to keep track of it
            Lfromring = Lfromring + dGm * ring%Iz(j) * ring%w(j)
         end do
         ! If there is still mass left, take it from the innermost feeding zone
         dGm = min(Gmleft,ring%Gm(seed_bin)) 
         ring%Gm(seed_bin) = ring%Gm(seed_bin) - dGm
         Lfromring = Lfromring + dGm * ring%Iz(seed_bin) * ring%w(seed_bin)
         ring%Gsigma(seed_bin) = ring%Gm(seed_bin) / ring%deltaA(seed_bin)
         Gmleft = Gmleft - dGm
         seeds%Gm(i) = seeds%Gm(i) + (dGmseeds(i) - Gmleft)
         
        
         ! Conserve angular momentum 
         seeds%a(i) = ((Lseed_original + Lfromring) / seeds%Gm(i))**2 / swifter_pl1P%mass

         !I'm hungry! What's there to eat?!
         do j = 1,seeds%N
            if (seeds%active(j).and.(j /= i)) then
               if ((seeds%a(j) > seeds%a(i) - fz_width).and.(seeds%a(j) < seeds%a(i) + fz_width)) then ! This one is in the feeding zone
                  ! conserve orbital angular momentum
                  seeds%a(i) = ((seeds%Gm(i) * sqrt(seeds%a(i)) + seeds%Gm(j) * sqrt(seeds%a(j))) / (seeds%Gm(i) + seeds%Gm(j)))**2
                  ! conserve mass
                  seeds%Gm(i) = seeds%Gm(i) + seeds%Gm(j)
                  ! deactivate particle for now and position it at the FRL to potentially activate later
                  seeds%Gm(j) = 0.0_DP
                  seeds%a(j) = ring%FRL
                  seeds%active(j) = .false.
               end if
            end if
         end do
              
         !Update the Hill's sphere
         seeds%Rhill(i) = seeds%a(i) * (seeds%Gm(i) / (3 * swifter_pl1P%mass))**(1.0_DP / 3.0_DP) 

         !Update the seed location
         do while (seeds%a(i) < ring%rinner(seeds%rbin(i)))
            seeds%rbin(i) = seeds%rbin(i) - 1
         end do
         do while (seeds%a(i) >= ring%router(seeds%rbin(i)))
            seeds%rbin(i) = seeds%rbin(i) + 1
         end do   
         
      end do   
      
   
      return
end subroutine ringmoons_seed_grow
