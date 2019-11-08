!**********************************************************************************************************************************
!
!  Unit Name   : ringmoons_seed_evolve
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
!  Invocation  : CALL ringmoons_seed_evolve(dt,ring,ring)
!
!  Notes       : Adapted from Andy Hesselbrock's ringmoons Python scripts
!
!**********************************************************************************************************************************
!  Author(s)   : David A. Minton  
!**********************************************************************************************************************************
subroutine ringmoons_seed_evolve(swifter_pl1P,ring,seeds,dt,stepfail)

! Modules
   use module_parameters
   use module_swifter
   use module_ringmoons
   use module_ringmoons_interfaces, EXCEPT_THIS_ONE => ringmoons_seed_evolve
   implicit none

! Arguments
   type(swifter_pl),pointer               :: swifter_pl1P
   type(ringmoons_ring), intent(inout)    :: ring
   type(ringmoons_seeds), intent(inout)   :: seeds
   real(DP), intent(in)                   :: dt
   logical(LGT),intent(out)               :: stepfail

! Internals
   integer(I4B)                              :: i, j, iRRL, nfz, seed_bin
   real(DP)                                  :: dadt, e, inc, sigavg
   real(DP),dimension(seeds%N)               :: Pa, Qa, Ra, Sa
   real(DP),dimension(seeds%N)               :: Pm, Qm, Rm, Sm
   real(DP),dimension(seeds%N)               :: daseeds,dGmseeds,ai,Gmi,fz_width
   type(ringmoons_ring)                      :: iring
   type(ringmoons_seeds)                     :: iseeds
   real(DP)                                  :: da,Gmleft,dGm,Gmdisk,Lfromring,Lseed_original,Ldiff
   real(DP),dimension(0:ring%N+1)            :: dTorque_ring

   

! Executable code

   e = 0.0_DP
   inc = 0.0_DP
   stepfail = .false.
   
   iring%N = ring%N
   iseeds%N = seeds%N
   call ringmoons_allocate(iring,iseeds)

   ! Save initial state of the seeds
   iring = ring 
   iseeds = seeds
   ai(:) = seeds%a(:)
   Gmi(:) = seeds%Gm(:)
   call ringmoons_update_seeds(swifter_pl1P,iring,iseeds)
   call ringmoons_calc_torques(swifter_pl1P,iring,iseeds)

   seeds%Ttide(:) = iseeds%Ttide(:)
   dTorque_ring(:) = iring%Torque(:)

   Pa(:) = 0.0_DP
   Pm(:) = 0.0_DP
   Qa(:) = 0.0_DP
   Qm(:) = 0.0_DP
   Ra(:) = 0.0_DP
   Rm(:) = 0.0_DP
   Sa(:) = 0.0_DP
   Sm(:) = 0.0_DP

   ! First RK step
   do i = 1, seeds%N
      if (seeds%active(i)) then
         Pa(i) = dt * ringmoons_seed_dadt(swifter_pl1P%mass,iseeds%Gm(i),iseeds%a(i),iseeds%Torque(i))
         nfz = iseeds%fz_bin_outer(i) - iseeds%fz_bin_inner(i) + 1
         sigavg = sum(iring%Gsigma(iseeds%fz_bin_inner(i):iseeds%fz_bin_outer(i))) / real(nfz, kind = DP)
         Pm(i)= dt * ringmoons_seed_dMdt(iring,swifter_pl1P%mass,sigavg,iseeds%Gm(i),iseeds%a(i))
      end if
   end do

   iseeds%a(:)  = ai(:)  + 0.5_DP * Pa(:)
   iseeds%Gm(:) = Gmi(:) + 0.5_DP * Pm(:)
   call ringmoons_update_seeds(swifter_pl1P,iring,iseeds)
   call ringmoons_calc_torques(swifter_pl1P,iring,iseeds)
   seeds%Ttide(:) = seeds%Ttide(:) + 2 * iseeds%Ttide(:)
   dTorque_ring(:) = dTorque_ring(:) + 2 * iring%Torque(:)

   ! Second RK step.
   do i = 1, seeds%N
      if (seeds%active(i)) then
         nfz = iseeds%fz_bin_outer(i) - iseeds%fz_bin_inner(i) + 1
         sigavg = sum(iring%Gsigma(iseeds%fz_bin_inner(i):iseeds%fz_bin_outer(i))) / real(nfz, kind = DP)

         Qa(i) = dt * ringmoons_seed_dadt(swifter_pl1P%mass,iseeds%Gm(i),iseeds%a(i),iseeds%Torque(i))
         Qm(i) = dt * ringmoons_seed_dMdt(iring,swifter_pl1P%mass,sigavg,iseeds%Gm(i),iseeds%a(i))
      end if
   end do

   iseeds%a(:)  = ai(:)  + 0.5_DP * Qa(:)
   iseeds%Gm(:) = Gmi(:) + 0.5_DP * Qm(:)
   call ringmoons_update_seeds(swifter_pl1P,iring,iseeds)
   call ringmoons_calc_torques(swifter_pl1P,iring,iseeds)
   seeds%Ttide(:) = seeds%Ttide(:) + 2 * iseeds%Ttide(:)
   dTorque_ring(:) = dTorque_ring(:) + 2 * iring%Torque(:)

   ! Third RK step.
   do i = 1, seeds%N
      if (seeds%active(i)) then
         nfz = iseeds%fz_bin_outer(i) - iseeds%fz_bin_inner(i) + 1
         sigavg = sum(iring%Gsigma(iseeds%fz_bin_inner(i):iseeds%fz_bin_outer(i))) / real(nfz, kind = DP)

         Ra(i) = dt * ringmoons_seed_dadt(swifter_pl1P%mass,iseeds%Gm(i),iseeds%a(i),iseeds%Torque(i))
         Rm(i) = dt * ringmoons_seed_dMdt(iring,swifter_pl1P%mass,sigavg,iseeds%Gm(i),iseeds%a(i))
      end if
   end do

   iseeds%a(:)  = ai(:)  + Ra(:)
   iseeds%Gm(:) = Gmi(:) + Rm(:)
   call ringmoons_update_seeds(swifter_pl1P,iring,iseeds)
   call ringmoons_calc_torques(swifter_pl1P,iring,iseeds)
   seeds%Ttide(:) = seeds%Ttide(:) + iseeds%Ttide(:)
   dTorque_ring(:) = dTorque_ring(:) + iring%Torque(:)

   ! Fourth RK step.
   do i = 1, seeds%N
      if (seeds%active(i)) then
         nfz = iseeds%fz_bin_outer(i) - iseeds%fz_bin_inner(i) + 1
         sigavg = sum(iring%Gsigma(iseeds%fz_bin_inner(i):iseeds%fz_bin_outer(i))) / real(nfz, kind = DP)

         Sa(i) = dt * ringmoons_seed_dadt(swifter_pl1P%mass,iseeds%Gm(i),iseeds%a(i),iseeds%Torque(i))
         Sm(i) = dt * ringmoons_seed_dMdt(iring,swifter_pl1P%mass,sigavg,iseeds%Gm(i),iseeds%a(i))
      end if
   end do

   daseeds(:)  = (Pa(:) + 2 * Qa(:) + 2 * Ra(:) + Sa(:)) / 6._DP
   dGmseeds(:) = (Pm(:) + 2 * Qm(:) + 2 * Rm(:) + Sm(:)) / 6._DP
   

   call ringmoons_deallocate(iring,iseeds)

   stepfail = .false.
   do i = 1,seeds%N
      if (seeds%active(i).and.((abs(daseeds(i)) / ai(i) > RK_FACTOR))) then
         stepfail = .true.
         return
      end if
   end do

   seeds%Ttide(:) = seeds%Ttide(:) / 6._DP
   ring%Torque(:) = ring%Torque(:) + dTorque_ring(:) / 6._DP

   !write(*,*) maxval(daseeds(:),seeds%active(:))
   swifter_pl1P%rot(3) = swifter_pl1P%rot(3) - dt * sum(seeds%Ttide(:),seeds%active(:)) / (swifter_pl1P%mass * swifter_pl1P%Ip(3) * swifter_pl1P%radius**2)
   seeds%Torque(:) = 0.0_DP

   ! Now move and grow the seeds
   ! Get the mass out of the feeding zone (if it's there), starting from the bin the seed is in and working its way outward
   do i = 1, seeds%N
      if (seeds%active(i)) then
         Lseed_original = seeds%Gm(i) * sqrt((swifter_pl1P%mass + seeds%Gm(i)) * seeds%a(i))
         seed_bin = seeds%rbin(i)
         Gmleft = dGmseeds(i)
         Lfromring = 0.0_DP
         nfz = seeds%fz_bin_outer(i) - seeds%fz_bin_inner(i) + 1
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
         da =  ((Lseed_original + Lfromring) / seeds%Gm(i))**2 / swifter_pl1P%mass - seeds%a(i)
         seeds%a(i) = ((Lseed_original + Lfromring) / seeds%Gm(i))**2 / (swifter_pl1P%mass + seeds%Gm(i))
      end if
   end do

   do i = 1, seeds%N
      if (seeds%active(i)) then
         seeds%a(i) = seeds%a(i) + daseeds(i)
         if (seeds%a(i) <= ring%RRL) then   ! Destroy the satellite!
            write(*,*) 'We are on our way to destruction!'
            DESTRUCTION_EVENT = .true.
            DESTRUCTION_COUNTER = 0
            seeds%active(i) = .false.
         end if
      end if
   end do

   !call ringmoons_update_seeds(swifter_pl1P,ring,seeds)

   fz_width(:) = FEEDING_ZONE_FACTOR * seeds%Rhill(:)

   !I'm hungry! What's there to eat?! Look for neighboring seeds
   do i = 1, seeds%N
      if (seeds%active(i)) then
         do j = 1,seeds%N
            if (seeds%active(j).and.(j /= i)) then
               if ((seeds%a(j) > seeds%a(i) - fz_width(i)).and.(seeds%a(j) < seeds%a(i) + fz_width(i))) then ! This one is in the feeding zone
                  ! conserve orbital angular momentum
                  seeds%a(i) = ((seeds%Gm(i) * sqrt(seeds%a(i)) + seeds%Gm(j) * sqrt(seeds%a(j))) &
                                 / (seeds%Gm(i) + seeds%Gm(j)))**2
                  ! conserve mass
                  seeds%Gm(i) = seeds%Gm(i) + seeds%Gm(j)
                  ! deactivate particle for now and position it at the FRL to potentially activate later
                  seeds%Gm(j) = 0.0_DP
                  seeds%active(j) = .false.
                  seeds%a(j) = 0.0_DP
               end if
            end if
         end do
      end if
   end do        
   call ringmoons_update_seeds(swifter_pl1P,ring,seeds)

   return


end subroutine ringmoons_seed_evolve
