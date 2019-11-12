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
   integer(I4B)                              :: i, j, iRRL, nfz, seed_bin,ilo,ihi, rkn
   real(DP)                                  :: dadt, e, inc, sigavg, sigrem, Tr_evol,Gmsdot, n, Li, Lj
   real(DP),dimension(count(seeds%active))               :: ai,af,Gmi,Gmf,fz_width, Ttide
   type(ringmoons_ring)                      :: iring
   type(ringmoons_seeds)                     :: iseeds
   real(DP)                                  :: da,Gmleft,dGm,Gmdisk
   real(DP),dimension(0:ring%N+1)            :: dTorque_ring,Gmrdot,Gmringi,Gmringf,dGmring
   real(DP),dimension(0:ring%N+1)            :: Tlind,Tring
   real(DP),dimension(2:4),parameter         :: rkh = (/0.5_DP, 0.5_DP, 1._DP/)
   integer(I4B),dimension(4),parameter       :: rkmult = (/1, 2, 2, 1/)
   real(DP),dimension(0:ring%N+1)            :: kr
   real(DP),dimension(count(seeds%active))               :: ka,km

   

! Executable code

   e = 0.0_DP
   inc = 0.0_DP
   stepfail = .false.

   iring%N = ring%N
   iseeds%N = count(seeds%active(:)) 
   call ringmoons_allocate(iring,iseeds)

   ! Save initial state of the seeds
   iring = ring 
   iseeds%a(:) = pack(seeds%a(:),seeds%active(:))
   iseeds%Gm(:) = pack(seeds%Gm(:),seeds%active(:))
   iseeds%Rhill(:) = pack(seeds%Rhill(:),seeds%active(:))
   iseeds%rbin(:) = pack(seeds%rbin(:),seeds%active(:))
   iseeds%fz_bin_inner(:) = pack(seeds%fz_bin_inner(:),seeds%active(:))
   iseeds%fz_bin_outer(:) = pack(seeds%fz_bin_outer(:),seeds%active(:))
   iseeds%Torque(:) = pack(seeds%Torque(:),seeds%active(:))
   iseeds%Ttide(:) = pack(seeds%Ttide(:),seeds%active(:))
   iseeds%active(1:seeds%N) = .true. 

   ai(:) = iseeds%a(:)
   Gmi(:) = iseeds%Gm(:)
   Gmringi(:) = ring%Gm(:)
   dTorque_ring(:) = 0.0_DP
   Ttide(:) = 0.0_DP
   ka(:) = 0._DP
   km(:) = 0._DP
   kr(:) = 0._DP
   af(:) = 0._DP 
   Gmf(:) = 0._DP
   Gmringf(:) = 0._DP

   do rkn = 1,4 ! Runge-Kutta steps 
      if (rkn > 1) then
         iseeds%a(:)  = ai(:)       + rkh(rkn) * ka(:)
         iseeds%Gm(:) = Gmi(:)      + rkh(rkn) * km(:)
         iring%Gm(:)  = Gmringi(:)  + rkh(rkn) * kr(:)
         iring%Gsigma(:) = iring%Gm(:) / iring%deltaA(:)
         kr(:) = 0._DP
      end if
      
      call ringmoons_update_seeds(swifter_pl1P,iring,iseeds)

      !$OMP PARALLEL DO DEFAULT(PRIVATE) SCHEDULE (static) &
      !$OMP SHARED(iseeds,iring,swifter_pl1P,dt,ka,km,e,inc,rkn) &
      !$OMP REDUCTION(+:dTorque_ring,Ttide,af,Gmf,kr)
      do i = 1, iseeds%N
         ihi = min(iseeds%fz_bin_outer(i),iring%N)
         ilo = max(iseeds%fz_bin_inner(i),1)
         nfz = ihi - ilo + 1

         ! Calculate torques
         Tlind(:) = ringmoons_lindblad_torque(swifter_pl1P,iring,iseeds%Gm(i),iseeds%a(i),e,inc)

         n = sqrt((swifter_pl1P%mass + iseeds%Gm(i)) / iseeds%a(i)**3)
         iseeds%Ttide(i) = ringmoons_tidal_torque(swifter_pl1P,iseeds%Gm(i),n,iseeds%a(i),e,inc) 
         iseeds%Torque(i) = iseeds%Ttide(i) - sum(Tlind(:)) 

         sigavg = sum(iring%Gsigma(ilo:ihi)) / real(nfz, kind = DP)
         Gmsdot = ringmoons_seed_dMdt(iring,swifter_pl1P%mass,sigavg,iseeds%Gm(i),iseeds%a(i))
         if (Gmsdot / iseeds%Gm(i) > tiny(1._DP)) then
             
            km(i) = dt * Gmsdot ! Grow the seed

            Gmrdot(ilo:ihi) = -Gmsdot * iring%Gsigma(ilo:ihi) / sigavg ! pull out of the ring proportional to its surface mass density
            kr(ilo:ihi) = kr(ilo:ihi) + dt * Gmrdot(ilo:ihi)  ! Remove mass from the ring

            ! Make sure we conserve angular momentum
            Tr_evol = sum(Gmrdot(ilo:ihi) * iring%Iz(ilo:ihi) * iring%w(ilo:ihi))
            
            iseeds%Torque(i) = iseeds%Torque(i) - Tr_evol
            Gmf(i) = Gmf(i) + rkmult(rkn) * km(i)
         else
            Gmsdot = 0.0_DP
         end if

         ka(i) = dt * ringmoons_seed_dadt(swifter_pl1P%mass,iseeds%Gm(i),iseeds%a(i),iseeds%Torque(i),Gmsdot)
         
         ! Accumulate RK solutions for seeds
         af(i)  = af(i)  + rkmult(rkn) * ka(i)

         ! Save weighted averages of torques for later
         dTorque_ring(:) = dTorque_ring(:) + rkmult(rkn) * Tlind(:)
         Ttide(i) = Ttide(i) + rkmult(rkn) * iseeds%Ttide(i)
      end do
      !$OMP END PARALLEL DO
      ! Accumulate RK solutions for ring

      Gmringf(:) = Gmringf(:) + rkmult(rkn) * kr(:)

   end do
   call ringmoons_deallocate(iring,iseeds)

   af(:) = ai(:) + af(:) / 6._DP

   stepfail = .false.
   if (any(abs(af(:) - ai(:)) / ai(:) > 2 * RK_FACTOR)) then
      !write(*,*) 'Failed the step: Migration too far'
      stepfail = .true.
      return
   end if 

   Gmf(:) = Gmi(:) + Gmf(:) / 6._DP
   Gmringf(:) = Gmringi(:) + Gmringf(:) / 6._DP

   if (any(Gmringf(:) < 0.0_DP)) then
      !write(*,*) 'Failed the step: Negative disk mass'
      stepfail = .true.
      return
   end if 

   seeds%a(:) = af(:)
   seeds%Gm(:) = Gmf(:)
   ring%Gm(:) = Gmringf(:)
   ring%Gsigma(:) = ring%Gm(:) / ring%deltaA(:)

   call ringmoons_update_seeds(swifter_pl1P,ring,seeds)

   !write(*,*) 'calculate total torques'
   seeds%Ttide(:) = Ttide(:) / 6._DP
   ring%Torque(:) = ring%Torque(:) + dTorque_ring(:) / 6._DP

   swifter_pl1P%rot(3) = swifter_pl1P%rot(3) - dt * sum(seeds%Ttide(:),seeds%active(:)) / (swifter_pl1P%mass * swifter_pl1P%Ip(3) * swifter_pl1P%radius**2)
   seeds%Torque(:) = 0.0_DP
   seeds%Ttide(:) = 0.0_DP


   fz_width(:) = 2 * FEEDING_ZONE_FACTOR * seeds%Rhill(:)

   !I'm hungry! What's there to eat?! Look for neighboring seeds
   !write(*,*) 'chomp'
   do i = 1, seeds%N
      if (seeds%active(i)) then
         do j = i+1,seeds%N
            if (seeds%active(j)) then
               if ((seeds%a(j) > seeds%a(i) - fz_width(i)).and.(seeds%a(j) < seeds%a(i) + fz_width(i))) then ! This one is in the feeding zone
                  ! conserve both mass and angular momentum
                  Li = seeds%Gm(i) * sqrt((swifter_pl1P%mass + seeds%Gm(i)) * seeds%a(i))
                  Lj = seeds%Gm(j) * sqrt((swifter_pl1P%mass + seeds%Gm(j)) * seeds%a(j))
                  seeds%Gm(i) = seeds%Gm(i) + seeds%Gm(j)
                  seeds%a(i) = ((Li + Lj) / seeds%Gm(i))**2 / (swifter_pl1P%mass + seeds%Gm(i))

                  ! deactivate particle for now and position it at the FRL to potentially activate later
                  seeds%Gm(j) = 0.0_DP
                  seeds%active(j) = .false.
               end if
           end if
         end do
      end if
   end do        
   !write(*,*) 'update_seeds'
   call ringmoons_update_seeds(swifter_pl1P,ring,seeds)


   do i = 1, seeds%N
      if (seeds%active(i)) then
         if (seeds%a(i) <= ring%RRL) then   ! Destroy the satellite!
            write(*,*) 'We are on our way to destruction!'
            DESTRUCTION_EVENT = .true.
            DESTRUCTION_COUNTER = 0
            seeds%active(i) = .false.
         end if
      end if
   end do

   return


end subroutine ringmoons_seed_evolve
