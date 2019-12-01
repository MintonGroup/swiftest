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
subroutine ringmoons_seed_evolve(swifter_pl1P,ring,seeds,dtin,stepfail)

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
   real(DP), intent(in)                   :: dtin
   logical(LGT),intent(out)               :: stepfail

! Internals
   integer(I4B)                              :: i, j, iRRL, nfz, seed_bin,ilo,ihi, rkn,rbin, loop, N
   real(DP)                                  :: dadt, e, inc, sigavg, sigsum, ns, Tr_evol,Gmsdot, Li, Lj, Ls,dt,dtleft,dtmin
   real(DP)                                  :: Mr0,Ms0,Mr1,Ms1
   type(ringmoons_ring)                      :: iring
   type(ringmoons_seeds)                     :: iseeds
   real(DP)                                  :: da,Gmleft,dGm,Gmdisk
   real(DP),dimension(0:ring%N+1)            :: dTorque_ring,Gmringi,Gmringf,Torquei
   real(DP),dimension(0:ring%N+1)            :: Tlind,Tring,L,dM1,dM2
   !real(DP),dimension(2:4),parameter         :: rkh = (/0.5_DP, 0.5_DP, 1._DP/)
   !integer(I4B),dimension(4),parameter       :: rkmult = (/1, 2, 2, 1/)
   real(DP),dimension(0:ring%N+1,rkfo)       :: kr,kL
   real(DP),dimension(seeds%N,rkfo)          :: ka,km,kT
   real(DP),dimension(0:ring%N+1)            :: Evr
   real(DP),dimension(seeds%N)               :: Ea, Em
   real(DP),dimension(seeds%N)               :: ai,af,Gmi,Gmf,fz_width, dTtide, Lseeds
   integer(I4B)                              :: Nactive 
   real(DP),dimension(0:ring%N+1)            :: Lring_orig,Lring_now
   real(DP),dimension(seeds%N)               :: Lseeds_orig,Lseeds_now,Lres
   real(DP)                                  :: Lr0,Ls0,Lp0,Lr1,Ls1,Lp1,Lorig,Lnow,Llind,Ltide,sarr,Tide
   logical(lgt)                              :: chomped,goodstep
   real(DP),parameter                        :: DTMIN_FAC = 1.0e-6_DP
   real(DP),parameter                        :: TOL = 1e-8_DP


! Executable code
   e = 0.0_DP
   inc = 0.0_DP
   stepfail = .false.
   if (seeds%N == 0) return

   Lr0 = sum(ring%Gm(:) * ring%Iz(:) * ring%w(:)) + sum(ring%Torque(:)) * dtin
   Ls0 = sum(seeds%Gm(:) * sqrt((swifter_pl1P%mass + seeds%Gm(:)) * seeds%a(:)),seeds%active(:))
   Lp0 = ring%LPi + ring%dLp !swifter_pl1P%Ip(3) * swifter_pl1P%rot(3) * swifter_pl1P%mass * swifter_pl1P%radius**2
   Lorig = Lr0 + Ls0 + Lp0

   Mr0 = sum(ring%Gm(:))
   Ms0 = sum(seeds%Gm(:),seeds%active(:))
   dtleft = dtin
   dt = dtin
   dtmin = DTMIN_FAC * dtin

   iring%N = ring%N
   iseeds%N = seeds%N
   N = ring%N
   call ringmoons_allocate(iring,iseeds)

   ! Save initial state of the seeds
   iring = ring 
   iseeds = seeds
   ai(1:seeds%N) = seeds%a(1:seeds%N)
   Gmi(1:seeds%N) = seeds%Gm(1:seeds%N)
   Gmringi(:) = ring%Gm(:)

   Torquei(:) = ring%Torque(:)
   dTorque_ring(:) = 0.0_DP
   dTtide(:) = 0.0_DP


   do loop = 1, LOOPMAX 
      if (loop == LOOPMAX) then
         stepfail = .true.
         return
      end if

      stepfail = .false.

      ka(:,:) = 0._DP
      km(:,:) = 0._DP
      kr(:,:) = 0._DP
      kT(:,:) = 0.0_DP
      kL(:,:) = 0.0_DP
      goodstep = .true.

      do rkn = 1,rkfo ! Runge-Kutta steps 

         iseeds%a(1:iseeds%N) = ai(1:iseeds%N) + matmul(ka(:,1:rkn-1), rkf45_btab(2:rkn,rkn-1))
         iseeds%Gm(1:iseeds%N) = Gmi(1:iseeds%N) + matmul(km(:,1:rkn-1), rkf45_btab(2:rkn,rkn-1))
         iring%Gm(:)  = Gmringi(:) + matmul(kr(:,1:rkn-1),rkf45_btab(2:rkn,rkn-1))
         if (any(iring%Gm(:) < 0.0_DP)) then
            !write(*,*) 'negative ring mass in rk step: ',rkn
            goodstep = .false.
            exit
         end if

         iring%Gsigma(:) = iring%Gm(:) / iring%deltaA(:)

         call ringmoons_update_seeds(swifter_pl1P,iring,iseeds)

         do i = 1, iseeds%N
            rbin = iseeds%rbin(i)

            ! Calculate torques
            Tlind(:) = ringmoons_lindblad_torque(swifter_pl1P,iring,iseeds%Gm(i),iseeds%a(i),e,inc)

            ns = sqrt((swifter_pl1P%mass + iseeds%Gm(i)) / iseeds%a(i)**3)
            Tide = ringmoons_tidal_torque(swifter_pl1P,iseeds%Gm(i),ns,iseeds%a(i),e,inc) 
            iseeds%Torque(i) = Tide - sum(Tlind(:)) 

            !if ((iring%Gm(rbin) / iseeds%Gm(i)) > 0.1_DP)  then
               Gmsdot = ringmoons_seed_dMdt(iring,swifter_pl1P%mass,iring%Gsigma(rbin),iseeds%Gm(i),iring%rho_pdisk(rbin),iseeds%a(i))
               Gmsdot = min(Gmsdot,iring%Gm(rbin) / dt)
               kr(rbin,rkn) = kr(rbin,rkn) - dt * Gmsdot  ! Remove mass from the ring
               ! Make sure we conserve angular momentum during growth
               Tr_evol = Gmsdot * iring%Iz(rbin) * iring%w(rbin)
            !else
            !   Tr_evol = 0.0_DP
            !   Gmsdot = 0.0_DP
            !end if

            km(i,rkn) = dt * Gmsdot ! Grow the seed
            ka(i,rkn) = dt * ringmoons_seed_dadt(swifter_pl1P%mass,iseeds%Gm(i),iseeds%a(i),iseeds%Torque(i) + Tr_evol,Gmsdot)
            kT(i,rkn) = dt * Tide
            kL(:,rkn) = kL(:,rkn) +  dt * Tlind(:)
            
         end do
      end do
      if (.not.goodstep) then
         dt = 0.5_DP * dt
         cycle 
      end if
    
      Gmringf(:) = Gmringi(:) + matmul(kr(:,1:rkfo-1), rkf4_coeff(1:rkfo-1))
      ! Prevent any bins from having negative mass by shifting mass upward from interior bins  
      do while (any(Gmringf(1:N) < 0.0_DP))
         where(Gmringf(:) < 0.0_DP)
            dM1(:) = Gmringf(:)
         elsewhere
            dM1(:) = 0.0_DP
         end where
         L(:) = iring%Iz(:) * iring%w(:)
         dM2(:) = dM1(:) * (L(:) - cshift(L(:),1)) / (cshift(L(:),1) - cshift(L(:),2)) 
        ! Make sure we conserve both mass and angular momentum
         Gmringf(1:N) = Gmringf(1:N) - dM1(1:N) + cshift(dM1(1:N),1) + &
            cshift(dM2(1:N),1)  - cshift(dM2(1:N),2)
      end do 
      if (any(Gmringf(:) < 0.0_DP)) then
         dt = 0.5_DP * dt
         cycle 
      end if

      af(1:seeds%N) = ai(1:seeds%N) + matmul(ka(1:seeds%N,1:rkfo-1), rkf4_coeff(1:rkfo-1))
      Gmf(1:seeds%N) = Gmi(1:seeds%N) + matmul(km(1:seeds%N,1:rkfo-1), rkf4_coeff(1:rkfo-1))

      Em(:) = abs(matmul(km(:,:), (rkf5_coeff(:) - rkf4_coeff(:))))
      Ea(:) = abs(matmul(ka(:,:), (rkf5_coeff(:) - rkf4_coeff(:))))
      sarr = (TOL / (2 * max(maxval(Em(:)),maxval(Ea(:)))))**(0.25_DP) 

      if ((sarr < 1._DP).and.(dt > dtmin)) then
         goodstep =.false.
         dt = 0.5_DP * sarr * dt
         cycle
      end if

      dTtide(1:iseeds%N) = dTtide(1:iseeds%N) + matmul(kT(1:iseeds%N,1:rkfo-1), rkf4_coeff(1:rkfo-1))
      dTorque_ring(:) = dTorque_ring(:) + matmul(kL(:,1:rkfo-1), rkf4_coeff(1:rkfo-1)) 

      ! save final state of seeds and ring
      ai(1:iseeds%N) = af(1:iseeds%N)
      Gmi(1:iseeds%N) = Gmf(1:iseeds%N)
      Gmringi(:) = Gmringf(:)

      dtleft = dtleft - dt
   
      if (dtleft <= 0.0_DP) exit
      dt = min(0.9_DP * sarr * dt,dtleft)

   end do

   seeds%a(1:seeds%N) = af(1:seeds%N)
   seeds%Gm(1:seeds%N) = Gmf(1:seeds%N)
   ring%Gm(:) = Gmringf(:)
   ring%Gsigma(:) = ring%Gm(:) / ring%deltaA(:)
   ring%dLP = ring%dLP - sum(dTtide(1:seeds%N))
   ring%Torque(:) = Torquei(:) + dTorque_ring(:) / dtin

   swifter_pl1P%rot(3) = (ring%LPi + ring%dLP) / (swifter_pl1P%Ip(3) * swifter_pl1P%mass * (swifter_pl1P%radius)**2) 
   seeds%Torque(:) = 0.0_DP
   seeds%Ttide(:) = 0.0_DP
   call ringmoons_deallocate(iring,iseeds)
   call ringmoons_update_seeds(swifter_pl1P,ring,seeds)

   Lr1 = sum(ring%Gm(:) * ring%Iz(:) * ring%w(:)) + sum(ring%Torque(:)) * dtin
   Ls1 = sum(seeds%Gm(:) * sqrt((swifter_pl1P%mass + seeds%Gm(:)) * seeds%a(:)),seeds%active(:))
   Lp1 = ring%Lpi + ring%dLp 
   Lnow = Lr1 + Ls1 + Lp1

   if (abs((Lnow - Lorig) / Lorig) > epsilon(1._DP)) then
      stepfail = .true.
      return
   end if 

   fz_width(:) = FEEDING_ZONE_FACTOR * seeds%Rhill(1:seeds%N)

   !I'm hungry! What's there to eat?! Look for neighboring seeds
   !write(*,*) 'chomp'
   chomped = .false.
   do i = 1, seeds%N
      if (seeds%active(i)) then
         do j = i + 1, seeds%N
            if (seeds%active(j)) then
               if ((seeds%a(j) > seeds%a(i) - fz_width(i)).and.(seeds%a(j) < seeds%a(i) + fz_width(i))) then ! This one is in the feeding zone
                  ! conserve both mass and angular momentum
                  
                  Li = seeds%Gm(i) * sqrt((swifter_pl1P%mass + seeds%Gm(i)) * seeds%a(i))
                  Lj = seeds%Gm(j) * sqrt((swifter_pl1P%mass + seeds%Gm(j)) * seeds%a(j))
                  seeds%Gm(i) = seeds%Gm(i) + seeds%Gm(j)
                  seeds%a(i) = ((Li + Lj) / seeds%Gm(i))**2 / (swifter_pl1P%mass + seeds%Gm(i))

                  ! deactivate particle 
                  seeds%Gm(j) = 0.0_DP
                  seeds%a(j) = 0.0_DP
                  seeds%active(j) = .false.
                  chomped = .true.
               end if
           end if
         end do
      end if
   end do      
   if (chomped) then  
      !write(*,*) 'chomped'
      Nactive = count(seeds%active(:))
      seeds%a(1:Nactive) = pack(seeds%a(:),seeds%active(:))
      seeds%Gm(1:Nactive) = pack(seeds%Gm(:),seeds%active(:))
      seeds%active(1:Nactive) = .true.
      if (size(seeds%active) > Nactive) seeds%active(Nactive+1:size(seeds%active)) = .false.
      seeds%N = Nactive
      call ringmoons_update_seeds(swifter_pl1P,ring,seeds)
   end if


       
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

   stepfail = .false.
   return

end subroutine ringmoons_seed_evolve
