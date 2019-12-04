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
   integer(I4B)                              :: i, j, iRRL, nfz, seed_bin,ilo,ihi, rkn,rbin, loop, N,nloops
   real(DP)                                  :: dadt, e, inc, sigavg, sigsum, ns, Tr_evol,Gmsdot, Li, Lj, Ls,dt,dtleft,dtmin
   real(DP)                                  :: Mr0,Ms0,Mr1,Ms1
   type(ringmoons_ring)                      :: iring
   type(ringmoons_seeds)                     :: iseeds
   real(DP)                                  :: da,Gmleft,dGm,Gmdisk
   real(DP),dimension(0:ring%N+1)            :: dTorque_ring,Gmringi,Gmringf,Torquei,Torquef
   real(DP),dimension(0:ring%N+1)            :: Tlind,Tring
   real(DP),dimension(0:ring%N+1,rkfo)       :: kr,kL
   real(DP),dimension(seeds%N,rkfo)          :: ka,km,kT
   real(DP),dimension(0:ring%N+1)            :: Er,rscale,rmdot
   real(DP),dimension(seeds%N)               :: Ea, Em,ascale,mscale,mdot,adot
   real(DP),dimension(seeds%N)               :: ai,af,Gmi,Gmf,fz_width, dTtide,Ttidef
   integer(I4B)                              :: Nactive 
   real(DP),dimension(0:ring%N+1)            :: Lring_orig,Lring_now
   real(DP),dimension(seeds%N)               :: Lseeds_orig,Lseeds_now,Lres
   real(DP)                                  :: Lr0,Ls0,Lp0,Lr1,Ls1,Lp1,Lorig,sarr,Ttide,maxE
   logical(lgt)                              :: chomped,goodstep
   real(DP),parameter                        :: DTMIN_FAC = 1e-4_DP
   real(DP),parameter                        :: TOL = 1e-8_DP 
   integer(I4B)                              :: Nnegative_seed,Nnegative_ring,Nbig_error


! Executable code
   e = 0.0_DP
   inc = 0.0_DP
   stepfail = .false.
   if (seeds%N == 0) return

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

   Nnegative_seed = 0
   Nnegative_ring = 0
   Nbig_error = 0

   steploop: do loop = 1, LOOPMAX 
     ! write(*,*) 'loop ',loop,dt,dtleft,dt/dtin
   !write(*,*) 'negative seed fails: ',Nnegative_seed
   !write(*,*) 'negative ring fails: ',Nnegative_ring
   !write(*,*) 'error too big fails: ',Nbig_error
      nloops = loop
      if (loop == LOOPMAX) then
         stepfail = .true.
         write(*,*) 'max loop reached in seed_evolve'
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
         !write(*,*) 'rkn ',rkn
         iseeds%a(1:iseeds%N) = ai(1:iseeds%N) + matmul(ka(:,1:rkn-1), rkf45_btab(2:rkn,rkn-1))
         if (any(iseeds%a(1:iseeds%N) < 0.0_DP)) then
            !write(*,*) 'negative a'
            Nnegative_seed = Nnegative_seed + 1 
            goodstep = .false.
            dt = 0.5_DP * dt
            cycle steploop
         end if 

         iseeds%Gm(1:iseeds%N) = Gmi(1:iseeds%N) + matmul(km(:,1:rkn-1), rkf45_btab(2:rkn,rkn-1))
         if (any(iseeds%Gm(1:iseeds%N) < 0.0_DP)) then
            !write(*,*) 'negative m'
            Nnegative_seed = Nnegative_seed + 1 
            goodstep = .false.
            dt = 0.5_DP * dt
            cycle steploop
         end if 
         iring%Gm(:)  = Gmringi(:) + matmul(kr(:,1:rkn-1),rkf45_btab(2:rkn,rkn-1))
         if (any(iring%Gm(:) < 0.0_DP)) then
            !write(*,*) 'negative ring m'
            Nnegative_ring = Nnegative_ring + 1 
            goodstep = .false.
            dt = 0.5_DP * dt
            cycle steploop
         end if

         iring%Gsigma(:) = iring%Gm(:) / iring%deltaA(:)

         !write(*,*) 'updating seeds'
         call ringmoons_update_seeds(swifter_pl1P,iring,iseeds)
         !write(*,*) 'updated'

         do i = 1, iseeds%N
          !  write(*,*) 'seed ',i
            rbin = iseeds%rbin(i)

            ! Calculate torques
            Tlind(:) = ringmoons_lindblad_torque(swifter_pl1P,iring,iseeds%Gm(i),iseeds%a(i),e,inc)

            ns = sqrt((swifter_pl1P%mass + iseeds%Gm(i)) / iseeds%a(i)**3)
            Ttide = ringmoons_tidal_torque(swifter_pl1P,iseeds%Gm(i),ns,iseeds%a(i),e,inc) 
            iseeds%Torque(i) = Ttide - sum(Tlind(:)) 

            if ((iring%Gm(rbin) / iseeds%Gm(i)) > epsilon(1._DP))  then
               Gmsdot = ringmoons_seed_dMdt(iring,swifter_pl1P%mass,iring%Gsigma(rbin),iseeds%Gm(i),iring%rho_pdisk(rbin),iseeds%a(i))
               Gmsdot = max(0.0_DP,min(Gmsdot,iring%Gm(rbin) / dt))
               kr(rbin,rkn) = kr(rbin,rkn) - dt * Gmsdot  ! Remove mass from the ring

               ! Make sure we conserve angular momentum during growth
               Tr_evol = Gmsdot * iring%Iz(rbin) * iring%w(rbin)
            else
               Tr_evol = 0.0_DP
               Gmsdot = 0.0_DP
            end if

            km(i,rkn) = dt * Gmsdot ! Grow the seed
            ka(i,rkn) = dt * ringmoons_seed_dadt(swifter_pl1P%mass,iseeds%Gm(i),iseeds%a(i),iseeds%Torque(i) + Tr_evol,Gmsdot)
            kT(i,rkn) = dt * Ttide
            kL(:,rkn) = kL(:,rkn) + dt * Tlind(:)
            !write(*,*) 'all ks calc' 
         end do
      end do
    
      Gmringf(:) = Gmringi(:) + matmul(kr(:,1:rkfo), rkf5_coeff(1:rkfo))
      if (any(Gmringf(:) < 0.0_DP)) then
         Nnegative_ring = Nnegative_ring + 1 
         dt = 0.5_DP * dt
         cycle steploop
      end if

      af(1:iseeds%N) = ai(1:iseeds%N) + matmul(ka(1:iseeds%N,1:rkfo), rkf5_coeff(1:rkfo))

      if (any(af(1:iseeds%N) < 0.0_DP)) then
         Nnegative_seed = Nnegative_seed + 1 
         dt = 0.5_DP * dt
         cycle steploop
      end if

      Gmf(1:iseeds%N) = Gmi(1:iseeds%N) + matmul(km(1:iseeds%N,1:rkfo), rkf5_coeff(1:rkfo))

      if (any(Gmf(1:iseeds%N) < 0.0_DP)) then
         Nnegative_seed = Nnegative_seed + 1 
         dt = 0.5_DP * dt
         cycle steploop
      end if

      ! save values for error estimation
      ascale(:) = abs(ai(:)) + abs(ka(i,1)) + epsilon(1._DP) !1e-3_DP

      Ea(:) = abs(matmul(ka(:,:), (rkf5_coeff(:) - rkf4_coeff(:))))
      maxE = maxval(Ea(1:iseeds%N) / ascale(1:iseeds%N))
      if (maxE > tiny(1._DP)) then
         sarr = (TOL / (2 * maxE))**(0.25_DP) 
      else
         sarr = 2._DP / 0.9_DP
      end if

      mscale(:) = abs(Gmi(:)) + abs(km(i,1)) + epsilon(1._DP) !1e-3_DP
      Em(:) = abs(matmul(km(:,:), (rkf5_coeff(:) - rkf4_coeff(:))))
      maxE = maxval(Em(1:iseeds%N) / mscale(1:iseeds%N))
      if (maxE > tiny(1._DP)) then
         sarr = min(sarr,(TOL / (2 * maxE)))**(0.25_DP) 
      end if

      if ((sarr < 1._DP - epsilon(1._DP))) then
         if (dt > dtmin) then
            goodstep =.false.
            Nbig_error = Nbig_error + 1
            dt = max(0.5_DP * sarr * dt,dtmin)
            cycle steploop
         else
            sarr = 1._DP / 0.90_DP
         end if
      end if

      Torquef(:) = matmul(kL(:,1:rkfo), rkf5_coeff(1:rkfo))
      Ttidef(1:iseeds%N) = matmul(kT(1:seeds%N,1:rkfo), rkf5_coeff(1:rkfo))

      ! save final state of seeds and ring
      ai(1:iseeds%N) = af(1:iseeds%N)
      Gmi(1:iseeds%N) = Gmf(1:iseeds%N)
      Gmringi(:) = Gmringf(:)
      dTtide(1:iseeds%N) = dTtide(1:iseeds%N) + Ttidef(1:seeds%N)
      dTorque_ring(:) = dTorque_ring(:) + Torquef(:)

      dtleft = dtleft - dt
   
      if (dtleft <= 0.0_DP) exit steploop
      dt = min(0.90_DP * sarr * dt,dtleft)

   end do steploop
   
   !write(*,*) 'seed_evolve steploop num: ',nloops
   !write(*,*) 'negative seed fails: ',Nnegative_seed
   !write(*,*) 'negative ring fails: ',Nnegative_ring
   !write(*,*) 'error too big fails: ',Nbig_error

   seeds%a(1:seeds%N) = af(1:seeds%N)
   seeds%Gm(1:seeds%N) = Gmf(1:seeds%N)
   ring%Gm(:) = Gmringf(:)
   ring%Gsigma(:) = ring%Gm(:) / ring%deltaA(:)
   ring%dLP = ring%dLP - sum(dTtide(1:seeds%N))
   ring%Torque(:) = Torquei(:) + dTorque_ring(:) / dtin

   
   do while(any(abs(ring%Torque(2:ring%N)) > 0.0_DP .and. ring%Gm(2:ring%N) < N_DISK_FACTOR * ring%Gm_pdisk(2:ring%N)))
      do i = 2,ring%N
         if (abs(ring%Torque(i)) > 0.0_DP .and. ring%Gm(i) < N_DISK_FACTOR * ring%Gm_pdisk(i)) then
            ring%Torque(i - 1) = ring%Torque(i-1) + ring%Torque(i)
            ring%Torque(i) = 0.0_DP
         end if
      end do
   end do  
      

   swifter_pl1P%rot(3) = (ring%LPi + ring%dLP) / (swifter_pl1P%Ip(3) * swifter_pl1P%mass * (swifter_pl1P%radius)**2) 
   seeds%Torque(:) = 0.0_DP
   seeds%Ttide(:) = 0.0_DP
   call ringmoons_deallocate(iring,iseeds)
   call ringmoons_update_seeds(swifter_pl1P,ring,seeds)

   fz_width(1:seeds%N) = FEEDING_ZONE_FACTOR * seeds%Rhill(1:seeds%N)

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


   stepfail = .false.
   return

end subroutine ringmoons_seed_evolve
