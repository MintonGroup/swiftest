! Copyright 2026 - The Minton Group at Purdue University
! This file is part of Swiftest.
! Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
! Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
! You should have received a copy of the GNU General Public License along with Swiftest. 
! If not, see: https://www.gnu.org/licenses. 

submodule(ringmoons) s_ringmoons_step
    use swiftest
contains

    module subroutine ringmoons_step_restructure_seed(self,cb,ring,param)
        use, intrinsic :: ieee_exceptions
        implicit none
        class(ringmoons_seed),      intent(inout) :: self
        class(ringmoons_cb),        intent(inout) :: cb
        class(ringmoons_ring),      intent(inout) :: ring
        class(swiftest_parameters), intent(in)    :: param

        integer(I4B)                        :: i,j,inner_outer_sign,nbin,Nactive,rbin
        real(DP)                            :: a, delta_mass, mass_left
        logical                             :: open_space,destructo,spawnbin
        real(DP), parameter                 :: NOSPAWN_THRESHOLD_SIGMA = 0.1_DP
            !! Threshold surface mass density below which seeds cannot spawn in units of kg/m^2
        real(DP), parameter                 :: dzone_width = 0.025_DP 
            !! Width of the destruction zone as a fraction of the RRL distance
        integer(I4B)                        :: dzone_inner,dzone_outer ! inner and outer destruction zone bins
        real(DP)                            :: ndz,Lseed_orig,c,b,dr,Lring,deltaL
        real(DP)                            :: mass_min,R_min,dt,dtleft
        logical                             :: stepfail
        class(ringmoons_ring), allocatable  :: seedring
        logical, dimension(:), allocatable  :: lactive
        logical, dimension(size(IEEE_ALL))      :: fpe_halting_modes

        ! Guard against underflow errors when rings surface mass density gets too small
        call ieee_get_halting_mode(IEEE_ALL,fpe_halting_modes)
        call ieee_set_halting_mode(ieee_underflow, .false.)
        
        ! First convert any recently destroyed satellites into ring material
        destructo = .false.
        associate(seed => self)
            if (seed%nbody > 0) then
                do i = 1, seed%nbody
                    if (seed%status(i) == ACTIVE) then
                        if (seed%a(i) <= ring%RRL) then   ! Destroy the satellite!
                            write(*,*) 'We are on our way to destruction!'
                            ! DESTRUCTION_EVENT = .true.
                            ! DESTRUCTION_COUNTER = 0
                            seed%status(i) = INACTIVE
                        end if
                    end if
                end do

                do i = 1,seed%nbody
                    if ((seed%status(i) == INACTIVE).and.(seed%mass(i) > 0.0_DP)) then
                        destructo = .true.
                        Lseed_orig = seed%mass(i) * sqrt(seed%mu(i) * seed%a(i)) 

                        allocate(seedring, source=ring)
                        seedring%mass(:) = 0.0_DP

                        mass_left = seed%mass(i)
                        c = dzone_width * seed%a(i) ! Create an approximately Gaussian distribution of mass
                        rbin = ring%find_bin(seed%a(i))
                        a = mass_left / (sqrt(2 * PI) * c)
                        do j = 0,(ring%nbins - rbin)
                            do inner_outer_sign = -1,1,2
                                nbin = rbin + inner_outer_sign * j
                                if ((nbin > seedring%inside).and.(nbin < seedring%nbins).and.(mass_left > 0.0_DP)) then
                                    dr = 0.5_DP * seedring%X(nbin) * seedring%deltaX
                                    delta_mass = min(mass_left,a * dr * exp(-(seedring%r(nbin) - seed%a(i))**2 / (2 * c**2)))
                                    seedring%mass(nbin) = seedring%mass(nbin) + delta_mass
                                    mass_left = mass_left - delta_mass
                                end if
                                if (j == 0) exit
                            end do
                            if (mass_left == 0.0_DP) exit
                        end do 
                        
                        ! Offset in angular momentum
                        Lring = sum(seedring%mass(:) * seedring%Iz(:) * seedring%nkep(:))
                        deltaL = Lseed_orig - Lring

                        ! Apply a torque to the temporary ring to bring it back to the seed's original angular momentum
                        dt = 1._DP
                        seedring%nu(1:ring%nbins) =  1._DP / (16 * 12 * dt / (seedring%deltaX)**2) / seedring%X2(1:ring%nbins)
                        seedring%sigma(:) = seedring%mass(:) / seedring%deltaA(:)
                        where (seedring%mass(:) > 0.0_DP)
                            seedring%Torque(:) = deltaL * seedring%mass(:) / sum(seedring%mass(:))
                        elsewhere
                            seedring%Torque(:) = 0.0_DP
                        end where
                        dtleft = dt
                        dt =  seedring%get_dt(dt) 
                        do 
                            call seedring%step(cb,dt,param,stepfail)
                            dtleft = dtleft - dt
                            if (dtleft <= 0.0_DP) exit
                            dt = min(dtleft,dt)
                        end do
                
                        ring%mass(:) = ring%mass(:) + seedring%mass(:) 
                        ring%sigma(:) = ring%mass(:) / ring%deltaA(:)
                        seed%mass(i) = 0.0_DP
                        call seedring%dealloc()
                    end if
                end do
                allocate(lactive(size(seed%status)))
                where(seed%status(:) == ACTIVE)
                    lactive(:) = .true.
                elsewhere
                    lactive(:) = .false.
                end where
                Nactive = count(lactive(:))
                seed%a(1:Nactive) = pack(seed%a(:),lactive(:))
                seed%mass(1:Nactive) = pack(seed%mass(:),lactive(:))
                seed%status(1:Nactive) = ACTIVE
                if (size(seed%status) > Nactive) then
                    seed%status(Nactive+1:size(seed%status)) = INACTIVE
                    seed%ringbin(Nactive+1:size(seed%status)) = -1
                    seed%a(Nactive+1:size(seed%status)) = 0.0_DP
                    seed%mass(Nactive+1:size(seed%status)) = 0.0_DP
                end if
                seed%nbody = Nactive
                seed%ringbin(1:seed%nbody) = ring%find_bin(seed%a(1:seed%nbody))
            end if
                
            ! Check to see if we can spawn seeds beyond the FRL
            do j = ring%iFRL,ring%nbins
                ! skip bins that don't have enough mass
                if (ring%sigma(j) * param%MU2KG/param%DU2M**2 < NOSPAWN_THRESHOLD_SIGMA) cycle
                    
                ! skip bins that already have a seed
                if (seed%nbody > 0) then
                    if (any(seed%ringbin(1:seed%nbody) == j)) cycle
                end if
                call seed%spawn(cb,ring,j,param)
            end do     
            if (seed%nbody > 0) then
                where (seed%status(:) == ACTIVE)
                    seed%ringbin(:) = ring%find_bin(seed%a(:))
                elsewhere
                    seed%ringbin(:) = -1
                end where
            end if

        end associate

        call ieee_set_halting_mode(IEEE_ALL, fpe_halting_modes)

        return
    end subroutine ringmoons_step_restructure_seed

    module subroutine ringmoons_step_ring(self,cb,dt,param,stepfail)
        use, intrinsic :: ieee_exceptions
        implicit none
        ! Arguments
        class(ringmoons_ring),      intent(inout) :: self
        class(ringmoons_cb),        intent(in)    :: cb
        real(DP),                   intent(in)    :: dt
        class(swiftest_parameters), intent(in)    :: param
        logical,                    intent(out)   :: stepfail
        ! Internals
        real(DP),dimension(0:self%nbins+1) :: S,Snew,Sn1,Sn2,fac,artnu,L,dM1,dM2
        integer(I4B)                       :: i,N,j,loop
        logical, dimension(size(IEEE_ALL))      :: fpe_halting_modes
        integer(I4B), parameter :: MAX_LOOP_SOLVER = 100000
        real(DP), parameter :: ARTNU_FAC = 1.0_DP / 16.0_DP 
            !! Artificial viscosity factor to prevent negative mass bins. 

        ! Guard against underflow errors when rings surface mass density gets too small
        call ieee_get_halting_mode(IEEE_ALL,fpe_halting_modes)
        call ieee_set_halting_mode(ieee_underflow, .false.)

        call self%update(cb,param)
        stepfail = .false.

        associate(ring => self, N => self%nbins)
            S(0) = 0.0_DP
            S(1:N) = ring%sigma(1:N) * ring%X(1:N)
            S(N+1) = 0.0_DP
            ring%Torque(0) = 0.0_DP
            ring%Torque(N+1) = 0.0_DP
            ring%mass(0) = 0.0_DP
            ring%mass(N+1) = 0.0_DP
            ring%sigma(0) = 0.0_DP
            ring%sigma(N+1) = 0.0_DP

            fac(:)  = 12 * dt / (ring%deltaX)**2  / ring%X2(:)
            Sn1(1:N) = ring%nu(2:N+1) * S(2:N+1) - 2 * ring%nu(1:N) * S(1:N) + ring%nu(0:N-1) * S(0:N-1)
            Sn2(1:N) = ring%Torque(2:N+1) - ring%Torque(0:N-1)
            Sn2(1:N) = Sn2(1:N) * (1._DP / (3 * PI * sqrt(cb%Gmass)))
            Snew(1:N) = S(1:N) + fac(1:N) * (Sn1(1:N) - Sn2(1:N))

            ! Prevent any bins from having negative mass by diffusing mass with an artificial viscosity
            loop = 1
            do while (any(Snew(1:N) < -epsilon(1._DP) * maxval(S(1:N))))
                Snew(0) = 0.0_DP
                Snew(N+1) = 0.0_DP
                artnu(:) = 0.0_DP
                where (Snew(1:N) < 0._DP)
                    artnu(1:N) = ARTNU_FAC / fac(1:N)
                    artnu(0:N-1) = ARTNU_FAC / fac(0:N-1)
                    artnu(2:N+1) = ARTNU_FAC / fac(2:N+1)
                end where
                S(1:N) = Snew(1:N) 
                Sn1(1:N) = artnu(2:N+1) * S(2:N+1) - 2 * artnu(1:N) * S(1:N) + artnu(0:N-1) * S(0:N-1)
                Snew(1:N) = S(1:N) + fac(1:N) * Sn1(1:N) 
                loop = loop + 1
                if (loop > MAX_LOOP_SOLVER) then
                    stepfail = .true.
                    exit
                end if
            end do
            ring%sigma(1:N) = Snew(1:N) / ring%X(1:N)
            ring%Gsigma(1:N) = param%GU * ring%sigma(1:N)
            ring%mass(1:N) = ring%sigma(1:N) * ring%deltaA(1:N)
            ring%Torque(:) = 0.0_DP
        end associate
        
        call ieee_set_halting_mode(IEEE_ALL, fpe_halting_modes)

        return
    end subroutine ringmoons_step_ring

    module subroutine ringmoons_step_seed(self, cb, ring, dt, param, stepfail)
        !! author: David A. Minton
        !!
        !! Step the seed forward in time and handle any accretion events that occur during the step.
        use, intrinsic :: ieee_exceptions
        implicit none
        ! Arguments
        class(ringmoons_seed),      intent(inout) :: self
        class(ringmoons_cb),        intent(inout) :: cb
        class(ringmoons_ring),      intent(inout) :: ring
        real(DP),                   intent(in)    :: dt
        class(swiftest_parameters), intent(in)    :: param
        logical,                    intent(out)   :: stepfail
        ! Internals
        !Runge-Kutta-Fehlberg parameters
        integer(I4B),parameter                      :: rkfo = 6
        real(DP),dimension(6,5),parameter           :: rkf45_btab = reshape( & ! Butcher tableau for Runge-Kutta-Fehlberg method
            [         1./4.,       1./4.,          0.,            0.,           0.,           0.,&
                        3./8.,      3./32.,      9./32.,            0.,           0.,           0.,&
                    12./13., 1932./2197., -7200./2197.,  7296./2197.,           0.,           0.,&
                        1.,   439./216.,          -8.,   3680./513.,   -845./4104.,          0.,&
                        1./2.,     -8./27.,           2., -3544./2565.,   1859./4104.,    -11./40.], shape(rkf45_btab))
        real(DP),dimension(6),parameter           :: rkf5_coeff =  [ 16./135., 0., 6656./12825., 28561./56430., -9./50., 2./55. ]
        real(DP),dimension(6),parameter           :: rkf4_coeff =  [ 25./216., 0., 1408./2565. ,  2197./4104. , -1./5. ,     0. ]
        integer(I4B)                              :: i, j, iRRL, nfz, ilo,ihi, rkn,rbin, loop, nloops, Ns
        real(DP)                                  :: dadt, sigavg, sigsum, Li, Lj, Ls,dti,dtleft,dtmin
        real(DP)                                  :: impact_b
        class(ringmoons_ring), allocatable        :: iring
        class(ringmoons_seed), allocatable        :: iseed
        real(DP)                                  :: da,mleft,dm,mdisk
        real(DP),dimension(0:ring%nbins+1)        :: dTorque_ring,mringi,mringf,Torquei,Torquef
        real(DP),dimension(0:ring%nbins+1)        :: Tlind,Tring
        real(DP),dimension(0:ring%nbins+1,rkfo)   :: kr,kL
        real(DP),dimension(self%nbody,rkfo)       :: ka,km,kT
        real(DP),dimension(0:ring%nbins+1)        :: Er,rscale,rmdot
        real(DP),dimension(self%nbody)            :: Ea, Em,ascale,mscale
        real(DP),dimension(self%nbody)            :: ai,af,mi,mf, dTtide,Ttidef
        integer(I4B)                              :: Nactive 
        real(DP),dimension(0:ring%nbins+1)        :: Lring_orig,Lring_now
        real(DP),dimension(self%nbody)            :: Lseeds_orig,Lseeds_now,Lres, mdot,adot
        logical, dimension(:), allocatable        :: lactive
        real(DP)                                  :: Lr0,Ls0,Lp0,Lr1,Ls1,Lp1,Lorig,sarr,maxE
        logical                                   :: chomped,goodstep
        real(DP),parameter                        :: DTMIN_FAC = 1e-16_DP
        integer(I4B)                              :: Nnegative_seed,Nnegative_ring,Nbig_error
        logical, dimension(size(IEEE_ALL))      :: fpe_halting_modes
        ! Guard against underflow errors when rings surface mass density gets too small
        call ieee_get_halting_mode(IEEE_ALL,fpe_halting_modes)
        call ieee_set_halting_mode(ieee_underflow, .false.)
        associate(seed => self, Ns => self%nbody)
            ! The first time through we will initialize the laplace coefficients in the torque calculation
            stepfail = .false.
            if (seed%nbody == 0) return
            dti = dt
            dtleft = dt
            dtmin = DTMIN_FAC * dt

            ! Save initial state of the seed
            allocate(iring, source=ring)
            allocate(iseed, source=self)
            ai(1:Ns) = seed%a(1:Ns)
            mi(1:Ns) = seed%mass(1:Ns)
            mringi(:) = ring%mass(:)

            Torquei(:) = ring%Torque(:)
            dTorque_ring(:) = 0.0_DP
            dTtide(:) = 0.0_DP

            Nnegative_seed = 0
            Nnegative_ring = 0
            Nbig_error = 0

            steploop: do loop = 1, LOOPMAX 
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
                    iseed%Torque(:) = 0.0_DP
                    iseed%Ttide(:) = 0.0_DP
                    ! Allow seeds to have negative mass/semimajor axis during intermediate steps, but we'll check at the end to make
                    ! sure that the final result isn't negative
                    iseed%a(1:Ns) = ai(1:Ns) + matmul(ka(1:Ns,1:rkn-1), rkf45_btab(2:rkn,rkn-1))
                    iseed%mass(1:Ns) = mi(1:Ns) + matmul(km(1:Ns,1:rkn-1), rkf45_btab(2:rkn,rkn-1))
                    iring%mass(:)  = mringi(:) + matmul(kr(:,1:rkn-1),rkf45_btab(2:rkn,rkn-1))
                    iring%sigma(:) = iring%mass(:) / iring%deltaA(:)
                    iring%Gsigma(:) = param%GU * iring%sigma(:)
                    iseed%ringbin(1:Ns) = iring%find_bin(seed%a(1:Ns))

                    call iseed%get_tidal_torque(cb,param) 
                    mdot(1:Ns) = ringmoons_dMdt_seed(iseed,iring,cb,param)
                    Tlind(:) = iring%get_lindblad_torque(cb,iseed,param)
                    ! Add the correction factor that comes from conservation of angular momentum between the seed and the ring mass 
                    ! that it consumes
                    where(mdot(1:Ns) > VSMALL)
                        iseed%Torque(1:Ns) = iseed%Torque(1:Ns) + mdot(1:Ns) * iring%Iz(iseed%ringbin(1:Ns)) * iring%nkep(iseed%ringbin(1:Ns))
                    end where
                    adot(:) = ringmoons_dadt_seed(iseed,cb,mdot)
                    ! Set the RKF45 coefficients for this step
                    kr(iseed%ringbin(1:Ns),rkn) = kr(iseed%ringbin(1:Ns),rkn) - dti * mdot(1:Ns) ! Eat the ring
                    km(:,rkn) = dti * mdot(:) ! Grow the seed
                    ka(:,rkn) = dti * adot(:) ! Migrate the seed
                    kT(:,rkn) = dti * iseed%Ttide(:) ! Tidal torque on the seed
                    kL(:,rkn) = kL(:,rkn) + dti * Tlind(:) ! Lindblad torque on the ring
                end do
               
                ! Allow ring mass to go negative, as it will get filled in by sigma_solver
                mringf(:) = mringi(:) + matmul(kr(:,1:rkfo), rkf5_coeff(1:rkfo))
                af(1:Ns) = ai(1:Ns) + matmul(ka(1:Ns,1:rkfo), rkf5_coeff(1:rkfo))

                !Don't let seed semimajor axes go negative
                if (any(af(1:Ns) < 0.0_DP)) then
                    Nnegative_seed = Nnegative_seed + 1 
                    dti = 0.5_DP * dti
                    cycle steploop
                end if

                mf(1:Ns) = mi(1:Ns) + matmul(km(1:Ns,1:rkfo), rkf5_coeff(1:rkfo))

                ! Don't let final masses go netative.
                if (any(mf(1:Ns) < 0.0_DP)) then
                    Nnegative_seed = Nnegative_seed + 1 
                    dti = 0.5_DP * dti
                    cycle steploop
                end if

                ! use the initial value and derivative for error scaling
                ascale(1:Ns) = abs(ai(1:Ns)) + abs(ka(1:Ns,1)) 

                Ea(1:Ns) = abs(matmul(ka(1:Ns,:), (rkf5_coeff(:) - rkf4_coeff(:))))
                maxE = maxval(Ea(1:Ns) / ascale(1:Ns)) / iseed%rkf_tol

                if ((maxE > 1.0_DP).and.(dti > dtmin)) then
                    ! seed a error too high
                    dti = 0.9_DP * dti / maxE**(0.25_DP)
                    goodstep =.false.
                    Nbig_error = Nbig_error + 1
                    cycle steploop
                end if

                mscale(1:Ns) = abs(mi(1:Ns)) + abs(km(1:Ns,1)) 
                Em(1:Ns) = abs(matmul(km(1:Ns,:), (rkf5_coeff(:) - rkf4_coeff(:))))
                maxE = max(maxE, maxval(Em(1:Ns) / mscale(1:Ns)) / iseed%rkf_tol)

                if (maxE > 1.0_DP) then
                    if (dti > dtmin) then
                        ! seed m error too high
                        dti = 0.9_DP * dti / maxE**(0.25_DP)
                        goodstep =.false.
                        Nbig_error = Nbig_error + 1
                        cycle steploop
                    else
                        ! already at minimum step size
                        sarr = 1.0_DP
                    end if
                else if (maxE < 2e-4_DP) then
                    ! error very low
                    sarr = 5._DP
                else
                    ! adjust step size based on error estimate
                    sarr = max(0.90_DP / maxE**(0.25_DP),1._DP)
                end if

                ! save final state of seed and ring and average of torques
                Torquef(:) = matmul(kL(:,1:rkfo), rkf5_coeff(1:rkfo))
                Ttidef(1:Ns) = matmul(kT(1:Ns,1:rkfo), rkf5_coeff(1:rkfo))
                ai(1:Ns) = af(1:Ns)
                mi(1:Ns) = mf(1:Ns)
                mringi(:) = mringf(:)
                dTtide(1:Ns) = dTtide(1:Ns) + Ttidef(1:Ns)
                dTorque_ring(:) = dTorque_ring(:) + Torquef(:)

                dtleft = dtleft - dti
            
                if (dtleft <= 0.0_DP) exit steploop
                dti = min(sarr * dti,dtleft)
            end do steploop

            seed%a(1:Ns) = af(1:Ns)
            seed%mass(1:Ns) = mf(1:Ns)
            seed%Gmass(1:Ns) = param%GU * seed%mass(1:Ns)
            seed%mu(1:Ns) = cb%Gmass + seed%Gmass(1:Ns)
            ring%mass(:) = mringf(:)
            ring%sigma(:) = ring%mass(:) / ring%deltaA(:)
            call ring%reset(seed,cb,param)
            call ring%update(cb,param)
            cb%dL(3) = cb%dL(3) - sum(dTtide(1:Ns))
            ring%Torque(:) = Torquei(:) + dTorque_ring(:) / dt
            
            cb%rot(3) = (cb%L0(3) + cb%dL(3)) / (cb%Ip(3) * cb%mass * cb%radius**2) * RAD2DEG
            seed%Torque(:) = 0.0_DP
            seed%Ttide(:) = 0.0_DP
            seed%ringbin(1:Ns) = ring%find_bin(seed%a(1:Ns))

            !I'm hungry! What's there to eat?! Look for neighboring seed
            chomped = .false.
            do i = 1, Ns
                if (seed%status(i) == ACTIVE) then
                    do j = i + 1, Ns
                        if (seed%status(j) == ACTIVE) then
                        impact_b =  0.5_DP*(seed%a(i)+seed%a(j)) * ((seed%mass(i)+seed%mass(j)) / (3 * cb%mass))**THIRD
                        impact_b = impact_b * seed%feeding_zone_factor 
                        if (abs(seed%a(i) - seed%a(j)) < impact_b) then
                            ! conserve both mass and angular momentum
                            Li = seed%mass(i) * sqrt(seed%mu(i) * seed%a(i))
                            Lj = seed%mass(j) * sqrt(seed%mu(j) * seed%a(j))
                            seed%mass(i) = seed%mass(i) + seed%mass(j)
                            seed%Gmass(i) = param%GU * seed%mass(i)
                            seed%mu(i) = cb%Gmass + seed%Gmass(i)
                            seed%a(i) = ((Li + Lj) / seed%mass(i))**2 / seed%mu(i)
                            seed%density(i) = 0.5_DP * (seed%density(i) + seed%density(j))
                            seed%radius(i) = (3*seed%mass(i) / (4 * PI * seed%density(i)))**THIRD
                            seed%rhill(i) = seed%a(i) * (seed%mass(i) / (3 * cb%mass))**THIRD
                            ! deactivate particle 
                            seed%status(j) = INACTIVE
                            seed%a(j) = 0.0_DP
                            seed%mass(j) = 0.0_DP
                            seed%mu(j) = 0.0_DP
                            seed%Gmass(j) = 0.0_DP
                            seed%radius(j) = 0.0_DP
                            seed%density(j) = 0.0_DP
                            seed%rhill(j) = 0.0_DP
                            chomped = .true.
                        end if
                    end if
                    end do
                end if
            end do      
            if (chomped) then
                allocate(lactive(size(seed%status)))
                lactive(:) = seed%status(:) == ACTIVE
                Nactive = count(lactive(:))
                seed%id(1:Nactive) = pack(seed%id(:),lactive(:))
                seed%a(1:Nactive) = pack(seed%a(:),lactive(:))
                seed%mass(1:Nactive) = pack(seed%mass(:),lactive(:))
                seed%Gmass(1:Nactive) = pack(seed%Gmass(:),lactive(:))
                seed%info(1:Nactive) = pack(seed%info(:), lactive(:))
                seed%status(1:Nactive) = pack(seed%status(:),lactive(:))
                seed%density(1:Nactive) = pack(seed%density(:),lactive(:))
                seed%radius(1:Nactive) = pack(seed%radius(:),lactive(:))
                seed%rhill(1:Nactive) = pack(seed%rhill(:),lactive(:))
                seed%ringbin(1:Nactive) = ring%find_bin(seed%a(1:Nactive))
                seed%nbody = Nactive
                if (size(seed%id) > Nactive) then
                    seed%status(Nactive+1:size(seed%id)) = INACTIVE
                    seed%id(Nactive+1:size(seed%id)) = 0
                    seed%a(Nactive+1:size(seed%id)) = 0.0_DP
                    seed%mass(Nactive+1:size(seed%id)) = 0.0_DP
                    seed%Gmass(Nactive+1:size(seed%id)) = 0.0_DP
                    seed%density(Nactive+1:size(seed%id)) = 0.0_DP
                    seed%rhill(Nactive+1:size(seed%id)) = 0.0_DP
                    seed%radius(Nactive+1:size(seed%id)) = 0.0_DP
                    seed%ringbin(Nactive+1:size(seed%id)) = 0
                end if

            end if
            stepfail = .false.
        end associate

        call ieee_set_halting_mode(IEEE_ALL, fpe_halting_modes)

        return
    end subroutine ringmoons_step_seed

    function ringmoons_dadt_seed(seed,cb,mdot) result(adot)
        use, intrinsic :: ieee_exceptions
        implicit none
        ! Arguments
        class(ringmoons_seed), intent(in) :: seed
        class(ringmoons_cb), intent(in) :: cb
        real(DP), dimension(:), intent(in) :: mdot
        real(DP), dimension(1:seed%nbody)  :: adot
        logical, dimension(size(IEEE_ALL))      :: fpe_halting_modes

        ! Guard against underflow errors when rings surface mass density gets too small
        call ieee_get_halting_mode(IEEE_ALL,fpe_halting_modes)
        call ieee_set_halting_mode(ieee_underflow, .false.)
        associate(Ns => seed%nbody)
            adot(1:Ns) = seed%Torque(1:Ns) / (sqrt(seed%mu(1:Ns) * seed%a(1:Ns))) &
                        - mdot(1:Ns) * (1._DP  + seed%Gmass(1:Ns) / (2 *seed%mu(1:Ns)))
            adot = 2 * adot * seed%a(1:Ns) / seed%mass(1:Ns)
        end associate

        call ieee_set_halting_mode(IEEE_ALL, fpe_halting_modes)
        return
    end function ringmoons_dadt_seed

    function ringmoons_dMdt_seed(seed,ring,cb,param) result(mdot)
        use, intrinsic :: ieee_exceptions
        implicit none
        class(ringmoons_seed),intent(in) :: seed
        class(ringmoons_ring), intent(in) :: ring
        class(ringmoons_cb), intent(in) :: cb 
        class(swiftest_parameters), intent(in) :: param
        real(DP), dimension(1:seed%nbody) :: mdot

        ! Internals
        integer(I4B)                           :: i
        real(DP), dimension(seed%nbody)        :: C
        real(DP),parameter                     :: eff2   = 1e-7_DP ! This term gets the growth rate to match up closely to Andy's
        real(DP),parameter                     :: growth_exp = 4._DP / 3._DP
        real(DP), parameter                    :: SIGLIMIT = 1e-100_DP
        logical, dimension(size(IEEE_ALL))     :: fpe_halting_modes

        ! Guard against underflow errors when rings surface mass density gets too small
        call ieee_get_halting_mode(IEEE_ALL,fpe_halting_modes)
        call ieee_set_halting_mode(ieee_underflow, .false.)

        associate(Ns => seed%nbody)
            ! Executable code
            where((seed%a(1:Ns) > VSMALL).and.ring%mass(seed%ringbin(1:Ns))/seed%mass(1:Ns) > epsilon(1.0_DP))
                C(1:Ns) = 12*PI**(2._DP/3._DP)*(3._DP/(4*seed%density(1:Ns)))**(1._DP/3._DP)*sqrt(param%GU)/sqrt(cb%mass)  
                mdot(1:Ns) = C(1:Ns)*ring%sigma(seed%ringbin(1:Ns))/(eff2*sqrt(seed%a(1:Ns)))*abs(seed%mass(1:Ns))**(growth_exp) 
            elsewhere
                mdot(1:Ns) = 0.0_DP
            end where
        end associate
        
        return
        call ieee_set_halting_mode(IEEE_ALL, fpe_halting_modes)
    end function ringmoons_dMdt_seed

    module subroutine ringmoons_step_system(self, param, t, dt)
        !! author: David A. Minton
        !!
        !! Step the full nbody system forward in time and call the SyMBA step
        !!
        !! Adapted from Andrew Hesselbrock's  RING-MOONS Python scripts
        implicit none
        ! Arguments
        class(ringmoons_nbody_system), intent(inout) :: self   
            !! Ringmoons nbody system object
        class(swiftest_parameters), intent(inout)    :: param  
            !! Current run configuration parameters
        real(DP),                   intent(in)       :: t      
            !! Simulation time
        real(DP),                   intent(in)       :: dt   
        ! Internals
        real(DP),     parameter             :: DTMIN_FAC = 1e-16_DP
        integer(I4B), parameter             :: SUBMAX = 2
        logical                             :: stepfail
        real(DP)                            :: dtring,dtleft
        integer(I4B)                        :: i,loop,subcount,loopcount
        class(ringmoons_nbody_system), allocatable  :: old_system


        allocate(old_system, mold=self)
        subcount = 0
        dtleft = dt
        dtring = dtleft
        
        allocate(old_system%ring, source=self%ring)
        allocate(old_system%seed, source=self%seed)
        allocate(old_system%cb,   source=self%cb)
        do loop = 1, LOOPMAX
            select type(cb => self%cb)
            class is (ringmoons_cb)
            associate(ring => self%ring, seed => self%seed)
                if (loop == LOOPMAX) then
                    write(*,*) 'LOOPMAX reached in seed evolution. Ringmoons_step failed'
                    call base_util_exit(FAILURE)
                end if

                ring%Torque(:) = 0.0_DP
                seed%Torque(:) = 0.0_DP
               
                call ring%update(cb,param)
                
                dtring = ring%get_dt(dtring)

                call cb%accrete(ring,seed,param,dtring)

                call seed%step(cb,ring,dtring,param,stepfail)

                if (stepfail) then
                    dtring = dtring / SUBMAX
                    subcount = 0
                    call ringmoons_step_restart(self,old_system)
                    cycle
                end if

                call ring%step(cb,dtring,param,stepfail)
                if (stepfail) then
                    if (dtring > dt * DTMIN_FAC) then
                        dtring = dtring / submax
                        subcount = 0
                        call ringmoons_step_restart(self,old_system)
                        cycle
                    end if
                end if

                call seed%restructure(cb,ring,param) ! Spawn new seed in any available bins outside the FRL 
                subcount = subcount + 1

                dtleft = dtleft - dtring
                ! Scale the change in the ring torques by the step size reduction in order to get the time-averaged Torque
                if (dtleft <= 0.0_DP) exit
                loopcount = loop
                if (subcount == 2 * submax) then
                    dtring = min(dtleft,submax * dtring)
                    subcount = 0
                end if
                dtring = min(dtleft,dtring)

                deallocate(old_system%cb, old_system%ring, old_system%seed)
                allocate(old_system%ring, source=self%ring)
                allocate(old_system%seed, source=self%seed)
                allocate(old_system%cb,   source=self%cb)
            end associate
            end select
        end do

        self%ring%t = t + dt

        ! Step the nbody system like normal
        call symba_step_system(self, param, t, dt)

        return
    end subroutine ringmoons_step_system

    subroutine ringmoons_step_restart(system, old_system)
        implicit none
        class(ringmoons_nbody_system), intent(inout) :: system
        class(ringmoons_nbody_system), intent(inout) :: old_system

        call move_alloc(old_system%ring, system%ring)
        call move_alloc(old_system%seed, system%seed)
        call move_alloc(old_system%cb, system%cb)
        allocate(old_system%ring, source=system%ring)
        allocate(old_system%seed, source=system%seed)
        allocate(old_system%cb, source=system%cb)

    end subroutine ringmoons_step_restart

end submodule s_ringmoons_step
