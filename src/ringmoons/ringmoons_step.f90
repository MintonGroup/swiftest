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
        implicit none
        class(ringmoons_seed), intent(inout) :: self
        class(ringmoons_cb),   intent(inout) :: cb
        class(ringmoons_ring), intent(inout) :: ring
        class(swiftest_parameters), intent(in) :: param

        integer(I4B)                        :: i,j,seed_bin,inner_outer_sign,nbin,Nactive,rbin
        real(DP)                            :: a, delta_mass, mass_left
        logical                             :: open_space,destructo,spawnbin
        real(DP), parameter                 :: dzone_width = 0.025_DP 
            !! Width of the destruction zone as a fraction of the RRL distance
        integer(I4B)                        :: dzone_inner,dzone_outer ! inner and outer destruction zone bins
        real(DP)                            :: ndz,Lseed_orig,c,b,dr,Lring,deltaL
        real(DP)                            :: mass_min,R_min,dt,dtleft
        logical                             :: stepfail
        class(ringmoons_ring), allocatable   :: seedring
        
        ! First convert any recently destroyed satellites into ring material
        destructo = .false.
        associate(seed => self)
            if (seed%nbody > 0) then
                do i = 1, seed%nbody
                    if (seed%lactive(i)) then
                        if (seed%a(i) <= ring%RRL) then   ! Destroy the satellite!
                            write(*,*) 'We are on our way to destruction!'
                            ! DESTRUCTION_EVENT = .true.
                            ! DESTRUCTION_COUNTER = 0
                            seed%lactive(i) = .false.
                        end if
                    end if
                end do

                do i = 1,seed%nbody
                    if ((.not.seed%lactive(i)).and.(seed%mass(i) > 0.0_DP)) then
                        write(*,*) 'Destruction activated!',i,seed%a(i),seed%mass(i)
                        destructo = .true.
                        Lseed_orig = seed%mass(i) * sqrt((cb%mass + seed%mass(i)) * seed%a(i)) 

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
                        !j = seedring%iRRL
                        ! Offset in angular momentum
                        Lring = sum(seedring%mass(:) * seedring%Iz(:) * seedring%wkep(:))
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
                            call seedring%step(cb,dt,stepfail)
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
                Nactive = count(seed%lactive(:))
                seed%a(1:Nactive) = pack(seed%a(:),seed%lactive(:))
                seed%mass(1:Nactive) = pack(seed%mass(:),seed%lactive(:))
                seed%lactive(1:Nactive) = .true.
                if (size(seed%lactive) > Nactive) seed%lactive(Nactive+1:size(seed%lactive)) = .false.
                seed%nbody = Nactive
                seed%ringbin(1:seed%nbody) = ring%find_bin(seed%a(1:seed%nbody))
            end if
                
            ! Make seed small enough to fit into each bin 
            do i = ring%iFRL,ring%nbins
                spawnbin = .true.
                if (ring%sigma(i)*(param%MU2KG*1e3)/(param%DU2M*1e2)**2/param%GU < 1e-2_DP) spawnbin = .false. 
                    ! no bins that don't have enough mass
                if (seed%nbody > 0) then
                    if (any(seed%ringbin(:) == i .and. seed%lactive(:))) spawnbin = .false. 
                end if
                ! no bins that already have a seed
                ! See Tajeddine et al. (2017) section 2.3. Spawn seed if aggregates are 1% the gap opening mass
                !R_min = 0.01_DP * (3.3e5 / DU2CM) *  (ring%nu(i) / (100 * TU2S / DU2CM**2))
                !if (ring%r_pdisk(i) > R_min) spawnbin = .true.
                if (spawnbin) then
                    a = ring%r(i)
                    delta_mass = ring%m_p(i) * 1
                    call seed%spawn(cb,ring,a,delta_mass,param)
                end if
            end do     
            if (seed%nbody > 0) then
                where (seed%lactive(:))
                    seed%ringbin(:) = ring%find_bin(seed%a(:))
                elsewhere
                    seed%ringbin(:) = 0
                end where
            end if

        end associate
        return
    end subroutine ringmoons_step_restructure_seed

    module subroutine ringmoons_step_ring(self,cb,dt,stepfail)
        implicit none
        class(ringmoons_ring), intent(inout) :: self
        class(ringmoons_cb),   intent(in)    :: cb
        real(DP),              intent(in)    :: dt
        logical,               intent(out)   :: stepfail

        call self%update(cb)
        stepfail = .false.
        return
    end subroutine ringmoons_step_ring

    module subroutine ringmoons_step_seed(self, cb, ring, dt, stepfail)
        implicit none
        class(ringmoons_seed), intent(inout) :: self
        class(ringmoons_cb),   intent(inout) :: cb
        class(ringmoons_ring), intent(inout) :: ring
        real(DP),              intent(in)    :: dt
        logical,               intent(out)   :: stepfail

        stepfail = .false.
        return
    end subroutine ringmoons_step_seed

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
        integer(I4B), parameter             :: LOOPMAX = 2147483646 
        integer(I4B), parameter             :: SUBMAX = 2
        logical                             :: stepfail
        real(DP)                            :: dtring,dtleft
        integer(I4B)                        :: i,loop,seedloop,subcount,loopcount
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
                
                call ring%update(cb)
                
                dtring = ring%get_dt(dtring)

                call cb%accrete(ring,seed,param,dtring)

                call seed%step(cb,ring,dtring,stepfail)

                if (stepfail) then
                    dtring = dtring / SUBMAX
                    subcount = 0
                    call ringmoons_step_restart(self,old_system)
                    cycle
                end if

                call ring%step(cb,dtring,stepfail)
                if (stepfail) then
                    if (dtring > dt * DTMIN_FAC) then
                        dtring = dtring / submax
                        subcount = 0
                        call ringmoons_step_restart(self,old_system)
                        cycle
                    end if
                end if

                call seed%restructure(cb,ring,param) ! Spawn new seed in any available bins outside the FRL 
                                                                            ! where there is ring material
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
            end associate
            end select

            deallocate(old_system%cb, old_system%ring, old_system%seed)
            allocate(old_system%ring, source=self%ring)
            allocate(old_system%seed, source=self%seed)
            allocate(old_system%cb,   source=self%cb)
        end do

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
