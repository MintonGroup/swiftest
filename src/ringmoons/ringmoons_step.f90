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
        real(DP),     parameter            :: DTMIN_FAC = 1e-16_DP
        integer(I4B), parameter            :: LOOPMAX = 2147483646 
        integer(I4B), parameter            :: SUBMAX = 2
        logical                            :: stepfail
        real(DP)                           :: dtring,dtleft
        integer(I4B)                       :: i,loop,seedloop,subcount,loopcount
        type(ringmoons_ring), allocatable  :: old_ring
        type(ringmoons_pl),   allocatable  :: old_seed
        type(ringmoons_cb),   allocatable  :: old_cb 

        select type(cb => self%cb)
        class is (ringmoons_cb)
            ! Step ring
            subcount = 0
            dtleft = dt
            dtring = dtleft
            allocate(old_ring, source=self%ring)
            allocate(old_seed, source=self%seed) 
            allocate(old_cb, source=cb)
            do loop = 1, LOOPMAX
                if (loop == LOOPMAX) then
                    write(*,*) 'LOOPMAX reached in seed evolution. Ringmoons_step failed'
                    call util_exit(FAILURE)
                end if

                self%ring%Torque(:) = 0.0_DP
                self%seed%Torque(:) = 0.0_DP
                
                call ringmoons_update_ring(cb,self%ring)
                
                dtring = ringmoons_ring_timestep(self%ring,dtring)

                call ringmoons_cb_accrete(cb,self%ring,self%seed,dtring)

                call ringmoons_seed_evolve(cb,self%ring,self%seed,dtring,stepfail)

                if (stepfail) then
                    dtring = dtring / SUBMAX
                    subcount = 0
                    call move_alloc(old_ring, self%ring)
                    call move_alloc(old_seed, self%seed)
                    cb%mass = old_cb%mass
                    cb%Gmass = old_cb%Gmass
                    cb%radius = old_cb%radius
                    cb%rot(:) = old_cb%rot(:)
                    cycle
                end if

                call ringmoons_update_ring(cb,self%ring)

                call ringmoons_sigma_solver(self%ring,cb%mass,dtring,stepfail)
                if (stepfail) then
                    if (dtring > dt * DTMIN_FAC) then
                        dtring = dtring / submax
                        subcount = 0
                        call move_alloc(old_ring, self%ring)
                        call move_alloc(old_seed, self%seed)
                        cb%mass = old_cb%mass
                        cb%Gmass = old_cb%Gmass
                        cb%radius = old_cb%radius
                        cb%rot(:) = old_cb%rot(:)
                        cycle
                    end if
                end if

        !write(*,*) 'seed_construct'

                call ringmoons_seed_construct(cb,self%ring,self%seed) ! Spawn new seed in any available bins outside the FRL 
                                                                            ! where there is ring material

                subcount = subcount + 1
                !if (DESTRUCTION_EVENT) then
                !   call ringmoons_io_write_frame(t + (dtin - dtleft - dtleft), ring, seed, ring_outfile, out_stat = "APPEND")
                !   DESTRUCTION_COUNTER = DESTRUCTION_COUNTER + 1
                !   if (DESTRUCTION_COUNTER > DESTRUCTION_SAVE_FRAMES) then
                !      DESTRUCTION_EVENT = .false.
                !      DESTRUCTION_COUNTER = 0
                !   end if
                !end if

                dtleft = dtleft - dtring
                ! Scale the change in the ring torques by the step size reduction in order to get the time-averaged Torque
                if (dtleft <= 0.0_DP) exit
                loopcount = loop
                if (subcount == 2 * submax) then
                    dtring = min(dtleft,submax * dtring)
                    subcount = 0
                end if
                dtring = min(dtleft,dtring)

                deallocate(old_ring, old_seed, old_cb)
                allocate(old_ring, source=self%ring)
                allocate(old_seed, source=self%seed) 
                allocate(old_cb, source=cb)
            end do
        end select

        ! Step the nbody system like normal
        call symba_step_system(self, param, t, dt)

        return
    end subroutine ringmoons_step_system


    function ringmoons_ring_timestep(ring,dtin) result(dtout)
        !! author: David A. Minton
        !!
        !! Calculates the maximum stable timestep for the surface mass density evolution
        !!
        !! Adapted from Andrew Hesselbrock's  RING-MOONS Python scripts
        implicit none

        ! Arguments
        type(ringmoons_ring), intent(in)  :: ring
        real(DP), intent(in)              :: dtin
        real(DP)                          :: dtout

        ! Internals
        integer(I4B)                           :: i
        real(DP)                               :: sig_max,nu_max
        real(DP)                               :: torque_term
        
        ! Start with viscous stability
        dtout = dtin

        nu_max = max(maxval(ring%nu(:) / ring%X2(:)),0.0_DP)


        if (nu_max > 0.0_DP) then
            sig_max = 16 * (12._DP / (ring%deltaX)**2) * nu_max
            dtout = min(dtin,(sig_max)**(-1))
        end if
        
        return
    end function ringmoons_ring_timestep



    subroutine ringmoons_cb_accrete(cb,ring,seed,dt)
        implicit none

    ! Arguments
        type(ringmoons_cb),  intent(inout) :: cb
        type(ringmoons_ring),intent(inout) :: ring
        type(ringmoons_pl),  intent(inout) :: seed
        real(DP),intent(in)                :: dt

    ! Internals
        integer(I4B) :: i,j,iin
        real(DP) :: rlo,rhi,GMP, RP,rhoP, Mratio, Rratio, Mratiosqrt,MratioHill,rfac
        real(DP) :: Lplanet, Lring, Ltot,Rnew,Mnew, Lorig,Mring,dMtot,Lnow
        real(DP),dimension(seed%nbody) :: afac
        real(DP),dimension(0:ring%nbins+1)        :: Gmtmp,Lring_orig,Lring_now,dL
        real(DP) :: Lp0,Ls0,Lp1,Ls1,Lr0,Lr1
        
    ! Executable code
        ring%inside = ringmoons_ring_bin_finder(ring,cb%radius)
        GMP = cb%mass_init + cb%mass_accreted
        dMtot = sum(ring%mass(0:ring%inside))
                
        !Add ring mass and angular momentum to planet
        Lring_orig(:) = ring%mass(:) * ring%Iz(:) * ring%w(:) 
        Lring = sum(Lring_orig(0:ring%inside))
        ring%dGMP = ring%dGMP + dMtot
        ring%dLP = ring%dLP + Lring
        swifter_pl1P%mass = ring%GMPi + ring%dGMP 
        swifter_pl1P%radius = ring%RPi * (1._DP + ring%dGMP / ring%GMPi)**(1.0_DP / 3.0_DP)
        swifter_pl1P%rot(3) = (ring%LPI + ring%dLP) / (swifter_pl1P%Ip(3) * swifter_pl1P%mass * (swifter_pl1P%radius)**2)

        ring%Gm(0:ring%inside) = 0.0_DP
        ring%Gsigma(0:ring%inside) = 0.0_DP

        ! update body-dependent parameters as needed
        rfac = 1._DP - dMtot / GMP
        afac(1:seed%N) = 1._DP - dMtot / (GMP + seed%Gm(1:seed%N))
        seed%a(1:seed%N) = seed%a(1:seed%N) * afac(1:seed%N)
        ring%r_F = ring%r_F * rfac
        ring%r_I = ring%r_I * rfac

        ! Save the mass so that we can correct for the change in geometry
        Gmtmp(:) = ring%Gm
        call ringmoons_ring_construct(swifter_pl1P,ring,seed)
        ring%Gm(:) = Gmtmp(:)
        ring%Gsigma(:) = ring%Gm(:) / ring%deltaA(:)
        
        ! Any difference in angular momentum in each ring bin will result in a torque in that bin
        Lring_now(:) = ring%Gm(:) * ring%Iz(:) * ring%w(:) 
        dL(0:ring%inside) = 0.0_DP
        dL(ring%N+1) = 0.0_DP
        dL(ring%inside+1:ring%N) = (Lring_now(ring%inside+1:ring%N) - Lring_orig(ring%inside+1:ring%N)) / dt 
        ring%Torque(:) = ring%Torque(:) - dL(:)

        return

    end subroutine ringmoons_cb_accrete



end submodule s_ringmoons_step
