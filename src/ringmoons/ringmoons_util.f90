! Copyright 2026 - The Minton Group at Purdue University
! This file is part of Swiftest.
! Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
! Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
! You should have received a copy of the GNU General Public License along with Swiftest. 
! If not, see: https://www.gnu.org/licenses. 

submodule(ringmoons) s_ringmoons_util
    use swiftest
contains

    module subroutine ringmoons_util_accrete_cb(self,ring,seed,param,dt)
        use, intrinsic :: ieee_exceptions
        implicit none

        ! Arguments
        class(ringmoons_cb),        intent(inout) :: self
            !! Ringmoons central body object
        class(ringmoons_ring), intent(inout) :: ring
            !! Ringmoons ring object
        class(ringmoons_seed), intent(inout) :: seed
            !! Ringmoons seed obje3ct
        class(swiftest_parameters), intent(in) :: param
            !! Current run configuration parameters
        real(DP),intent(in)                :: dt
            !! Current time step size

        ! Internals
        integer(I4B) :: i,j,iin
        real(DP) :: Lplanet, Lring, Ltot,Rnew,Mnew, Lorig,Mring,dGMtot,Lnow,rfac
        real(DP),dimension(seed%nbody) :: afac
        real(DP),dimension(0:ring%nbins+1)        :: mtmp,Lring_orig,Lring_now,dL
        real(DP) :: Lp0,Ls0,Lp1,Ls1,Lr0,Lr1,drot0,drot1
        logical, dimension(size(IEEE_ALL))      :: fpe_halting_modes

        ! Guard against underflow errors when rings surface mass density gets too small
        call ieee_get_halting_mode(IEEE_ALL,fpe_halting_modes)
        call ieee_set_halting_mode(ieee_underflow, .false.)
        
        associate(cb => self)
            ring%inside = ring%find_bin(cb%radius)
            dGMtot = sum(param%GU*ring%mass(0:ring%inside))
                    
            !Add ring mass and angular momentum to planet
            Lring_orig(:) = ring%mass(:) * ring%Iz(:) * ring%nkep(:) 
            Lring = sum(Lring_orig(0:ring%inside))
            cb%dGM = cb%dGM + dGMtot 
            cb%dL(3) = cb%dL(3) + Lring
            cb%Gmass = cb%GM0 + cb%dGM 
            cb%mass = cb%Gmass / param%GU
            cb%dR = cb%R0 * (cb%dGM / cb%GM0)**(1.0_DP / 3.0_DP) 
            cb%radius = cb%R0 + cb%dR
            drot0 = cb%L0(3) * RAD2DEG / (cb%Ip(3) * cb%mass * cb%radius**2)  
            drot1 = cb%dL(3) * RAD2DEG / (cb%Ip(3) * cb%mass * cb%radius**2)
            cb%rot(3) = drot0 + drot1

            ring%mass(0:ring%inside) = 0.0_DP
            ring%sigma(0:ring%inside) = 0.0_DP

            if (seed%nbody > 0) then
                afac(1:seed%nbody) = 1._DP - dGMtot / seed%mu(1:seed%nbody)
                seed%a(1:seed%nbody) = seed%a(1:seed%nbody) * afac(1:seed%nbody)
            end if

            ! update body-dependent parameters as needed
            rfac = 1._DP - dGMtot / cb%Gmass
            ring%r_outer = ring%r_outer * rfac
            ring%r_inner = ring%r_inner * rfac

            ! Save the mass so that we can correct for the change in geometry
            mtmp(:) = ring%mass(:)
            call ring%reset(seed,cb,param)
            ring%mass(:) = mtmp(:)
            ring%sigma(:) = ring%mass(:) / ring%deltaA(:)
            ring%Gsigma(:) = param%GU * ring%sigma(:)
            
            ! Any difference in angular momentum in each ring bin will result in a torque in that bin
            Lring_now(:) = ring%mass(:) * ring%Iz(:) * ring%nkep(:) 
            dL(0:ring%inside) = 0.0_DP
            dL(ring%nbins+1) = 0.0_DP
            dL(ring%inside+1:ring%nbins) = (Lring_now(ring%inside+1:ring%nbins) - Lring_orig(ring%inside+1:ring%nbins)) / dt
            ring%Torque(:) = ring%Torque(:) - dL(:) 
        end associate

        call ieee_set_halting_mode(IEEE_ALL, fpe_halting_modes)
        return
    end subroutine ringmoons_util_accrete_cb

    module subroutine ringmoons_util_dealloc_ring(self)
        !! author: David A. Minton
        !!
        !! Deallocates all allocatabale arrays
        implicit none
        ! Arguments
        class(ringmoons_ring),  intent(inout) :: self 
            !! Ringmoons ring object

        self%nbins = 0
        if (allocated(self%r))       deallocate(self%r)
        if (allocated(self%X))       deallocate(self%X)
        if (allocated(self%X2))      deallocate(self%X2)
        if (allocated(self%r_hstar)) deallocate(self%r_hstar)
        if (allocated(self%deltaA))  deallocate(self%deltaA)
        if (allocated(self%sigma))   deallocate(self%sigma)
        if (allocated(self%Gsigma))  deallocate(self%Gsigma)
        if (allocated(self%mass))    deallocate(self%mass)
        if (allocated(self%tau))     deallocate(self%tau)
        if (allocated(self%nu))      deallocate(self%nu)
        if (allocated(self%Q))       deallocate(self%Q)
        if (allocated(self%Iz))      deallocate(self%Iz)
        if (allocated(self%nkep))    deallocate(self%nkep)
        if (allocated(self%Torque))  deallocate(self%Torque)
        if (allocated(self%r_p))     deallocate(self%r_p)
        if (allocated(self%m_p))     deallocate(self%m_p)
        if (allocated(self%rho_p))   deallocate(self%rho_p)
        if (allocated(self%vrel_p))  deallocate(self%vrel_p)

        return
    end subroutine ringmoons_util_dealloc_ring

    module subroutine ringmoons_util_dealloc_seed(self)
        !! author: David A. Minton
        !!
        !! Deallocates all allocatabale arrays
        implicit none
        ! Arguments
        class(ringmoons_seed),  intent(inout) :: self 
            !! Ringmoons ring object

        self%nbody = 0
        if (allocated(self%status))  deallocate(self%status)
        if (allocated(self%id))      deallocate(self%id)
        if (allocated(self%info))    deallocate(self%info)
        if (allocated(self%a))       deallocate(self%a)
        if (allocated(self%mu))      deallocate(self%mu)
        if (allocated(self%mass))    deallocate(self%mass)
        if (allocated(self%Gmass))   deallocate(self%Gmass)
        if (allocated(self%rhill))   deallocate(self%rhill)
        if (allocated(self%radius))  deallocate(self%radius)
        if (allocated(self%density)) deallocate(self%density)
        if (allocated(self%ringbin)) deallocate(self%ringbin)
        if (allocated(self%Torque))  deallocate(self%Torque)
        if (allocated(self%Ttide))   deallocate(self%Ttide)

        return
    end subroutine ringmoons_util_dealloc_seed

    module subroutine ringmoons_util_dealloc_system(self)
        !! author: David A. Minton
        !!
        !! Deallocates all allocatables and resets all values to defaults. Acts as a base for a finalizer
        implicit none
        class(ringmoons_nbody_system), intent(inout) :: self
            !! Ringmoons nbody system object to deallocate
        if (allocated(self%ring)) deallocate(self%ring)
        if (allocated(self%seed)) deallocate(self%seed)
        call self%symba_nbody_system%dealloc()

        return
    end subroutine ringmoons_util_dealloc_system

    pure elemental module function ringmoons_util_find_bin(self,r) result(bin)
        !! author: David A. Minton
        !!
        !! Returns the bin containing radius r from the input ring.
        implicit none
        class(ringmoons_ring), intent(in)      :: self
            !! Ringmoons ring object
        real(DP), intent(in)                   :: r
            !! Radial distance at which to search for the bin
        integer(I4B)                           :: bin
            !! The bin containing radial distance r

        if (r > self%r_outer) then
            bin = self%nbins+1
        else if (r < self%r_inner) then
            bin = 0
        else
            bin = ceiling(2 * (sqrt(r) - sqrt(self%r_inner)) / self%deltaX) 
        end if
        
        return
    end function ringmoons_util_find_bin

    module function ringmoons_util_get_dt_ring(self,dtin) result(dtout)
        !! author: David A. Minton
        !!
        !! Calculates the maximum stable timestep for the surface mass density evolution that is not larger than dtin.
        !!
        !! Adapted from Andrew Hesselbrock's  RING-MOONS Python scripts
        implicit none

        ! Arguments
        class(ringmoons_ring), intent(in)  :: self
        real(DP), intent(in)               :: dtin
        real(DP)                           :: dtout

        ! Internals
        integer(I4B)                           :: i
        real(DP)                               :: sig_max,nu_max
        real(DP)                               :: torque_term
       
        associate(ring => self)
            ! Start with viscous stability
            dtout = dtin
    
            nu_max = max(maxval(ring%nu(:) / ring%X2(:)),0.0_DP)
    
            if (nu_max > 0.0_DP) then
                sig_max = 16 * (12._DP / (ring%deltaX)**2) * nu_max
                dtout = min(dtin,(sig_max)**(-1))
            end if
        end associate
        
        return
    end function ringmoons_util_get_dt_ring

    module subroutine ringmoons_util_update_ring(self,cb,param)
        !! author: David A. Minton
        !!
        !! Updates the ring velocity dispersion, Toomre parameter, and viscosity values
        implicit none
        class(ringmoons_ring), intent(inout) :: self
        class(ringmoons_cb), intent(in) :: cb
        class(swiftest_parameters), intent(in) :: param
        ! Internals
        integer(I4B) :: i

        associate(ring => self)       
            call ring%set_velocity_dispersion(cb,param)
            where(ring%sigma(:) > 1000 * VSMALL) 
                ring%Q(:) = ring%nkep(:) * ring%vrel_p(:) / (3.36_DP * ring%Gsigma(:))
                ring%tau(:) = PI * ring%r_p(:)**2 * ring%sigma(:) / ring%m_p(:)
                ring%nu(:) = ringmoons_viscosity(ring%Gsigma(:), ring%m_p(:), (ring%vrel_p(:))**2, &
                                                ring%r_p(:), ring%r_hstar(:), ring%Q(:), ring%tau(:), ring%nkep(:))
            elsewhere
                ring%Q(:) = huge(1._DP) / 10._DP
                ring%tau(:) = 0.0_DP
                ring%nu(:) = 0.0_DP
            end where
        end associate

        return
    end subroutine ringmoons_util_update_ring

    module subroutine ringmoons_util_reset_ring(self,seed,cb,param)
        !! author: David A. Minton
        !!
        !!  Resets ring torques and recomputes all dimensional quantities, such as ring extent and limits based on the current
        !! surface mass density and central body properties.
        use, intrinsic :: ieee_exceptions
        implicit none
        ! Arguments
        class(ringmoons_ring), intent(inout) :: self
        class(ringmoons_seed), intent(inout) :: seed
        class(ringmoons_cb), intent(in) :: cb
        class(swiftest_parameters), intent(in) :: param
        ! Internals
        integer(I4B)                        :: i, iFRL, iRRL
        real(DP)                            :: Xlo, rho_p
        logical, dimension(size(IEEE_ALL))      :: fpe_halting_modes

        ! Guard against underflow errors when rings surface mass density gets too small
        call ieee_get_halting_mode(IEEE_ALL,fpe_halting_modes)
        call ieee_set_halting_mode(ieee_underflow, .false.)
        associate(ring => self)
            ring%nu = 0.0_DP
            ring%X_inner = 2 * sqrt(ring%r_inner)
            ring%X_outer = 2 * sqrt(ring%r_outer)
            ring%deltaX = (ring%X_outer - ring%X_inner) / ring%nbins
            where (ring%r_p(:) > tiny(1.0_DP))
                ring%rho_p(:) = ring%m_p(:) / ((4.0_DP / 3.0_DP) * PI * ring%r_p(:)**3)
            end where
            ring%iFRL = ring%nbins / 2
            ring%iRRL = ring%nbins / 2
            do i = 1, 10
                if (ring%rho_p(ring%iFRL) > tiny(1.0_DP)) then
                    rho_p = ring%rho_p(ring%iFRL)
                else
                    rho_p = ring%rho_p(ring%nbins / 2)
                end if
                ring%FRL = 2.456_DP * cb%radius * (cb%density / rho_p)**(1._DP / 3._DP)
                if (ring%rho_p(ring%iRRL) > tiny(1.0_DP)) then
                    rho_p = ring%rho_p(ring%iRRL)
                else
                    rho_p = ring%rho_p(ring%nbins / 2)
                end if
                ring%RRL = 1.44_DP  * cb%radius * (cb%density / rho_p)**(1._DP / 3._DP)
                iFRL = ring%find_bin(ring%FRL)
                iRRL = ring%find_bin(ring%RRL)
                if ((iFRL == ring%iFRL) .and. (iRRL == ring%iRRL)) exit
                ring%iFRL = iFRL
                ring%iRRL = iRRL
            end do

            do i = 0,ring%nbins + 1
                ! Set up X coordinate system (see Bath & Pringle 1981)
                Xlo = ring%X_inner + ring%deltaX * (i - 1)
        
                ring%X(i) = Xlo + 0.5_DP * ring%deltaX
            end do
            ring%X2(:) = ring%X(:)**2
            ring%deltaA(:) = 0.25_DP * PI * ring%X(:)**3 * ring%deltaX 
            ring%r(:) = 0.25_DP * (ring%X(:))**2

            ! Specific moment of inertia of the ring bin and ring angular velocity
            ring%Iz(:) = (ring%r(:))**2
            ring%nkep(:) = sqrt(cb%Gmass / ring%r(:)**3)
            ring%mass(:) = ring%sigma(:) * ring%deltaA(:)
            ring%Gsigma(:) = param%GU * ring%sigma(:)
            
            ring%Torque(:) = 0.0_DP

            if (seed%nbody > 0) then
                seed%rhill(:) = seed%a(:) * (seed%mass(:) / cb%mass / 3)**THIRD
                seed%mu(:) = cb%Gmass + seed%Gmass(:)
                where (seed%status(:) == ACTIVE)
                    seed%ringbin(:)   = ring%find_bin(seed%a(:))
                elsewhere
                    seed%ringbin(:)   = 0
                end where
            end if

        end associate
        call ieee_set_halting_mode(IEEE_ALL, fpe_halting_modes)
        return
    end subroutine ringmoons_util_reset_ring

    module subroutine ringmoons_util_setup_ring(self, n, param)
        !! author: David A. Minton
        !!
        !! Constructor for the ringmoons_ring class. 
        !! Allocates space for all bins and initializes all components with a value. 
        implicit none
        class(ringmoons_ring),      intent(inout) :: self  
            !! Ringmoons ring object
        integer(I4B),               intent(in)    :: n     
            !! Number of bins to allocate space for
        class(swiftest_parameters), intent(in)    :: param 
            !! Current run configuration parameters

        if (n < 0) return

        call self%dealloc()
        self%nbins = n

        if (n == 0) return

        allocate(self%r(0:n+1))
        allocate(self%X(0:n+1))
        allocate(self%X2(0:n+1))
        allocate(self%r_hstar(0:n+1))
        allocate(self%deltaA(0:n+1))
        allocate(self%sigma(0:n+1))
        allocate(self%Gsigma(0:n+1))
        allocate(self%mass(0:n+1))
        allocate(self%tau(0:n+1))
        allocate(self%nu(0:n+1))
        allocate(self%Q(0:n+1))
        allocate(self%Iz(0:n+1))
        allocate(self%nkep(0:n+1))
        allocate(self%Torque(0:n+1))
        allocate(self%r_p(0:n+1))
        allocate(self%m_p(0:n+1))
        allocate(self%rho_p(0:n+1))
        allocate(self%vrel_p(0:n+1))

        self%r(:) = 0.0_DP
        self%X(:) = 0.0_DP
        self%X2(:) = 0.0_DP
        self%r_hstar(:) = 0.0_DP
        self%deltaA(:) = 0.0_DP
        self%sigma(:) = 0.0_DP
        self%Gsigma(:) = 0.0_DP
        self%mass(:) = 0.0_DP
        self%tau(:) = 0.0_DP
        self%nu(:) = 0.0_DP
        self%Q(:) = 0.0_DP
        self%Iz(:) = 0.0_DP
        self%nkep(:) = 0.0_DP
        self%Torque(:) = 0.0_DP
        self%r_p(:) = 0.0_DP
        self%m_p(:) = 0.0_DP
        self%rho_p(:) = 0.0_DP
        self%vrel_p(:) = 0.0_DP

        return
    end subroutine ringmoons_util_setup_ring

    module subroutine ringmoons_util_setup_seed(self, n, param)
        !! author: David A. Minton
        !!
        !! Constructor for the ringmoons_seed class. 
        !! Allocates space for all bins and initializes all components with a value. 
        implicit none
        class(ringmoons_seed),     intent(inout) :: self  
            !! Ringmoons seed object
        integer(I4B),               intent(in)    :: n
            !! Number of seeds to allocate space for
        class(swiftest_parameters), intent(in)    :: param 
            !! Current run configuration parameters
        integer(I4B) :: i

        if (n < 0) return
        call self%dealloc()

        self%nbody = n
        if (n == 0) return

        allocate(self%id(n))
        allocate(swiftest_particle_info :: self%info(n))
        allocate(self%status(n))
        allocate(self%a(n))
        allocate(self%mass(n))
        allocate(self%Gmass(n))
        allocate(self%mu(n))
        allocate(self%rhill(n))
        allocate(self%radius(n))
        allocate(self%density(n))
        allocate(self%ringbin(n))
        allocate(self%Torque(n))
        allocate(self%Ttide(n))

        do i = 1, n
            call self%info(i)%set_value(&
                name = "UNNAMED", &
                particle_type = "UNKNOWN", &
                status = "INACTIVE", & 
                origin_type = "UNKNOWN", &
                collision_id = 0, &
                origin_time = -huge(1.0_DP), & 
                origin_rh = [0.0_DP, 0.0_DP, 0.0_DP], &
                origin_vh = [0.0_DP, 0.0_DP, 0.0_DP], &
                discard_time = huge(1.0_DP), & 
                discard_rh = [0.0_DP, 0.0_DP, 0.0_DP], &
                discard_vh = [0.0_DP, 0.0_DP, 0.0_DP], &
                discard_body_id = -1  &
            )
        end do
        self%id(:) = 0
        self%status(:) = INACTIVE
        self%a(:) = 0.0_DP
        self%mass(:) = 0.0_DP
        self%Gmass(:) = 0.0_DP
        self%mu(:) = 0.0_DP
        self%rhill(:) = 0.0_DP
        self%radius(:) = 0.0_DP
        self%density(:) = 0.0_DP
        self%ringbin(:) = 0
        self%Torque(:) = 0.0_DP
        self%Ttide(:) = 0.0_DP

        return
    end subroutine ringmoons_util_setup_seed

    module subroutine ringmoons_util_setup_initialize_system(self, system_history, param)
        !! author: David A. Minton
        !!
        !! Initialize an SyMBA nbody system from files and sets up the encounter and collision structures
        !! 
        implicit none
        ! Arguments
        class(ringmoons_nbody_system), intent(inout) :: self 
            !! SyMBA nbody_system object
        class(swiftest_storage),allocatable, intent(inout) :: system_history 
            !! Stores the system history between output dumps
        class(swiftest_parameters), intent(inout) :: param 
            !! Current run configuration parameters 

        ! Call parent method
        select type(cb => self%cb)
        class is (ringmoons_cb)
        associate(nbody_system => self, ring=>self%ring, seed=>self%seed)
            call symba_util_setup_initialize_system(nbody_system, system_history, param)
            ring%nc%file_name = param%ring_file
            call ring%read_frame(self%t,param)
            call seed%read_frame(self%t,system_history%nc,param)
            call ring%reset(seed,cb,param)
        end associate
        end select

        return
    end subroutine ringmoons_util_setup_initialize_system

    module subroutine ringmoons_util_spawn_seed(self, cb, ring, ring_bin, param)
        implicit none
        ! Arguments
        class(ringmoons_seed),          intent(inout) :: self
        class(ringmoons_ring),          intent(inout) :: ring
        class(ringmoons_cb),            intent(in)    :: cb
        integer(I4B),                   intent(in)    :: ring_bin
        class(swiftest_parameters),     intent(in)    :: param
        ! Internals
        integer(I4B)        :: i,j,seed_bin,nfz
        type(ringmoons_seed)  :: new_seed
        character(*), parameter :: SEEDFMT = '("Seed",I0.7)'
        character(NAMELEN) :: newname

        associate(seed => self, maxid => self%maxid)
            seed_bin = seed%nbody + 1 
            call new_seed%setup(seed_bin, param)
            if (seed%nbody > 0) then
                ! Copy over the old seed properties to the new
                new_seed%status(1:seed%nbody) = seed%status(1:seed%nbody)
                new_seed%id(1:seed%nbody) = seed%id(1:seed%nbody)
                new_seed%info(1:seed%nbody) = seed%info(1:seed%nbody)
                new_seed%a(1:seed%nbody) = seed%a(1:seed%nbody)
                new_seed%mu(1:seed%nbody) = seed%mu(1:seed%nbody)
                new_seed%mass(1:seed%nbody) = seed%mass(1:seed%nbody)
                new_seed%Gmass(1:seed%nbody) = seed%Gmass(1:seed%nbody)
                new_seed%density(1:seed%nbody) = seed%density(1:seed%nbody)
                new_seed%radius(1:seed%nbody) = seed%radius(1:seed%nbody)
                new_seed%rhill(1:seed%nbody) = seed%rhill(1:seed%nbody)
                new_seed%ringbin(1:seed%nbody) = seed%ringbin(1:seed%nbody)
                new_seed%Torque(1:seed%nbody) = seed%Torque(1:seed%nbody)
                new_seed%Ttide(1:seed%nbody) = seed%Ttide(1:seed%nbody)
                new_seed%ringbin(1:seed%nbody) = seed%ringbin(1:seed%nbody)
            end if
            call seed%setup(seed_bin, param)
            seed%id = new_seed%id
            seed%status = new_seed%status
            seed%info = new_seed%info
            seed%a = new_seed%a
            seed%mass = new_seed%mass
            seed%Gmass = new_seed%Gmass
            seed%mu = new_seed%mu
            seed%density = new_seed%density
            seed%radius = new_seed%radius
            seed%rhill = new_seed%rhill
            seed%ringbin = new_seed%ringbin
            seed%Torque = new_seed%Torque
            seed%Ttide = new_seed%Ttide
            seed%ringbin = new_seed%ringbin

            i = seed_bin 
            j = ring_bin
            seed%ringbin(i) = ring_bin
            seed%id(i) = -1 ! Assign a temporary id to the seed and only update it it survives long enough to be recorded
            seed%status(i) = ACTIVE
            seed%mass(i) = min(ring%mass(j),ring%m_p(j))
            seed%Gmass(i) = param%GU * seed%mass(i)
            seed%mu(i) = cb%Gmass + seed%Gmass(i)
            seed%a(i) = (ring%Iz(j) * ring%nkep(j))**2 / seed%mu(i)
            seed%density(i) = ring%rho_p(j)
            seed%radius(i) = (3 * seed%mass(i) / (4 * PI * seed%density(i)))**(1._DP / 3._DP)
            seed%rhill(i) = seed%a(i) * (seed%mass(i) / cb%mass / 3)**THIRD 
            write(newname, SEEDFMT) seed%id(i)
            call seed%info(i)%set_value(name=newname, particle_type=SEED_TYPE_NAME, status="ACTIVE")

            ! Take away the mass from the ring
            ring%mass(j) = ring%mass(j) - seed%mass(i)
            ring%sigma(j) = ring%mass(j) / ring%deltaA(j)
            ring%Gsigma(j) = param%GU * ring%sigma(j)
            seed%ringbin(i) = ring%find_bin(seed%a(i))

        end associate

        return
    end subroutine ringmoons_util_spawn_seed

    module subroutine ringmoons_util_velocity_dispersion_ring(self,cb,param)
        implicit none
        class(ringmoons_ring), intent(inout) :: self
        class(ringmoons_cb),   intent(in)    :: cb
        class(swiftest_parameters), intent(in) :: param
        ! Internals
        integer(I4B)                         :: i
        real(DP),dimension(0:self%nbins+1)   :: kappa_rhstar,eta_rhstar

        associate(ring => self)
            where(ring%mass(:) > N_DISK_FACTOR * ring%m_p(:))
                ring%r_hstar(:) = ring%r(:) * (2 * ring%m_p(:) /(3._DP * cb%mass))**(1._DP/3._DP) / (2 * ring%r_p(:)) 
                ! See Salmon et al. 2010 for this
                kappa_rhstar(:) = ringmoons_transition_function(ring%r_hstar(:))
                eta_rhstar(:) = 1._DP - kappa_rhstar(:)
                ring%vrel_p(:) = kappa_rhstar(:) * sqrt(param%GU * ring%m_p(:) / ring%r_p(:)) + eta_rhstar(:) * &
                                                    (2 * ring%r_p(:) * ring%nkep(:))
            end where
        end associate

        return
    end subroutine ringmoons_util_velocity_dispersion_ring

    elemental pure function ringmoons_viscosity(Gsigma, m_p, v2_p, r_p, r_hstar, Q, tau, w) result(nu)
        ! Arguments
        real(DP),intent(in) :: Gsigma, m_p, v2_p, r_p, r_hstar, Q, tau, w
        real(DP) :: nu
        ! Internals
        real(DP)       :: kappa_Q,eta_Q,y
        real(DP)       :: nu_trans_stable,nu_grav_stable,nu_trans_unstable,nu_grav_unstable
        real(DP)       :: nu_trans,nu_grav,nu_coll

        nu_trans_unstable = 13 * r_hstar**5 * Gsigma**2 / w**3
        nu_trans_stable = (v2_p / (2 * w)) * (0.46_DP * tau / (1._DP + tau**2))

        nu_grav_stable = 0.0_DP
        nu_grav_unstable = nu_trans_unstable

        y = 0.25_DP * Q 
        kappa_Q = ringmoons_transition_function(y)
        eta_Q = 1._DP - kappa_Q

        nu_trans = kappa_Q * nu_trans_stable + eta_Q * nu_trans_unstable
        nu_grav  = kappa_Q * nu_grav_stable  + eta_Q * nu_grav_unstable
        nu_coll = r_p**2 * w * tau

        nu = nu_trans + nu_grav + nu_coll
        return
    end function ringmoons_viscosity


    elemental pure module function ringmoons_transition_function(yin) result(kappa)
        implicit none
        ! Arguments
        real(DP),intent(in) ::yin
        real(DP) :: kappa
        ! Internals
        real(DP) :: y
        y = min(max(yin,epsilon(1._DP)),1.0_DP-epsilon(1._DP))
        kappa =  0.5_DP * (1._DP + tanh((2 * y - 1._DP) / (y * (1._DP - y)))) 
        return

    end function ringmoons_transition_function

end submodule s_ringmoons_util