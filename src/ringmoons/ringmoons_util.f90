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
        real(DP) :: rlo,rhi,GMcb, Mcb, Rcb,rho_cb, Mratio, Rratio, Mratiosqrt,MratioHill,rfac
        real(DP) :: Lplanet, Lring, Ltot,Rnew,Mnew, Lorig,Mring,dMtot,Lnow
        real(DP),dimension(seed%nbody) :: afac
        real(DP),dimension(0:ring%nbins+1)        :: mtmp,Lring_orig,Lring_now,dL
        real(DP) :: Lp0,Ls0,Lp1,Ls1,Lr0,Lr1
        
        associate(cb => self)
            ring%inside = ring%find_bin(cb%radius)
            GMcb = cb%GM0 + cb%dGM
            Mcb = GMcb / param%GU
            dMtot = sum(ring%mass(0:ring%inside))
                    
            !Add ring mass and angular momentum to planet
            Lring_orig(:) = ring%mass(:) * ring%Iz(:) * ring%wkep(:) 
            Lring = sum(Lring_orig(0:ring%inside))
            cb%dGM = cb%dGM + dMtot * param%GU
            cb%dL(3) = cb%dL(3) + Lring
            cb%Gmass = cb%GM0 + cb%dGM 
            cb%radius = cb%R0 * (1._DP + cb%dGM / cb%GM0)**(1.0_DP / 3.0_DP)
            cb%rot(3) = (cb%L0(3) + cb%dL(3)) / (cb%Ip(3) * cb%mass * (cb%radius)**2)

            ring%mass(0:ring%inside) = 0.0_DP
            ring%sigma(0:ring%inside) = 0.0_DP

            if (seed%nbody > 0) then
                afac(1:seed%nbody) = 1._DP - dMtot / (Mcb+ seed%mass(1:seed%nbody))
                seed%a(1:seed%nbody) = seed%a(1:seed%nbody) * afac(1:seed%nbody)
            end if

            ! update body-dependent parameters as needed
            rfac = 1._DP - dMtot / Mcb
            ring%r_outer = ring%r_outer * rfac
            ring%r_inner = ring%r_inner * rfac

            ! Save the mass so that we can correct for the change in geometry
            mtmp(:) = ring%mass(:)
            call ring%reset(seed,cb)
            ring%mass(:) = mtmp(:)
            ring%sigma(:) = ring%mass(:) / ring%deltaA(:)
            
            ! Any difference in angular momentum in each ring bin will result in a torque in that bin
            Lring_now(:) = ring%mass(:) * ring%Iz(:) * ring%wkep(:) 
            dL(0:ring%inside) = 0.0_DP
            dL(ring%nbins+1) = 0.0_DP
            dL(ring%inside+1:ring%nbins) = (Lring_now(ring%inside+1:ring%nbins) - Lring_orig(ring%inside+1:ring%nbins)) 
            ring%Torque(:) = ring%Torque(:) - dL(:) / dt
        end associate

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
        if (allocated(self%tau))     deallocate(self%tau)
        if (allocated(self%nu))      deallocate(self%nu)
        if (allocated(self%Q))       deallocate(self%Q)
        if (allocated(self%Iz))      deallocate(self%Iz)
        if (allocated(self%wkep))    deallocate(self%wkep)
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
        if (allocated(self%lactive))  deallocate(self%lactive)
        if (allocated(self%a))       deallocate(self%a)
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

    module subroutine ringmoons_util_dealloc_storage(self)
        !! author: David A. Minton
        !!
        !! Resets a storage object by deallocating all items and resetting the frame counter to 0
        use base, only : base_util_dealloc_storage
        implicit none
        ! Arguments
        class(ringmoons_storage), intent(inout) :: self 
            !! Swiftest storage object

        if (allocated(self%nc)) deallocate(self%nc)

        call base_util_dealloc_storage(self)

        return
    end subroutine ringmoons_util_dealloc_storage

    pure elemental module function ringmoons_util_find_bin(self,r) result(bin)
        !! author: David A. Minton
        !!
        !! Returns the bin containing radius r from the input ring
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
        !! Calculates the maximum stable timestep for the surface mass density evolution
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

    module subroutine ringmoons_util_update_ring(self,cb)
        !! author: David A. Minton
        !!
        !! Updates the ring velocity dispersion, Toomre parameter, and viscosity values
        implicit none
        class(ringmoons_ring), intent(inout) :: self
        class(ringmoons_cb), intent(in) :: cb

        return
    end subroutine ringmoons_util_update_ring

    module subroutine ringmoons_util_reset_ring(self,seed,cb)
        !! author: David A. Minton
        !!
        !!  Resets ring torques and recomputes all dimensional quantities, such as ring extent and limits based on the current
        !! surface mass density and central body properties.
        implicit none
        class(ringmoons_ring), intent(inout) :: self
        class(ringmoons_seed), intent(inout) :: seed
        class(ringmoons_cb), intent(in) :: cb

        integer(I4B)                        :: i
        real(DP)                            :: Xlo

        associate(ring => self)
            ring%nu = 0.0_DP
            ring%X_inner = 2 * sqrt(ring%r_inner)
            ring%X_outer = 2 * sqrt(ring%r_outer)
            ring%deltaX = (ring%X_outer - ring%X_inner) / ring%nbins
            ring%rho_p(:) = ring%m_p(:) / ((4.0_DP / 3.0_DP) * PI * ring%r_p(:)**3)
            ring%FRL = 2.456_DP * cb%radius * (cb%density / ring%rho_p(ring%nbins))**(1._DP / 3._DP)
            ring%RRL = 1.44_DP  * cb%radius * (cb%density / ring%rho_p(ring%nbins))**(1._DP / 3._DP)
            ring%iFRL = ring%find_bin(ring%FRL)
            ring%iRRL = ring%find_bin(ring%RRL)

            do i = 0,ring%nbins + 1
                ! Set up X coordinate system (see Bath & Pringle 1981)
                Xlo = ring%X_inner + ring%deltaX * (i - 1)
        
                ring%X(i) = Xlo + 0.5_DP * ring%deltaX
            end do
            ring%X2(:) = ring%X(:)**2

            ! Convert X to r
            ring%r(:) = 0.25_DP * (ring%X(:))**2
                
            ! Factors to convert surface mass density into mass 
            ring%deltaA(:) = 0.25_DP * PI * ring%X(:)**3 * ring%deltaX !2 * PI * deltar * ring%r(i)
            ring%mass(:) = ring%sigma(:) * ring%deltaA(:)
            
            ! Specific moment of inertia of the ring bin
            ring%Iz(:) = (ring%r(:))**2
            ring%wkep(:) = sqrt(cb%Gmass / ring%r(:)**3)

            ring%Torque(:) = 0.0_DP

            where (seed%lactive(:))
                seed%ringbin(:)   = ring%find_bin(seed%a(:))
            elsewhere
                seed%ringbin(:)   = 0
            end where

        end associate
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
        allocate(self%tau(0:n+1))
        allocate(self%nu(0:n+1))
        allocate(self%Q(0:n+1))
        allocate(self%Iz(0:n+1))
        allocate(self%wkep(0:n+1))
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
        self%tau(:) = 0.0_DP
        self%nu(:) = 0.0_DP
        self%Q(:) = 0.0_DP
        self%Iz(:) = 0.0_DP
        self%wkep(:) = 0.0_DP
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

        if (n < 0) return
        call self%dealloc()

        self%nbody = n
        if (n == 0) return

        allocate(self%lactive(n))
        allocate(self%a(n))
        allocate(self%mass(n))
        allocate(self%Gmass(n))
        allocate(self%rhill(n))
        allocate(self%radius(n))
        allocate(self%density(n))
        allocate(self%ringbin(n))
        allocate(self%Torque(n))
        allocate(self%Ttide(n))

        self%lactive(:) = .false.
        self%a(:) = 0.0_DP
        self%mass(:) = 0.0_DP
        self%Gmass(:) = 0.0_DP
        self%rhill(:) = 0.0_DP
        self%radius(:) = 0.0_DP
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
        associate(nbody_system => self)
            call symba_util_setup_initialize_system(nbody_system, system_history, param)
            nbody_system%ring%nc%file_name = param%ring_file
            call nbody_system%ring%read_frame(self%t,param)
            call nbody_system%seed%setup(0, param)
        end associate

        return
    end subroutine ringmoons_util_setup_initialize_system

    module subroutine ringmoons_util_snapshot(self, param, nbody_system, t, arg)
        !! author: David A. Minton
        !!
        !! Takes a minimal snapshot of the state of the system during an ringmoons so that the trajectories
        !! can be played back through the ringmoons
        implicit none
        ! Internals
        class(ringmoons_storage),  intent(inout)        :: self         
            !! Swiftest storage object
        class(swiftest_parameters),    intent(inout)        :: param        
            !! Current run configuration parameters
        class(swiftest_nbody_system),  intent(inout)        :: nbody_system 
            !! Swiftest nbody system object to store
        real(DP),                  intent(in), optional :: t            
            !! Time of snapshot if different from system time
        character(*),              intent(in), optional :: arg          
            !! Optional argument (needed for extended storage type used in collision snapshots)
        ! Arguments

        if (.not.present(t)) then
            write(*,*) "ringmoons_util_snapshot_ringmoons requires `t` to be passed"
            return
        end if

        ! if (.not.present(arg)) then
        !     write(*,*) "ringmoons_util_snapshot_ringmoons requires `arg` to be passed"
        !     return
        ! end if
    end subroutine ringmoons_util_snapshot

    module subroutine ringmoons_util_spawn_seed(self, cb, ring, a, delta_mass, param)
        implicit none
        ! Arguments
        class(ringmoons_seed),          intent(inout) :: self
        class(ringmoons_ring),          intent(inout) :: ring
        class(ringmoons_cb),            intent(in)    :: cb
        real(DP),                       intent(in)    :: a
        real(DP),                       intent(in)    :: delta_mass
        class(swiftest_parameters),     intent(in)    :: param
        ! Internals
        integer(I4B)          :: i,j,seed_bin,nfz
        type(ringmoons_seed)  :: new_seed

        associate(seed => self)
            seed_bin = seed%nbody + 1 
            if (seed_bin > size(seed%lactive)) then 
                ! If no previously generated inactive seed, we'll take advantage of Fortran 2003 automatic allocation 
                ! and tack it on to the end  
                call new_seed%setup(seed_bin, param)
                new_seed%lactive(1:seed%nbody) = seed%lactive(1:seed%nbody)
                new_seed%a(1:seed%nbody) = seed%a(1:seed%nbody)
                new_seed%mass(1:seed%nbody) = seed%mass(1:seed%nbody)
                new_seed%ringbin(1:seed%nbody) = seed%ringbin(1:seed%nbody)
                new_seed%Torque(1:seed%nbody) = seed%Torque(1:seed%nbody)
                new_seed%Ttide(1:seed%nbody) = seed%Ttide(1:seed%nbody)
                call seed%setup(seed_bin, param)
                seed%lactive = new_seed%lactive
                seed%a = new_seed%a
                seed%mass = new_seed%mass
                seed%ringbin = new_seed%ringbin
                seed%Torque = new_seed%Torque
                seed%Ttide = new_seed%Ttide
                seed%nbody = new_seed%nbody
            end if

            i = seed_bin 
            seed%lactive(i) = .true.
            seed%nbody = i
            j = ring%find_bin(a)
            seed%ringbin(i) = j 
            seed%mass(i) = min(delta_mass,ring%mass(j))

            ! Adjust the semimajor axis in order to conserve angular momentum 
            seed%a(i) = (ring%Iz(j) * ring%wkep(j))**2 / (cb%mass + seed%mass(i))

            ! Take away the mass from the ring
            ring%mass(j) = ring%mass(j) - seed%mass(i)
            ring%sigma(j) = ring%mass(j) / ring%deltaA(j)
            seed%ringbin(i) = ring%find_bin(seed%a(i))

        end associate

        return
    end subroutine ringmoons_util_spawn_seed
end submodule s_ringmoons_util