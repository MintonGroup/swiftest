! Copyright 2025 - The Minton Group at Purdue University
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
        if (allocated(self%vkep))    deallocate(self%vkep)
        if (allocated(self%Torque))  deallocate(self%Torque)
        if (allocated(self%r_p))     deallocate(self%r_p)
        if (allocated(self%m_p))     deallocate(self%m_p)
        if (allocated(self%rho_p))   deallocate(self%rho_p)
        if (allocated(self%vrel_p))  deallocate(self%vrel_p)

        return
    end subroutine ringmoons_util_dealloc_ring

    module subroutine ringmoons_util_dealloc_pl(self)
        !! author: David A. Minton
        !!
        !! Deallocates all allocatabale arrays
        implicit none
        ! Arguments
        class(ringmoons_pl),  intent(inout) :: self 
            !! Ringmoons ring object

        if (allocated(self%is_seed))      deallocate(self%is_seed)
        if (allocated(self%ringbin))      deallocate(self%ringbin)
        if (allocated(self%Torque))    deallocate(self%Torque)
        if (allocated(self%Ttide))     deallocate(self%Ttide)

        call self%symba_pl%dealloc()

        return
    end subroutine ringmoons_util_dealloc_pl

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

        allocate(self%r(n))
        allocate(self%X(n))
        allocate(self%X2(n))
        allocate(self%r_hstar(n))
        allocate(self%deltaA(n))
        allocate(self%sigma(n))
        allocate(self%tau(n))
        allocate(self%nu(n))
        allocate(self%Q(n))
        allocate(self%Iz(n))
        allocate(self%vkep(n))
        allocate(self%Torque(n))
        allocate(self%r_p(n))
        allocate(self%m_p(n))
        allocate(self%rho_p(n))
        allocate(self%vrel_p(n))

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
        self%vkep(:) = 0.0_DP
        self%Torque(:) = 0.0_DP
        self%r_p(:) = 0.0_DP
        self%m_p(:) = 0.0_DP
        self%rho_p(:) = 0.0_DP
        self%vrel_p(:) = 0.0_DP

        return
    end subroutine ringmoons_util_setup_ring


    module subroutine ringmoons_util_setup_pl(self, n, param)
        !! author: David A. Minton
        !!
        !! Constructor for the ringmoons_pl class. 
        !! Allocates space for all bins and initializes all components with a value. 
        implicit none
        class(ringmoons_pl),     intent(inout) :: self  
            !! Ringmoons seeds object
        integer(I4B),               intent(in)    :: n     
            !! Number of bins to allocate space for
        class(swiftest_parameters), intent(in)    :: param 
            !! Current run configuration parameters

        call self%symba_pl%setup(n, param)
        if (n == 0) return

        allocate(self%is_seed(n))
        allocate(self%ringbin(n))
        allocate(self%Torque(n))
        allocate(self%Ttide(n))

        self%is_seed(:) = .false.
        self%ringbin(:) = 0
        self%Torque(:) = 0.0_DP
        self%Ttide(:) = 0.0_DP

        return
    end subroutine ringmoons_util_setup_pl


    ! module subroutine ringmoons_util_setup_initialize_system(self, system_history, param)
    !     !! author: David A. Minton
    !     !!
    !     !! Initialize a Ringmoons nbody system from files
    !     !!
    !     implicit none
    !     ! Arguments
    !     class(ringmoons_nbody_system),                 intent(inout) :: self            
    !         !! Ringmoons nbody system object
    !     class(swiftest_storage),    allocatable, intent(inout) :: system_history  
    !         !! Stores the system history between output dumps
    !     class(swiftest_parameters),              intent(inout) :: param           
    !         !! Current run configuration parameters 
    
    !     call symba_util_setup_initialize_system(self, system_history, param)

    !   return
    ! end subroutine ringmoons_util_setup_initialize_system

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

end submodule s_ringmoons_util