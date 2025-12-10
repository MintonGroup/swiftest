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
        if (allocated(self%r))          deallocate(self%r)
        if (allocated(self%X))          deallocate(self%X)
        if (allocated(self%X2))         deallocate(self%X2)
        if (allocated(self%r_hstar))    deallocate(self%r_hstar)
        if (allocated(self%deltaA))     deallocate(self%deltaA)
        if (allocated(self%sigma))      deallocate(self%sigma)
        if (allocated(self%tau))        deallocate(self%tau)
        if (allocated(self%nu))         deallocate(self%nu)
        if (allocated(self%Q))          deallocate(self%Q)
        if (allocated(self%Iz))         deallocate(self%Iz)
        if (allocated(self%w))          deallocate(self%w)
        if (allocated(self%Torque))     deallocate(self%Torque)
        if (allocated(self%r_p))       deallocate(self%r_p)
        if (allocated(self%m_p))       deallocate(self%m_p)
        if (allocated(self%rho_p))     deallocate(self%rho_p)
        if (allocated(self%vrel_p))    deallocate(self%vrel_p)

        return
    end subroutine ringmoons_util_dealloc_ring

    module subroutine ringmoons_util_dealloc_seeds(self)
        !! author: David A. Minton
        !!
        !! Deallocates all allocatabale arrays
        implicit none
        ! Arguments
        class(ringmoons_seeds),  intent(inout) :: self 
            !! Ringmoons ring object

        self%nbody = 0
        if (allocated(self%active))    deallocate(self%active)
        if (allocated(self%a))         deallocate(self%a)
        if (allocated(self%mass))      deallocate(self%mass)
        if (allocated(self%Rhill))     deallocate(self%Rhill)
        if (allocated(self%rbin))      deallocate(self%rbin)
        if (allocated(self%Torque))    deallocate(self%Torque)
        if (allocated(self%Ttide))     deallocate(self%Ttide)

        return
    end subroutine ringmoons_util_dealloc_seeds


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
        allocate(self%w(n))
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
        self%w(:) = 0.0_DP
        self%Torque(:) = 0.0_DP
        self%r_p(:) = 0.0_DP
        self%m_p(:) = 0.0_DP
        self%rho_p(:) = 0.0_DP
        self%vrel_p(:) = 0.0_DP

        return
    end subroutine ringmoons_util_setup_ring

    module subroutine ringmoons_util_setup_seeds(self, n, param)
        !! author: David A. Minton
        !!
        !! Constructor for the ringmoons_seeds class. 
        !! Allocates space for all bins and initializes all components with a value. 
        implicit none
        class(ringmoons_seeds),     intent(inout) :: self  
            !! Ringmoons seeds object
        integer(I4B),               intent(in)    :: n     
            !! Number of bins to allocate space for
        class(swiftest_parameters), intent(in)    :: param 
            !! Current run configuration parameters

        if (n < 0) return

        call self%dealloc()
        self%nbody = n

        if (n == 0) return

        allocate(self%active(n))
        allocate(self%a(n))
        allocate(self%mass(n))
        allocate(self%Rhill(n))
        allocate(self%rbin(n))
        allocate(self%Torque(n))
        allocate(self%Ttide(n))

        self%active(:) = .false.
        self%a(:) = 0.0_DP
        self%mass(:) = 0.0_DP
        self%Rhill(:) = 0.0_DP
        self%rbin(:) = 0
        self%Torque(:) = 0.0_DP
        self%Ttide(:) = 0.0_DP

        return
    end subroutine ringmoons_util_setup_seeds

end submodule s_ringmoons_util