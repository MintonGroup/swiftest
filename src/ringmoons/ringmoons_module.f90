! Copyright 2025 - The Minton Group at Purdue University
! This file is part of Swiftest.
! Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
! Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
! You should have received a copy of the GNU General Public License along with Swiftest. 
! If not, see: https://www.gnu.org/licenses. 

module ringmoons
    !! author: Andrew Hesselbrock and David A. Minton
    !!
    !! Definition of classes and methods specific to the Ringmoons integrator
    use swiftest
    use symba
    implicit none
    public

    !> Ringmoons central body particle class
    type, extends(symba_cb) :: ringmoons_cb
    end type ringmoons_cb


    !> Ringmoons massive body class
    type, extends(symba_pl) :: ringmoons_pl
    end type ringmoons_pl


    !> Ringmoons test particle class
    type, extends(symba_tp) :: ringmoons_tp
    contains
    end type ringmoons_tp

    type, extends(base_object) :: ringmoons_ring
        integer(I4B) :: nbins               
            !! number of bins in ring
        integer(I4B) :: inside = 1          
            !! bin id of innermost ring bin (can increase if primary accretes a lot mass through updates)
        real(DP)     :: r_F                 
            !! outside radius of ring in simulation length units
        real(DP)     :: X_F                 
            !! outside radius of ring in X units (see Bath & Pringle 1981)
        real(DP)     :: r_I                 
            !! inside radius of ring in simulation length units
        real(DP)     :: X_I                 
            !! inside radius of ring in X units (see Bath & Pringle 1981)
        real(DP)     :: deltaX              
            !! variable changed bin width used for viscosity calculations in X units
        real(DP)     :: RRL,FRL             
            !! Rigid and fluid Roche limits
        integer(I4B) :: iRRL,iFRL           
            !! Indexes of Roche limit bins
        real(DP), dimension(:), allocatable :: r                 
            !! radial distance of center of bin in simulation length units
        real(DP), dimension(:), allocatable :: r_hstar           
            !! normalized ring Hill's radius in simulation length units
        real(DP), dimension(:), allocatable :: X                 
            !! distance variable X bin center used for viscosity calculations
        real(DP), dimension(:), allocatable :: X2                
            !! distance variable X**2 bin center  used for viscosity calculations
        real(DP), dimension(:), allocatable :: deltaA            
            !! differential surface area of ring
        real(DP), dimension(:), allocatable :: mass              
            !! mass of ring particles in bin
        real(DP), dimension(:), allocatable :: sigma 
            !! surface mass density of ring bin
        real(DP), dimension(:), allocatable :: tau               
            !! ring optical depth
        real(DP), dimension(:), allocatable :: nu                
            !! viscocity of the ring bin
        real(DP), dimension(:), allocatable :: Q                 
            !! Toomre parameter of the ring bin
        real(DP), dimension(:), allocatable :: Iz                
            !! polar moment of inertia of ring bin
        real(DP), dimension(:), allocatable :: w                 
            !! Keplerian angular velocity of ring bin
        real(DP), dimension(:), allocatable :: Torque            
            !! total satellite torque density acting on the ring bin
        real(DP), dimension(:), allocatable :: r_p
            !! ring particle radius
        real(DP), dimension(:), allocatable :: m_p
            !! ring particle mass
        real(DP), dimension(:), allocatable :: rho_p
            !! ring particle mass density
        real(DP), dimension(:), allocatable :: vrel_p
            !! ring particle relative velocity
    contains
        procedure :: setup    => ringmoons_util_setup_ring
        procedure :: dealloc  => ringmoons_util_dealloc_ring
        final     ::             ringmoons_final_ring
            !! Finalizes the ringmoons ring object - deallocates all allocatables
    end type ringmoons_ring

    type, extends(base_object) :: ringmoons_seeds
        integer(I4B)                              :: nbody            
            !! Number of satellite seeds
        real(DP)                                  :: feeding_zone_factor 
            !! Width of feeding zone for seed mergers in units of mutual Hill's sphere
        real(DP)                                  :: rkf_tol      
            !! Error tolerance for Runge-Kutta-Fehlberg integrator for seed evolution
        real(DP)                                  :: mass_init       
            !! initial mass of seeds
        logical, dimension(:), allocatable   :: active       
            !! Flag to determine whether this body is active or not
        real(DP), dimension(:), allocatable       :: a            
            !! Semimajor axis of seed
        real(DP), dimension(:), allocatable       :: mass           
            !! Mass of seed
        real(DP), dimension(:), allocatable       :: Rhill        
            !! Hill's sphere radius of seed
        integer(I4B), dimension(:), allocatable   :: rbin         
            !! Ring bin location of seed
        real(DP), dimension(:), allocatable       :: Torque       
            !! Total torque acting on the seed
        real(DP), dimension(:), allocatable       :: Ttide        
            !! Tidal torque acting on the seed
    contains
        procedure :: setup    => ringmoons_util_setup_seeds
        procedure :: dealloc  => ringmoons_util_dealloc_seeds
        final     ::             ringmoons_final_seeds
            !! Finalizes the ringmoons seeds object - deallocates all allocatables
    end type ringmoons_seeds

    type, extends(symba_nbody_system) :: ringmoons_nbody_system
        class(ringmoons_ring),         allocatable :: ring
            !! Ringmoons ring object
        class(ringmoons_seeds),        allocatable :: seeds
            !! Ringmoons seeds object
    end type ringmoons_nbody_system

    interface
        module subroutine ringmoons_util_dealloc_ring(self)
            !! author: David A. Minton
            !!
            !! Deallocates all allocatabale arrays
            implicit none
            ! Arguments
            class(ringmoons_ring),  intent(inout) :: self 
                !! Ringmoons ring object
        end subroutine ringmoons_util_dealloc_ring

        module subroutine ringmoons_util_dealloc_seeds(self)
            !! author: David A. Minton
            !!
            !! Deallocates all allocatabale arrays
            implicit none
            ! Arguments
            class(ringmoons_seeds),  intent(inout) :: self 
                !! Ringmoons ring object
            end subroutine ringmoons_util_dealloc_seeds

        module subroutine ringmoons_util_setup_ring(self, n, param)
            implicit none
            class(ringmoons_ring),      intent(inout) :: self  
                !! Ringmoons ring object
            integer(I4B),               intent(in)    :: n     
                !! Number of bins to allocate space for
            class(swiftest_parameters), intent(in)    :: param 
                !! Current run configuration parameters
        end subroutine ringmoons_util_setup_ring

        module subroutine ringmoons_util_setup_seeds(self, n, param)
            implicit none
            class(ringmoons_seeds),     intent(inout) :: self  
                !! Ringmoons seeds object
            integer(I4B),               intent(in)    :: n     
                !! Number of bins to allocate space for
            class(swiftest_parameters), intent(in)    :: param 
                !! Current run configuration parameters
        end subroutine ringmoons_util_setup_seeds
    end interface

    contains
        subroutine ringmoons_final_ring(self)
            !! author: David A. Minton
            !!
            !! Finalize the ringmoons ring object - deallocates all allocatables
            implicit none
            ! Argument
            type(ringmoons_ring),  intent(inout) :: self 
                !! Ringmoons ring object
            call self%dealloc()
            return
        end subroutine ringmoons_final_ring

        subroutine ringmoons_final_seeds(self)
            !! author: David A. Minton
            !!
            !! Finalize the ringmoons seeds object - deallocates all allocatables
            implicit none
            ! Argument
            type(ringmoons_seeds),  intent(inout) :: self
                !! Ringmoons seeds object
            call self%dealloc()
            return
        end subroutine ringmoons_final_seeds

end module ringmoons