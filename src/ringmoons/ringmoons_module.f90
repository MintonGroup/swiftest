! Copyright 2026 - The Minton Group at Purdue University
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
        real(DP) :: mass_init
            !! initial mass of central body
        real(DP) :: mass_accreted
            !! mass accreted by central body through ring updates
    end type ringmoons_cb


    !> Ringmoons massive body class
    type, extends(symba_pl) :: ringmoons_pl
    contains
        procedure :: setup    => ringmoons_util_setup_pl
        procedure :: dealloc  => ringmoons_util_dealloc_pl
    end type ringmoons_pl

    type, extends(ringmoons_pl) :: ringmoons_seed
        real(DP)                                  :: feeding_zone_factor 
            !! Width of feeding zone for seed mergers in units of mutual Hill's sphere
        real(DP)                                  :: rkf_tol      
            !! Error tolerance for Runge-Kutta-Fehlberg integrator for seed evolution
        real(DP)                                  :: mass_init       
            !! initial mass of seeds
        integer(I4B), dimension(:), allocatable   :: ringbin         
            !! Ring bin location of seed
        real(DP), dimension(:), allocatable       :: Torque       
            !! Total torque acting on the seed
        real(DP), dimension(:), allocatable       :: Ttide        
            !! Tidal torque acting on the seed
    contains
        procedure :: setup    => ringmoons_util_setup_seed
        procedure :: dealloc  => ringmoons_util_dealloc_seed
    end type ringmoons_seed


    !> Ringmoons test particle class
    type, extends(symba_tp) :: ringmoons_tp
    contains
    end type ringmoons_tp

    type, extends(base_object) :: ringmoons_ring
        integer(I4B) :: nbins               
            !! number of bins in ring
        integer(I4B) :: inside = 1          
            !! bin id of innermost ring bin (can increase if primary accretes a lot mass through updates)
        real(DP)     :: r_outer             
            !! outside radius of ring in simulation length units
        real(DP)     :: X_outer             
            !! outside radius of ring in X units (see Bath & Pringle 1981)
        real(DP)     :: r_inner             
            !! inside radius of ring in simulation length units
        real(DP)     :: X_inner            
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
        real(DP), dimension(:), allocatable :: vkep             
            !! Keplerian angular velocity of ring bin
        real(DP), dimension(:), allocatable :: Torque            
            !! total satellite torque density acting on the ring bin
        real(DP), dimension(:), allocatable :: r_p
            !! ring particle radius per bin 
        real(DP), dimension(:), allocatable :: m_p
            !! ring particle mass per bin
        real(DP), dimension(:), allocatable :: rho_p
            !! ring particle mass density per bin
        real(DP), dimension(:), allocatable :: vrel_p
            !! ring particle relative velocity per bin
    contains
        procedure :: setup    => ringmoons_util_setup_ring
        procedure :: dealloc  => ringmoons_util_dealloc_ring
        final     ::             ringmoons_final_ring
            !! Finalizes the ringmoons ring object - deallocates all allocatables
    end type ringmoons_ring


    type, extends(symba_nbody_system) :: ringmoons_nbody_system
        class(ringmoons_ring),         allocatable :: ring
            !! Ringmoons ring object
        class(ringmoons_pl),           allocatable :: seeds
            !! Ringmoons seeds object
    contains
        procedure :: step             => ringmoons_step_system      
    end type ringmoons_nbody_system


    !> NetCDF dimension and variable names for the ringmoons objects
    type, extends(netcdf_parameters) :: ringmoons_netcdf_parameters
    contains
        procedure :: initialize => ringmoons_io_netcdf_initialize_output 
            !! Initialize a set of parameters used to identify a NetCDF output object
        final     ::               ringmoons_final_netcdf_parameters 
            !! Finalizer will close the NetCDF file
    end type ringmoons_netcdf_parameters 

    type, extends(base_storage) :: ringmoons_storage
        class(ringmoons_netcdf_parameters), allocatable :: nc             
            !! NetCDF object attached to this storage object
    contains
        procedure :: dump             => ringmoons_io_netcdf_dump        
            !! Dumps contents of ringmoons history to file
        procedure :: dealloc          => ringmoons_util_dealloc_storage  
            !! Deallocates all allocatables
        procedure :: take_snapshot => ringmoons_util_snapshot
            !! Take a snapshot of the ring to save to file
        final     ::                     ringmoons_final_storage
    end type ringmoons_storage


    interface
        module subroutine ringmoons_io_netcdf_dump(self, param)
            implicit none
            class(ringmoons_storage), intent(inout)        :: self   
                !! ringmoons storage object
            class(base_parameters),   intent(inout)        :: param  
                !! Current run configuration parameters 
        end subroutine ringmoons_io_netcdf_dump

        module subroutine ringmoons_io_netcdf_initialize_output(self, param)
            implicit none
            class(ringmoons_netcdf_parameters), intent(inout) :: self    
                !! Parameters used to identify a particular NetCDF dataset
            class(base_parameters),             intent(in)    :: param   
        end subroutine ringmoons_io_netcdf_initialize_output

        module subroutine ringmoons_io_netcdf_open(self, param, readonly)
            implicit none
            class(ringmoons_netcdf_parameters), intent(inout) :: self     
                !! Parameters used to identify a particular NetCDF dataset
            class(base_parameters),             intent(in)    :: param    
                !! Current run configuration parameters
            logical, optional,                  intent(in)    :: readonly 
                !! Logical flag indicating that this should be open read only
        end subroutine ringmoons_io_netcdf_open


        module subroutine ringmoons_step_system(self, param, t, dt)
            implicit none
            class(ringmoons_nbody_system),  intent(inout) :: self   
                !! Ringmoons nbody system object
            class(swiftest_parameters), intent(inout) :: param  
                !! Current run configuration parameters
            real(DP),                   intent(in)    :: t      
                !! Simulation time
            real(DP),                   intent(in)    :: dt   
        end subroutine ringmoons_step_system


        module subroutine ringmoons_util_dealloc_ring(self)
            !! author: David A. Minton
            !!
            !! Deallocates all allocatabale arrays
            implicit none
            ! Arguments
            class(ringmoons_ring),  intent(inout) :: self 
                !! Ringmoons ring object
        end subroutine ringmoons_util_dealloc_ring

        module subroutine ringmoons_util_dealloc_seed(self)
            !! author: David A. Minton
            !!
            !! Deallocates all allocatabale arrays
            implicit none
            ! Arguments
            class(ringmoons_seed),  intent(inout) :: self 
                !! Ringmoons seed object
        end subroutine ringmoons_util_dealloc_seed

        module subroutine ringmoons_util_dealloc_storage(self)
            implicit none
            class(ringmoons_storage), intent(inout) :: self 
                !! Ringmoons storage object
        end subroutine ringmoons_util_dealloc_storage

        module subroutine ringmoons_util_setup_ring(self, n, param)
            implicit none
            class(ringmoons_ring),      intent(inout) :: self  
                !! Ringmoons ring object
            integer(I4B),               intent(in)    :: n     
                !! Number of bins to allocate space for
            class(swiftest_parameters), intent(in)    :: param 
                !! Current run configuration parameters
        end subroutine ringmoons_util_setup_ring

        module subroutine ringmoons_util_setup_pl(self, n, param)
            implicit none
            class(ringmoons_pl),     intent(inout) :: self  
                !! Ringmoons massive body object
            integer(I4B),               intent(in)    :: n     
                !! Number of bins to allocate space for
            class(swiftest_parameters), intent(in)    :: param 
                !! Current run configuration parameters
        end subroutine ringmoons_util_setup_pl

        module subroutine ringmoons_util_setup_seed(self, n, param)
            implicit none
            class(ringmoons_seed),     intent(inout) :: self  
                !! Ringmoons seed object
            integer(I4B),               intent(in)    :: n     
                !! Number of bins to allocate space for
            class(swiftest_parameters), intent(in)    :: param 
                !! Current run configuration parameters
        end subroutine ringmoons_util_setup_seed

        module subroutine ringmoons_util_snapshot(self, param, nbody_system, t, arg)
            implicit none
            class(ringmoons_storage),      intent(inout)        :: self            
                !! Swiftest storage object
            class(swiftest_parameters),   intent(inout)        :: param           
                !! Current run configuration parameters
            class(swiftest_nbody_system), intent(inout)        :: nbody_system    
                !! Swiftest nbody system object to store
            real(DP),                     intent(in), optional :: t               
                !! Time of snapshot if different from nbody_system time
            character(*),                 intent(in), optional :: arg             
                !! Optional argument 
        end subroutine ringmoons_util_snapshot


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

        subroutine ringmoons_final_netcdf_parameters(self)
            !! author: David A. Minton
            !!
            !! Finalize the NetCDF by closing the file
            implicit none
            ! Arguments
            type(ringmoons_netcdf_parameters), intent(inout) :: self

            call self%close()

            return
        end subroutine ringmoons_final_netcdf_parameters

        subroutine ringmoons_final_storage(self)
            !! author: David A. Minton
            !!
            !! Deallocates allocatable arrays in an ringmoons ring snapshot
            implicit none
            ! Arguments
            type(ringmoons_storage),  intent(inout) :: self 
                !! Ringmoons storage object
            call self%dealloc()
            return
        end subroutine ringmoons_final_storage

end module ringmoons