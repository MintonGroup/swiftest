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

    integer(I4B), parameter        :: LOOPMAX = 2147483646 
    integer(I4B),public,parameter  :: N_DISK_FACTOR = 10000 

    !> NetCDF dimension and variable names for the ringmoons objects
    type, extends(netcdf_parameters) :: ringmoons_netcdf_parameters
        character(NAMELEN) :: ringbin_dimname = "ringbin"
            !! name of the ring bin dimension
        integer(I4B) :: ringbin_dimid
            !! ID for the ring bin dimension
        integer(I4B) :: ringbin_varid
            !! ID for the ring bin variable
        integer(I4B) :: nbins
            !! Number of elements in the ring bins
        character(NAMELEN) :: r_varname = "r"
            !! name of the radial distance of bin center variable
        integer(I4B) :: r_varid
            !! ID for the radial distance of bin center variable
        character(NAMELEN) :: sigma_varname = "sigma"
            !! name of the surface mass density of ring bin variable
        integer(I4B) :: sigma_varid
            !! ID for the surface mass density of ring bin variable
        character(NAMELEN) :: r_p_varname = "r_p"
            !! name of the ring particle radius per bin variable
        integer(I4B) :: r_p_varid
            !! ID for the ring particle radius per bin variable
        character(NAMELEN) :: m_p_varname = "m_p"
            !! name of the ring particle mass per bin variable
        integer(I4B) :: m_p_varid
            !! ID for the ring particle mass per bin variable
        character(NAMELEN) :: r_outer_varname = "r_outer"
            !! name of the outside radius of ring variable
        integer(I4B) :: r_outer_varid
            !! ID for the outside radius of ring variable
        character(NAMELEN) :: r_inner_varname = "r_inner"
            !! name of the inside radius of ring variable
        integer(I4B) :: r_inner_varid
            !! ID for the inside radius of ring variable
        character(NAMELEN) :: tau_varname = "tau"
            !! name of the ring optical depth variable
        integer(I4B) :: tau_varid
            !! ID for the ring optical depth variable
        character(NAMELEN) :: nu_varname = "nu"
            !! name of the viscosity of the ring bin variable
        integer(I4B) :: nu_varid
            !! ID for the viscosity of the ring bin variable
        character(NAMELEN) :: toomre_varname = "Q_toomre"
            !! name of the Toomre parameter of the ring bin variable
        integer(I4B) :: toomre_varid
            !! ID for the Toomre parameter of the ring bin variable
        character(NAMELEN) :: vrel_p_varname = "vrel_p"
            !! name of the ring particle relative velocity per bin variable
        integer(I4B) :: vrel_p_varid
            !! ID for the ring particle relative velocity per bin variable
    contains
        procedure :: open       => ringmoons_io_netcdf_open
            !! Open a Ringmoons NetCDF file
        procedure :: flush      => ringmoons_io_netcdf_flush
            !! Flushes a NetCDF file by closing it then opening it again
        final     ::               ringmoons_final_netcdf_parameters 
            !! Finalizer will close the NetCDF file
    end type ringmoons_netcdf_parameters 

    !> Ringmoons central body particle class
    type, extends(symba_cb) :: ringmoons_cb
    contains
        procedure :: accrete => ringmoons_util_accrete_cb
    end type ringmoons_cb

    !> Ringmoons massive body class
    type, extends(symba_pl) :: ringmoons_pl
    contains
    end type ringmoons_pl

    type, extends(base_object) :: ringmoons_seed
        integer(I4B)                                            :: nbody = 0       
            !! Number of seed bodies
        integer(I4B),                 dimension(:), allocatable :: id              
            !! Identifier 
        type(swiftest_particle_info), dimension(:), allocatable :: info            
            !! Particle metadata information
        integer(I4B),                 dimension(:), allocatable :: status          
            !! Active or inactive status indicator
        real(DP),                     dimension(:), allocatable :: a               
            !! Semimajor axis 
        real(DP),                     dimension(:), allocatable :: mass    
            !! Body mass (units MU)
        real(DP),                     dimension(:), allocatable :: Gmass   
            !! Mass gravitational term G * mass (units GU * MU)
        real(DP),                     dimension(:), allocatable :: mu
            !! Gravitational cononstant 
        real(DP),                     dimension(:), allocatable :: rhill   
            !! Hill's radius (units DU)
        real(DP),                     dimension(:), allocatable :: radius  
            !! Body radius (units DU)
        real(DP),                     dimension(:), allocatable :: density 
            !! Body mass density - calculated internally (units MU / DU**3)
        integer(I4B),                 dimension(:), allocatable :: ringbin         
            !! Ring bin location of seed
        real(DP),                     dimension(:), allocatable :: Torque       
            !! Total torque acting on the seed
        real(DP),                     dimension(:), allocatable :: Ttide        
            !! Tidal torque acting on the seed
        real(DP)                                                :: feeding_zone_factor = 10._DP
            !! Width of feeding zone for seed mergers in units of mutual Hill's sphere
        real(DP)                                                :: rkf_tol   = 1.0e-8_DP
            !! Error tolerance for Runge-Kutta-Fehlberg integrator for seed evolution
        real(DP)                                                :: mass_init       
            !! initial mass of seeds
        integer(I4B)                                            :: maxid
            !! Current maximum id of the system, set automatically from the nbody_system object
    contains
        procedure :: setup       => ringmoons_util_setup_seed
            !! Sets up a new seed system from an input file
        procedure :: dealloc     => ringmoons_util_dealloc_seed
            !! Deallocates all allocatable arrays
        procedure :: read_frame  => ringmoons_io_read_frame_seed
            !! Read seed data from file
        procedure :: restructure => ringmoons_step_restructure_seed
            !! Restructures the seed system by merging seeds that are within each other's feeding zones and reassigning ring bins 
            !! based on current semimajor axes
        procedure :: step        => ringmoons_step_seed
            !! Advances the evolution of the seeds by one time step, including accretion and spawning events
        procedure :: spawn       => ringmoons_util_spawn_seed
            !! Spawn new seeds from the ring at the FRL
        procedure :: get_tidal_torque => ringmoons_torque_tidal_seed
            !! Calculates the tidal torque on the seed from the central body
        procedure :: write_frame  => ringmoons_io_write_frame_seed
            !! Writes seed data to file
        final     ::                ringmoons_final_seed
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
        real(DP), dimension(:), allocatable :: Gsigma 
            !! G*surface mass density of ring bin
        real(DP), dimension(:), allocatable :: tau               
            !! ring optical depth
        real(DP), dimension(:), allocatable :: nu                
            !! viscocity of the ring bin
        real(DP), dimension(:), allocatable :: Q                 
            !! Toomre parameter of the ring bin
        real(DP), dimension(:), allocatable :: Iz                
            !! polar moment of inertia of ring bin
        real(DP), dimension(:), allocatable :: wkep             
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
        real(DP)                        :: t = -1.0_DP            
            !! Integration current time (set internally from the nbody_system)
        type(ringmoons_netcdf_parameters) :: nc
            !! NetCDF file object associated with this ring stucture
    contains
        procedure :: setup                   => ringmoons_util_setup_ring
            !! Sets up a new ring system from an input file
        procedure :: reset                   => ringmoons_util_reset_ring
            !!  Resets ring torques and recomputes all dimensional quantities, such as ring extent and limits based on the current
            !! surface mass density and central body properties.
        procedure :: update                  => ringmoons_util_update_ring
            !! Updates the ring velocity dispersion, Toomre parameter, and viscosity values
        procedure :: step                    => ringmoons_step_ring
            !! Adnances the evolution of the ring by one time step.
        procedure :: find_bin                => ringmoons_util_find_bin
            !! Returns the bin containing radius r from the input ring.
        procedure :: get_dt                  => ringmoons_util_get_dt_ring
            !! Calculates the maximum stable timestep for the surface mass density evolution that is not larger than dtin.
        procedure :: set_velocity_dispersion => ringmoons_util_velocity_dispersion_ring
            !! Calculates the velocity 
        procedure :: get_lindblad_torque     => ringmoons_torque_lindblad_ring
            !! Calculates the lindblad torques between each ring element and an input body.
        procedure :: dealloc                 => ringmoons_util_dealloc_ring
            !! Deallocates allocatable arrays
        procedure :: read_frame              => ringmoons_io_read_frame_ring
            !! Read in ring data from file
        procedure :: write_frame             => ringmoons_io_write_frame_ring
            !! Writes ring data to file
        final     ::                            ringmoons_final_ring
            !! Finalizes the ringmoons ring object - deallocates all allocatables
    end type ringmoons_ring

    type, extends(symba_nbody_system) :: ringmoons_nbody_system
        class(ringmoons_ring),         allocatable :: ring
            !! Ringmoons ring object
        class(ringmoons_seed),         allocatable :: seed
            !! Ringmoons seed object
    contains
        procedure :: dealloc    => ringmoons_util_dealloc_system          
            !! Deallocates all allocatables
        procedure :: initialize => ringmoons_util_setup_initialize_system 
            !! Performs ringmoons-specific initilization steps
        procedure :: step       => ringmoons_step_system                  
            !! Advance the ringmoons nbody system forward in time by one step
    end type ringmoons_nbody_system

    interface
        module subroutine ringmoons_io_netcdf_flush(self, param)
            implicit none
            class(ringmoons_netcdf_parameters), intent(inout) :: self 
                !! Parameters used to identify a particular NetCDF dataset
            class(swiftest_parameters),         intent(inout) :: param 
                !! Current run configuration parameters 
        end subroutine ringmoons_io_netcdf_flush

        module subroutine ringmoons_io_netcdf_open(self, param, readonly)
            implicit none
            class(ringmoons_netcdf_parameters), intent(inout) :: self     
                !! Parameters used to identify a particular NetCDF dataset
            class(swiftest_parameters),         intent(in)    :: param    
                !! Current run configuration parameters
            logical, optional,                  intent(in)    :: readonly 
                !! Logical flag indicating that this should be open read only
        end subroutine ringmoons_io_netcdf_open

        module subroutine ringmoons_io_read_frame_ring(self, t, param) 
            implicit none
            class(ringmoons_ring), intent(inout) :: self
            real(DP), intent(in)                  :: t  
            class(swiftest_parameters), intent(in) :: param
        end subroutine ringmoons_io_read_frame_ring

        module subroutine ringmoons_io_read_frame_seed(self, t, nc, param) 
            implicit none
            class(ringmoons_seed),             intent(inout) :: self
            real(DP),                          intent(in)    :: t  
            class(swiftest_netcdf_parameters), intent(inout) :: nc
            class(swiftest_parameters),        intent(inout) :: param
        end subroutine ringmoons_io_read_frame_seed

        module subroutine ringmoons_io_write_frame_ring(self, param) 
            implicit none
            class(ringmoons_ring),      intent(inout) :: self
            class(swiftest_parameters), intent(in) :: param
        end subroutine ringmoons_io_write_frame_ring

        module subroutine ringmoons_io_write_frame_seed(self, nc, param) 
            implicit none
            class(ringmoons_seed),             intent(inout) :: self
            class(swiftest_netcdf_parameters), intent(inout) :: nc
            class(swiftest_parameters),        intent(in)    :: param
        end subroutine ringmoons_io_write_frame_seed

        module subroutine ringmoons_step_restructure_seed(self,cb,ring,param)
            implicit none
            class(ringmoons_seed),      intent(inout) :: self
            class(ringmoons_cb),        intent(inout) :: cb
            class(ringmoons_ring),      intent(inout) :: ring
            class(swiftest_parameters), intent(in)    :: param
        end subroutine ringmoons_step_restructure_seed

        module subroutine ringmoons_step_ring(self,cb,dt,param,stepfail)
            implicit none
            class(ringmoons_ring),      intent(inout) :: self
            class(ringmoons_cb),        intent(in)    :: cb
            real(DP),                   intent(in)    :: dt
            class(swiftest_parameters), intent(in)    :: param
            logical,                    intent(out)   :: stepfail
        end subroutine ringmoons_step_ring

        module subroutine ringmoons_step_seed(self, cb, ring, dt, param, stepfail)
            implicit none
            class(ringmoons_seed),      intent(inout) :: self
            class(ringmoons_cb),        intent(inout) :: cb
            class(ringmoons_ring),      intent(inout) :: ring
            real(DP),                   intent(in)    :: dt
            class(swiftest_parameters), intent(in)    :: param
            logical,                    intent(out)   :: stepfail
        end subroutine ringmoons_step_seed

        module subroutine ringmoons_step_system(self, param, t, dt)
            implicit none
            class(ringmoons_nbody_system), intent(inout) :: self   
                !! Ringmoons nbody system object
            class(swiftest_parameters),    intent(inout) :: param  
                !! Current run configuration parameters
            real(DP),                      intent(in)    :: t      
                !! Simulation time
            real(DP),                      intent(in)    :: dt   
        end subroutine ringmoons_step_system

        module function ringmoons_torque_lindblad_ring(self,cb,asat,esat,isat,msat,param) result(Torque)
            implicit none
            class(ringmoons_ring),      intent(inout) :: self
            class(swiftest_cb),         intent(in)    :: cb 
            real(DP),                   intent(in)    :: asat,esat,isat,msat
            class(swiftest_parameters), intent(in)    :: param
            real(DP),dimension(0:self%nbins+1)        :: Torque
        end function ringmoons_torque_lindblad_ring

        module subroutine ringmoons_torque_tidal_seed(self,cb,param)
            implicit none
            class(ringmoons_seed),      intent(inout) :: self
            class(swiftest_cb),         intent(in)    :: cb
            class(swiftest_parameters), intent(in)    :: param
        end subroutine ringmoons_torque_tidal_seed

        module subroutine ringmoons_util_accrete_cb(self,ring,seed,param,dt)
            implicit none
            class(ringmoons_cb),        intent(inout) :: self
                !! Ringmoons central body object
            class(ringmoons_ring),      intent(inout) :: ring
                !! Ringmoons ring object
            class(ringmoons_seed),      intent(inout) :: seed
                !! Ringmoons seed obje3ct
            class(swiftest_parameters), intent(in)    :: param
                !! Current run configuration parameters
            real(DP),intent(in)                       :: dt
                !! Current time step size
        end subroutine ringmoons_util_accrete_cb

        module subroutine ringmoons_util_dealloc_system(self)
            implicit none
            class(ringmoons_nbody_system), intent(inout) :: self
                !! Ringmoons nbody system object to deallocate
        end subroutine ringmoons_util_dealloc_system

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

        pure elemental module function ringmoons_util_find_bin(self,r) result(bin)
            import ringmoons_ring, DP, I4B
            implicit none
            class(ringmoons_ring), intent(in)      :: self
                !! Ringmoons ring object
            real(DP), intent(in)                   :: r
                !! Radial distance at which to search for the bin
            integer(I4B)                           :: bin
                !! The bin containing radial distance r
        end function ringmoons_util_find_bin

        module function ringmoons_util_get_dt_ring(self,dtin) result(dtout)
            implicit none
            class(ringmoons_ring), intent(in)  :: self
                !! Ringmoons ring object
            real(DP), intent(in)               :: dtin
                !! Input time step size, which serves as a maximum for dtout
            real(DP)                           :: dtout
                !! Output time step size, where dtout <= dtin
        end function ringmoons_util_get_dt_ring

        module subroutine ringmoons_util_update_ring(self,cb)
            implicit none
            class(ringmoons_ring), intent(inout) :: self
            class(ringmoons_cb), intent(in) :: cb
        end subroutine ringmoons_util_update_ring

        module subroutine ringmoons_util_reset_ring(self,seed,cb,param)
            implicit none
            class(ringmoons_ring), intent(inout) :: self
            class(ringmoons_seed), intent(inout) :: seed
            class(ringmoons_cb), intent(in) :: cb
            class(swiftest_parameters), intent(in) :: param
        end subroutine ringmoons_util_reset_ring

        module subroutine ringmoons_util_setup_initialize_system(self, system_history, param)
            implicit none
            class(ringmoons_nbody_system), intent(inout) :: self 
                !! SyMBA nbody_system object
            class(swiftest_storage),allocatable, intent(inout) :: system_history 
                !! Stores the system history between output dumps
            class(swiftest_parameters), intent(inout) :: param 
                !! Current run configuration parameters 
        end subroutine ringmoons_util_setup_initialize_system

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

        module subroutine ringmoons_util_spawn_seed(self, cb, ring, a, mass, param)
            implicit none
            class(ringmoons_seed),          intent(inout) :: self
            class(ringmoons_ring),          intent(inout) :: ring
            class(ringmoons_cb),            intent(in)    :: cb
            real(DP),                       intent(in)    :: a
            real(DP),                       intent(in)    :: mass
            class(swiftest_parameters),     intent(in)    :: param
        end subroutine ringmoons_util_spawn_seed

        module subroutine ringmoons_util_velocity_dispersion_ring(self,cb)
            implicit none
            class(ringmoons_ring), intent(inout) :: self
            class(ringmoons_cb),   intent(in)    :: cb
        end subroutine ringmoons_util_velocity_dispersion_ring

        elemental pure module function ringmoons_transition_function(yin) result(kappa)
            implicit none
            real(DP),intent(in) ::yin
            real(DP) :: kappa
        end function ringmoons_transition_function

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

        subroutine ringmoons_final_seed(self)
            !! author: David A. Minton
            !!
            !! Finalize the ringmoons seed object - deallocates all allocatables
            implicit none
            ! Argument
            type(ringmoons_seed),  intent(inout) :: self 
                !! Ringmoons seed object
            call self%dealloc()
            return
        end subroutine ringmoons_final_seed

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

end module ringmoons