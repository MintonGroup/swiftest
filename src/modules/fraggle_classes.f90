module fraggle_classes
   !! author: The Purdue Swiftest Team - David A. Minton, Carlisle A. Wishard, Jennifer L.L. Pouplin, and Jacob R. Elliott
   !!
   !! Definition of classes and methods specific to Fraggel: The Fragment Generation Model
   use swiftest_globals
   use swiftest_classes, only : swiftest_parameters, swiftest_nbody_system, swiftest_cb, swiftest_pl
   implicit none
   public

   integer(I4B),     parameter :: FRAGGLE_NMASS_DIST = 3             !! Number of mass bins returned by the regime calculation (largest fragment, second largest, and remainder)  
   character(len=*), parameter :: FRAGGLE_LOG_OUT    = "fraggle.log" !! Name of log file for Fraggle diagnostic information
   integer(I4B),     parameter :: FRAGGLE_LOG_UNIT   = 76            !! Unit number for Fraggle log file

   !********************************************************************************************************************************
   !                                    fraggle_colliders class definitions and method interfaces
   !*******************************************************************************************************************************
   !> Class definition for the variables that describe the bodies involved in the collision
   type :: fraggle_colliders
      integer(I4B)                                 :: ncoll   !! Number of bodies involved in the collision
      integer(I4B), dimension(:),      allocatable :: idx     !! Index of bodies involved in the collision
      real(DP),     dimension(NDIM,2)              :: xb      !! Two-body equivalent position vectors of the collider bodies prior to collision
      real(DP),     dimension(NDIM,2)              :: vb      !! Two-body equivalent velocity vectors of the collider bodies prior to collision
      real(DP),     dimension(NDIM,2)              :: rot     !! Two-body equivalent principal axes moments of inertia the collider bodies prior to collision
      real(DP),     dimension(NDIM,2)              :: L_spin  !! Two-body equivalent spin angular momentum vectors of the collider bodies prior to collision
      real(DP),     dimension(NDIM,2)              :: L_orbit !! Two-body equivalent orbital angular momentum vectors of the collider bodies prior to collision
      real(DP),     dimension(NDIM,2)              :: Ip      !! Two-body equivalent principal axes moments of inertia the collider bodies prior to collision
      real(DP),     dimension(2)                   :: mass    !! Two-body equivalent mass of the collider bodies prior to the collision
      real(DP),     dimension(2)                   :: radius  !! Two-body equivalent radii of the collider bodies prior to the collision
   contains
      procedure :: regime => fraggle_regime_colliders !! Determine which fragmentation regime the set of colliders will be
   end type fraggle_colliders

   !********************************************************************************************************************************
   !                                    fraggle_fragments class definitions and method interfaces
   !*******************************************************************************************************************************
   !> Class definition for the variables that describe a collection of fragments by Fraggle barycentric coordinates
   type, extends(swiftest_pl) :: fraggle_fragments
      real(DP)                                :: mtot      !! Total mass of fragments       
      real(DP)                                :: Qloss     !! Energy lost during the collision
      real(DP), dimension(FRAGGLE_NMASS_DIST) :: mass_dist !! Distribution of fragment mass determined by the regime calculation (largest fragment, second largest, and remainder)    
      integer(I4B)                            :: regime    !! Collresolve regime code for this collision

      ! Values in a coordinate frame centered on the collider barycenter and collisional system unit vectors (these are used internally by the fragment generation subroutine)
      real(DP), dimension(NDIM)              :: xbcom       !! Center of mass position vector of the collider system in system barycentric coordinates
      real(DP), dimension(NDIM)              :: vbcom       !! Velocity vector of the center of mass of the collider system in system barycentric coordinates
      real(DP), dimension(NDIM)              :: x_coll_unit !! x-direction unit vector of collisional system
      real(DP), dimension(NDIM)              :: y_coll_unit !! y-direction unit vector of collisional system
      real(DP), dimension(NDIM)              :: z_coll_unit !! z-direction unit vector of collisional system
      real(DP), dimension(:,:), allocatable  :: x_coll      !! Array of fragment position vectors in the collisional coordinate frame
      real(DP), dimension(:,:), allocatable  :: v_coll      !! Array of fragment velocity vectors in the collisional coordinate frame
      real(DP), dimension(:,:), allocatable  :: v_r_unit    !! Array of radial direction unit vectors of individual fragments in the collisional coordinate frame
      real(DP), dimension(:,:), allocatable  :: v_t_unit    !! Array of tangential direction unit vectors of individual fragments in the collisional coordinate frame
      real(DP), dimension(:,:), allocatable  :: v_n_unit    !! Array of normal direction unit vectors of individual fragments in the collisional coordinate frame
      real(DP), dimension(:),   allocatable  :: rmag        !! Array of radial distance magnitudes of individual fragments in the collisional coordinate frame 
      real(DP), dimension(:),   allocatable  :: rotmag      !! Array of rotation magnitudes of individual fragments 
      real(DP), dimension(:),   allocatable  :: v_r_mag     !! Array of radial direction velocity magnitudes of individual fragments 
      real(DP), dimension(:),   allocatable  :: v_t_mag     !! Array of tangential direction velocity magnitudes of individual fragments

      ! Energy and momentum book-keeping variables that characterize the whole system of fragments
      real(DP)                  :: ke_orbit  !! Current orbital kinetic energy of the system of fragments in the collisional frame
      real(DP)                  :: ke_spin   !! Current spin kinetic energy of the system of fragments in the collisional frame
      real(DP), dimension(NDIM) :: L_orbit   !! Current orbital angular momentum of the system of fragments in the collisional frame  
      real(DP), dimension(NDIM) :: L_spin    !! Current spin angular momentum of the system of fragments in the collisional frame
      real(DP)                  :: ke_budget !! Total kinetic energy budget for the system of fragmens in the collisional frame
      real(DP), dimension(NDIM) :: L_budget  !! Total angular momentum budget for the system of fragmens in the collisional frame

      ! For the following variables, "before" refers to the *entire* n-body system in its pre-collisional state and "after" refers to the system in its post-collisional state
      real(DP), dimension(NDIM) :: Lorbit_before, Lorbit_after     !! Before/after orbital angular momentum 
      real(DP), dimension(NDIM) :: Lspin_before, Lspin_after       !! Before/after spin angular momentum 
      real(DP), dimension(NDIM) :: Ltot_before, Ltot_after         !! Before/after total system angular momentum 
      real(DP)                  :: ke_orbit_before, ke_orbit_after !! Before/after orbital kinetic energy
      real(DP)                  :: ke_spin_before, ke_spin_after   !! Before/after spin kinetic energy
      real(DP)                  :: pe_before, pe_after             !! Before/after potential energy
      real(DP)                  :: Etot_before, Etot_after         !! Before/after total system energy

      ! Scale factors used to scale dimensioned quantities to a more "natural" system where important quantities (like kinetic energy, momentum) are of order ~1
      real(DP) :: dscale !! Distance dimension scale factor
      real(DP) :: mscale !! Mass scale factor
      real(DP) :: tscale !! Time scale factor
      real(DP) :: vscale !! Velocity scale factor (a convenience unit that is derived from dscale and tscale)
      real(DP) :: Escale !! Energy scale factor (a convenience unit that is derived from dscale, tscale, and mscale)
      real(DP) :: Lscale !! Angular momentum scale factor (a convenience unit that is derived from dscale, tscale, and mscale)
   contains
      procedure :: generate_fragments      => fraggle_generate_fragments         !! Generates a system of fragments in barycentric coordinates that conserves energy and momentum.
      procedure :: accel                   => fraggle_placeholder_accel          !! Placeholder subroutine to fulfill requirement for an accel method
      procedure :: kick                    => fraggle_placeholder_kick           !! Placeholder subroutine to fulfill requirement for a kick method
      procedure :: step                    => fraggle_placeholder_step           !! Placeholder subroutine to fulfill requirement for a step method
      procedure :: set_budgets             => fraggle_set_budgets_fragments      !! Sets the energy and momentum budgets of the fragments based on the collider value
      procedure :: set_coordinate_system   => fraggle_set_coordinate_system      !! Defines the collisional coordinate system, including the unit vectors of both the system and individual fragments. 
      procedure :: set_mass_dist           => fraggle_set_mass_dist_fragments    !! Sets the distribution of mass among the fragments depending on the regime type
      procedure :: set_natural_scale       => fraggle_set_natural_scale_factors  !! Scales dimenional quantities to ~O(1) with respect to the collisional system.  
      procedure :: set_original_scale      => fraggle_set_original_scale_factors !! Restores dimenional quantities back to the original system units
      procedure :: setup                   => fraggle_setup_fragments            !! Allocates arrays for n fragments in a Fraggle system. Passing n = 0 deallocates all arrays.
      procedure :: reset                   => fraggle_setup_reset_fragments      !! Resets all position and velocity-dependent fragment quantities in order to do a fresh calculation (does not reset mass, radius, or other values that get set prior to the call to fraggle_generate)
      procedure :: get_ang_mtm             => fraggle_util_ang_mtm               !! Calcualtes the current angular momentum of the fragments
      procedure :: get_energy_and_momentum => fraggle_util_get_energy_momentum   !! Calculates total system energy in either the pre-collision outcome state (lbefore = .true.) or the post-collision outcome state (lbefore = .false.)
      procedure :: restructure             => fraggle_util_restructure           !! Restructure the inputs after a failed attempt failed to find a set of positions and velocities that satisfy the energy and momentum constraints
   end type fraggle_fragments

   interface
      module subroutine fraggle_generate_fragments(self, colliders, system, param, lfailure)
         use swiftest_classes, only : swiftest_nbody_system, swiftest_parameters
         implicit none
         class(fraggle_fragments),     intent(inout) :: self      !! Fraggle fragment system object 
         class(fraggle_colliders),     intent(inout) :: colliders !! Fraggle colliders object containing the two-body equivalent values of the colliding bodies 
         class(swiftest_nbody_system), intent(inout) :: system    !! Swiftest nbody system object
         class(swiftest_parameters),   intent(inout) :: param     !! Current run configuration parameters 
         logical,                      intent(out)   :: lfailure  !! Answers the question: Should this have been a merger instead?
      end subroutine fraggle_generate_fragments

      module subroutine fraggle_io_log_generate(frag)
         implicit none
         class(fraggle_fragments),   intent(in) :: frag
      end subroutine fraggle_io_log_generate

      module subroutine fraggle_io_log_pl(pl, param)
         implicit none
         class(swiftest_pl),         intent(in) :: pl    !! Swiftest massive body object (only the new bodies generated in a collision)
         class(swiftest_parameters), intent(in) :: param !! Current swiftest run configuration parameters
      end subroutine fraggle_io_log_pl

      module subroutine fraggle_io_log_regime(colliders, frag)
         implicit none
         class(fraggle_colliders),   intent(in) :: colliders
         class(fraggle_fragments),   intent(in) :: frag
      end subroutine fraggle_io_log_regime

      !> The following interfaces are placeholders intended to satisfy the required abstract methods given by the parent class
      module subroutine fraggle_placeholder_accel(self, system, param, t, lbeg)
         use swiftest_classes, only : swiftest_nbody_system, swiftest_parameters
         implicit none
         class(fraggle_fragments),     intent(inout) :: self      !! Fraggle fragment system object 
         class(swiftest_nbody_system), intent(inout) :: system !! Swiftest nbody system object
         class(swiftest_parameters),   intent(inout) :: param  !! Current run configuration parameters 
         real(DP),                     intent(in)    :: t      !! Current simulation time
         logical,                      intent(in)    :: lbeg   !! Optional argument that determines whether or not this is the beginning or end of the step
      end subroutine fraggle_placeholder_accel

      module subroutine fraggle_placeholder_kick(self, system, param, t, dt, lbeg)
         use swiftest_classes, only :  swiftest_nbody_system, swiftest_parameters
         implicit none
         class(fraggle_fragments),     intent(inout) :: self   !! Fraggle fragment system object
         class(swiftest_nbody_system), intent(inout) :: system !! Swiftest nbody system objec
         class(swiftest_parameters),   intent(inout) :: param  !! Current run configuration parameters 
         real(DP),                     intent(in)    :: t      !! Current time
         real(DP),                     intent(in)    :: dt     !! Stepsize
         logical,                      intent(in)    :: lbeg   !! Logical flag indicating whether this is the beginning of the half step or not. 
      end subroutine fraggle_placeholder_kick

      module subroutine fraggle_placeholder_step(self, system, param, t, dt)
         use swiftest_classes, only : swiftest_nbody_system, swiftest_parameters
         implicit none
         class(fraggle_fragments),     intent(inout) :: self   !! Helio massive body particle object
         class(swiftest_nbody_system), intent(inout) :: system !! Swiftest nboody system
         class(swiftest_parameters),   intent(inout) :: param  !! Current run configuration parameters 
         real(DP),                     intent(in)    :: t      !! Current simulation time
         real(DP),                     intent(in)    :: dt     !! Stepsiz
      end subroutine fraggle_placeholder_step

      module subroutine fraggle_regime_colliders(self, frag, system, param)
         implicit none 
         class(fraggle_colliders),     intent(inout) :: self   !! Fraggle colliders object
         class(fraggle_fragments),     intent(inout) :: frag   !! Fraggle fragment system object
         class(swiftest_nbody_system), intent(in)    :: system !! Swiftest nbody system object
         class(swiftest_parameters),   intent(in)    :: param  !! Current Swiftest run configuration parameters
      end subroutine fraggle_regime_colliders

      module subroutine fraggle_set_budgets_fragments(self, colliders)
         implicit none
         class(fraggle_fragments), intent(inout) :: self      !! Fraggle fragment system object
         class(fraggle_colliders), intent(inout) :: colliders !! Fraggle collider system object
      end subroutine  fraggle_set_budgets_fragments

      module subroutine fraggle_set_coordinate_system(self, colliders)
         implicit none
         class(fraggle_fragments), intent(inout) :: self      !! Fraggle fragment system object
         class(fraggle_colliders), intent(inout) :: colliders !! Fraggle collider system object
      end subroutine fraggle_set_coordinate_system

      module subroutine fraggle_set_mass_dist_fragments(self, colliders, param)
         implicit none
         class(fraggle_fragments),     intent(inout) :: self      !! Fraggle fragment system object
         class(fraggle_colliders),     intent(inout) :: colliders !! Fraggle collider system object
         class(swiftest_parameters),   intent(in)    :: param     !! Current Swiftest run configuration parameters
      end subroutine fraggle_set_mass_dist_fragments

      module subroutine fraggle_set_natural_scale_factors(self, colliders)
         implicit none
         class(fraggle_fragments), intent(inout) :: self      !! Fraggle fragment system object
         class(fraggle_colliders), intent(inout) :: colliders !! Fraggle collider system object
      end subroutine fraggle_set_natural_scale_factors

      module subroutine fraggle_set_original_scale_factors(self, colliders)
         implicit none
         class(fraggle_fragments), intent(inout) :: self      !! Fraggle fragment system object
         class(fraggle_colliders), intent(inout) :: colliders !! Fraggle collider system object
      end subroutine fraggle_set_original_scale_factors

      module subroutine fraggle_setup_fragments(self, n, param)
         implicit none
         class(fraggle_fragments),   intent(inout) :: self  !! Fraggle fragment system object
         integer(I4B),               intent(in)    :: n     !! Number of fragments
         class(swiftest_parameters), intent(in)    :: param !! Current swiftest run configuration parameters
      end subroutine fraggle_setup_fragments

      module subroutine fraggle_setup_reset_fragments(self)
         implicit none
         class(fraggle_fragments), intent(inout) :: self
      end subroutine fraggle_setup_reset_fragments

      module subroutine fraggle_util_add_fragments_to_system(frag, colliders, system, param)
         use swiftest_classes, only : swiftest_nbody_system, swiftest_parameters
         implicit none
         class(fraggle_fragments),     intent(in)    :: frag      !! Fraggle fragment system object
         class(fraggle_colliders),     intent(in)    :: colliders !! Fraggle collider system object
         class(swiftest_nbody_system), intent(inout) :: system    !! Swiftest nbody system object
         class(swiftest_parameters),   intent(in)    :: param     !! Current swiftest run configuration parameters
      end subroutine fraggle_util_add_fragments_to_system

      module subroutine fraggle_util_ang_mtm(self) 
         implicit none
         class(fraggle_fragments), intent(inout) :: self !! Fraggle fragment system object
      end subroutine fraggle_util_ang_mtm

      module subroutine fraggle_util_construct_temporary_system(frag, system, param, tmpsys, tmpparam)
         use swiftest_classes, only : swiftest_nbody_system, swiftest_parameters
         implicit none
         class(fraggle_fragments),                   intent(in)  :: frag     !! Fraggle fragment system object
         class(swiftest_nbody_system),               intent(in)  :: system   !! Original swiftest nbody system object
         class(swiftest_parameters),                 intent(in)  :: param    !! Current swiftest run configuration parameters
         class(swiftest_nbody_system), allocatable,  intent(out) :: tmpsys   !! Output temporary swiftest nbody system object
         class(swiftest_parameters),   allocatable,  intent(out) :: tmpparam !! Output temporary configuration run parameters
      end subroutine fraggle_util_construct_temporary_system

      module subroutine fraggle_util_get_energy_momentum(self, colliders, system, param, lbefore)
         use swiftest_classes, only : swiftest_nbody_system, swiftest_parameters
         implicit none
         class(fraggle_fragments),     intent(inout) :: self      !! Fraggle fragment system object
         class(fraggle_colliders),     intent(inout) :: colliders !! Fraggle collider system object
         class(swiftest_nbody_system), intent(inout) :: system    !! Swiftest nbody system object
         class(swiftest_parameters),   intent(inout) :: param     !! Current swiftest run configuration parameters
         logical,                      intent(in)    :: lbefore   !! Flag indicating that this the "before" state of the system, with colliders included and fragments excluded or vice versa
      end subroutine fraggle_util_get_energy_momentum

      module subroutine fraggle_util_restructure(self, colliders, try, f_spin, r_max_start)
         implicit none
         class(fraggle_fragments), intent(inout) :: self        !! Fraggle fragment system object
         class(fraggle_colliders), intent(in)    :: colliders   !! Fraggle collider system object
         integer(I4B),             intent(in)    :: try         !! The current number of times Fraggle has tried to find a solution
         real(DP),                 intent(inout) :: f_spin      !! Fraction of energy/momentum that goes into spin. This decreases ater a failed attempt
         real(DP),                 intent(inout) :: r_max_start !! The maximum radial distance that the position calculation starts with. This increases after a failed attempt
      end subroutine fraggle_util_restructure

      module subroutine fraggle_util_shift_vector_to_origin(m_frag, vec_frag)
         implicit none
         real(DP), dimension(:),   intent(in)    :: m_frag    !! Fragment masses
         real(DP), dimension(:,:), intent(inout) :: vec_frag  !! Fragment positions or velocities in the center of mass frame
      end subroutine

      module function fraggle_util_vmag_to_vb(v_r_mag, v_r_unit, v_t_mag, v_t_unit, m_frag, vcom) result(vb) 
         implicit none
         real(DP), dimension(:),   intent(in)  :: v_r_mag   !! Unknown radial component of fragment velocity vector
         real(DP), dimension(:),   intent(in)  :: v_t_mag   !! Tangential component of velocity vector set previously by angular momentum constraint
         real(DP), dimension(:,:), intent(in)  :: v_r_unit, v_t_unit !! Radial and tangential unit vectors for each fragment
         real(DP), dimension(:),   intent(in)  :: m_frag    !! Fragment masses
         real(DP), dimension(:),   intent(in)  :: vcom      !! Barycentric velocity of collisional system center of mass
         real(DP), dimension(:,:), allocatable   :: vb
      end function fraggle_util_vmag_to_vb
   end interface

end module fraggle_classes