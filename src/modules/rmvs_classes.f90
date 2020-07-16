module rmvs_classes
   !! author: David A. Minton
   !!
   !! Definition of classes and methods specific to the Democratic Heliocentric Method
   !! Partially adapted from David E. Kaufmann's Swifter module: module_whm.f90
   use swiftest_globals
   use swiftest_classes
   implicit none

   !********************************************************************************************************************************
   ! rmvs_configuration class definitions and method interfaces
   !*******************************************************************************************************************************
   type, public, extends(whm_configuration) :: rmvs_configuration
      integer(I4B)         :: integrator     = WHM   !! Symbolic name of the nbody integrator  used
   contains
   end type

   !********************************************************************************************************************************
   ! rmvs_cb class definitions and method interfaces
   !*******************************************************************************************************************************
   !> WHM central body particle class
   type, public, extends(whm_cb) :: rmvs_cb
   contains
   end type rmvs_cb

   !********************************************************************************************************************************
   !                                    rmvs_pl class definitions and method interfaces
   !*******************************************************************************************************************************

   !> WHM massive body particle class
   type, public, extends(whm_pl) :: rmvs_pl
   
      !! Note to developers: If you add componenets to this class, be sure to update methods and subroutines that traverse the
      !!    component list, such as rmvs_setup_pl and rmvs_discard_spill_pl
   contains
      procedure, public :: h2j          => rmvs_coord_h2j_pl        !! Convert position and velcoity vectors from heliocentric to Jacobi coordinates 
      procedure, public :: j2h          => rmvs_coord_j2h_pl        !! Convert position and velcoity vectors from Jacobi to helliocentric coordinates 
      procedure, public :: vh2vj        => rmvs_coord_vh2vj_pl      !! Convert velocity vectors from heliocentric to Jacobi coordinates 
      procedure, public :: setup        => rmvs_setup_pl            !! Constructor method - Allocates space for number of particles
      procedure, public :: getacch      => rmvs_getacch_pl          !! Compute heliocentric accelerations of massive bodies
      procedure, public :: gr_getacch   => rmvs_gr_getacch_pl       !! Acceleration term arising from the post-Newtonian correction
      procedure, public :: gr_p4        => rmvs_gr_p4_pl            !! Position kick due to p**4 term in the post-Newtonian correction
      procedure, public :: gr_vh2pv     => rmvs_gr_vh2pv_pl         !! Converts from heliocentric velocity to psudeovelocity for GR calculations
      procedure, public :: gr_pv2vh     => rmvs_gr_pv2vh_pl         !! Converts from psudeovelocity to heliocentric velocity for GR calculations
      procedure, public :: set_mu       => rmvs_setup_set_mu_eta_pl !! Sets the Jacobi mass value for all massive bodies.
      procedure, public :: user_getacch => rmvs_user_getacch_pl     !! User-defined acceleration
      procedure, public :: drift        => rmvs_drift_pl            !! Loop through massive bodies and call Danby drift routine
   end type rmvs_pl

   interface
      !> WHM massive body constructor method
      module subroutine rmvs_setup_pl(self,n)
         implicit none
         class(rmvs_pl), intent(inout)    :: self !! Swiftest test particle object
         integer, intent(in)             :: n    !! Number of test particles to allocate
      end subroutine rmvs_setup_pl

      module subroutine rmvs_setup_set_mu_eta_pl(self, cb)
         implicit none
         class(rmvs_pl),                intent(inout) :: self    !! Swiftest system object
         class(swiftest_cb), intent(inout) :: cb     !! WHM central body particle data structure
      end subroutine rmvs_setup_set_mu_eta_pl

      !> Get heliocentric accelration of massive bodies
      module subroutine rmvs_getacch_pl(self, cb, config, t)
         implicit none
         class(rmvs_pl),                 intent(inout) :: self   !! WHM massive body particle data structure
         class(rmvs_cb),       intent(inout) :: cb     !! WHM central body particle data structure
         class(swiftest_configuration), intent(in)    :: config !! Input collection of user-defined parameter
         real(DP),                      intent(in)    :: t         !! Current time
      end subroutine rmvs_getacch_pl

      module subroutine rmvs_drift_pl(self, cb, config, dt)
         implicit none
         class(rmvs_pl),                 intent(inout) :: self   !! WHM massive body particle data structure
         class(swiftest_cb),  intent(inout) :: cb     !! WHM central body particle data structur
         class(swiftest_configuration), intent(in)    :: config !! Input collection of user-defined parameter
         real(DP),                      intent(in)    :: dt     !! Stepsize
      end subroutine rmvs_drift_pl

      module subroutine rmvs_getacch_int_pl(self, cb)
         implicit none
         class(rmvs_pl),                 intent(inout) :: self   !! WHM massive body particle data structure
         class(rmvs_cb),       intent(inout) :: cb     !! WHM central body particle data structure
      end subroutine rmvs_getacch_int_pl

      module subroutine rmvs_user_getacch_pl(self, cb, config, t)
         class(rmvs_pl),                 intent(inout) :: self   !! WHM massive body particle data structure
         class(rmvs_cb),       intent(inout) :: cb     !! WHM central body particle data structuree
         class(swiftest_configuration), intent(in)    :: config !! Input collection of user-defined parameter
         real(DP),                      intent(in)    :: t         !! Current time
      end subroutine rmvs_user_getacch_pl

      module subroutine rmvs_coord_h2j_pl(self, cb)
         implicit none
         class(rmvs_pl),                 intent(inout) :: self   !! WHM massive body particle data structure
         class(rmvs_cb),       intent(inout) :: cb     !! WHM central body particle data structuree
      end subroutine rmvs_coord_h2j_pl

      module subroutine rmvs_coord_j2h_pl(self, cb)
         implicit none
         class(rmvs_pl),                 intent(inout) :: self   !! WHM massive body particle data structure
         class(rmvs_cb),       intent(inout) :: cb     !! WHM central body particle data structuree
      end subroutine rmvs_coord_j2h_pl

      module subroutine rmvs_coord_vh2vj_pl(self, cb)
         implicit none
         class(rmvs_pl),       intent(inout) :: self   !! WHM massive body particle data structure
         class(rmvs_cb),       intent(inout) :: cb     !! WHM central body particle data structuree
      end subroutine rmvs_coord_vh2vj_pl

      module subroutine rmvs_gr_getacch_pl(self, cb, config)
         implicit none
         class(rmvs_pl),       intent(inout) :: self   !! WHM massive body particle data structure
         class(swiftest_cb),  intent(inout) :: cb     !! WHM central body particle data structuree
         class(swiftest_configuration), intent(in)    :: config !! Input collection of user-defined parameter
      end subroutine rmvs_gr_getacch_pl

      module pure subroutine rmvs_gr_p4_pl(self, config, dt)
         implicit none
         class(rmvs_pl),                 intent(inout) :: self   !! Swiftest particle object
         class(swiftest_configuration), intent(in)    :: config !! Input collection of user-defined configuration parameters 
         real(DP),                      intent(in)    :: dt     !! Step size
      end subroutine rmvs_gr_p4_pl

      module pure subroutine rmvs_gr_vh2pv_pl(self, config)
         implicit none
         class(rmvs_pl),                 intent(inout) :: self   !! Swiftest particle object
         class(swiftest_configuration), intent(in)    :: config !! Input collection of user-defined configuration parameters 
      end subroutine rmvs_gr_vh2pv_pl

      module pure subroutine rmvs_gr_pv2vh_pl(self, config)
         implicit none
         class(rmvs_pl),                 intent(inout) :: self   !! Swiftest particle object
         class(swiftest_configuration), intent(in)    :: config !! Input collection of user-defined configuration parameters 
      end subroutine rmvs_gr_pv2vh_pl
   end interface

   !********************************************************************************************************************************
   !  rmvs_tp class definitions and method interfaces
   !*******************************************************************************************************************************

   !! WHM test particle class
   type, public, extends(swiftest_tp) :: rmvs_tp
      !! Note to developers: If you add componenets to this class, be sure to update methods and subroutines that traverse the
      !!    component list, such as rmvs_setup_tp and rmvs_discard_spill_tp
   contains
      private
      procedure, public :: setup        => rmvs_setup_tp        !! Allocates new components of the whm class and recursively calls parent allocations
      procedure, public :: getacch      => rmvs_getacch_tp      !! Compute heliocentric accelerations of test particles
      procedure, public :: gr_getacch   => rmvs_gr_getacch_tp   !! Acceleration term arising from the post-Newtonian correction
      procedure, public :: gr_p4        => rmvs_gr_p4_tp        !! Position kick due to p**4 term in the post-Newtonian correction
      procedure, public :: gr_vh2pv     => rmvs_gr_vh2pv_tp     !! Converts from heliocentric velocity to psudeovelocity for GR calculations
      procedure, public :: gr_pv2vh     => rmvs_gr_pv2vh_tp     !! Converts from psudeovelocity to heliocentric velocity for GR calculations
      procedure, public :: user_getacch => rmvs_user_getacch_tp !! User-defined acceleration
      procedure, public :: drift        => rmvs_drift_tp        !! Loop through test particles and call Danby drift routine
   end type rmvs_tp

   interface
      !> WHM test particle constructor 
      module subroutine rmvs_setup_tp(self,n)
         implicit none
         class(rmvs_tp),                 intent(inout) :: self   !! WHM test particle data structure
         integer,                       intent(in)    :: n      !! Number of test particles to allocate
      end subroutine rmvs_setup_tp

      module subroutine rmvs_drift_tp(self, cb, config, dt)
         implicit none
         class(rmvs_tp),                 intent(inout) :: self   !! WHM test particle data structure
         class(swiftest_cb),  intent(inout) :: cb     !! WHM central body particle data structuree
         class(swiftest_configuration), intent(in)    :: config !! Input collection of user-defined parameter
         real(DP),                      intent(in)    :: dt     !! Stepsize
      end subroutine rmvs_drift_tp

      !> Get heliocentric accelration of the test particle
      module subroutine rmvs_getacch_tp(self, cb, pl, config, t)
         implicit none
         class(rmvs_tp),                 intent(inout) :: self   !! WHM test particle data structure
         class(rmvs_cb),       intent(inout) :: cb     !! WHM central body particle data structuree 
         class(rmvs_pl),                 intent(inout) :: pl     !! WHM massive body particle data structure. 
         class(swiftest_configuration), intent(in)    :: config !! Input collection of user-defined parameter
         real(DP),                      intent(in)    :: t         !! Current time
      end subroutine rmvs_getacch_tp

      module subroutine rmvs_user_getacch_tp(self, cb, config, t)
         implicit none
         class(rmvs_tp),                 intent(inout) :: self   !! WHM test particle data structure
         class(rmvs_cb),       intent(inout) :: cb     !! WHM central body particle data structuree
         class(swiftest_configuration), intent(in)    :: config    !! Input collection of user-defined parameter
         real(DP),                      intent(in)    :: t         !! Current time
      end subroutine rmvs_user_getacch_tp

      module subroutine rmvs_gr_getacch_tp(self, cb, config)
         implicit none
         class(rmvs_tp),                 intent(inout) :: self   !! WHM massive body particle data structure
         class(swiftest_cb),  intent(inout) :: cb     !! WHM central body particle data structuree
         class(swiftest_configuration), intent(in)    :: config !! Input collection of user-defined parameter
      end subroutine rmvs_gr_getacch_tp

      module pure subroutine rmvs_gr_p4_tp(self, config, dt)
         implicit none
         class(rmvs_tp),                 intent(inout) :: self   !! Swiftest particle object
         class(swiftest_configuration), intent(in)    :: config !! Input collection of user-defined configuration parameters 
         real(DP),                      intent(in)    :: dt     !! Step size
      end subroutine rmvs_gr_p4_tp

      module pure subroutine rmvs_gr_vh2pv_tp(self, config)
         implicit none
         class(rmvs_tp),                 intent(inout)    :: self   !! Swiftest particle object
         class(swiftest_configuration), intent(in)    :: config !! Input collection of user-defined configuration parameters 
      end subroutine rmvs_gr_vh2pv_tp

      module pure subroutine rmvs_gr_pv2vh_tp(self, config)
         implicit none
         class(rmvs_tp),                 intent(inout)    :: self   !! Swiftest particle object
         class(swiftest_configuration), intent(in)    :: config !! Input collection of user-defined configuration parameters 
      end subroutine rmvs_gr_pv2vh_tp
   end interface

   !********************************************************************************************************************************
   !  rmvs_nbody_system class definitions and method interfaces
   !********************************************************************************************************************************
   !> An abstract class for the WHM integrator nbody system 
   type, public, extends(swiftest_nbody_system) :: rmvs_nbody_system
      !> In the WHM integrator, only test particles are discarded
      class(swiftest_tp), allocatable :: tp_discards
   contains
      private
      !> Replace the abstract procedures with concrete ones
      procedure, public :: initialize    => rmvs_setup_system  !! Performs WHM-specific initilization steps, like calculating the Jacobi masses
      
   end type rmvs_nbody_system

   interface
      module subroutine rmvs_setup_system(self, config)
         implicit none
         class(rmvs_nbody_system),       intent(inout) :: self    !! Swiftest system object
         class(swiftest_configuration), intent(inout) :: config  !! Input collection of user-defined configuration parameters 
      end subroutine rmvs_setup_system

   end interface

   !> Interfaces for all non-type bound whm methods that are implemented in separate submodules 
   interface
      !> Move spilled (discarded) Swiftest basic body components from active list to discard list
      module subroutine rmvs_discard_spill(keeps, discards, lspill_list)
         implicit none
         class(swiftest_body),  intent(inout) :: keeps       !! WHM test particle object
         class(swiftest_body),  intent(inout) :: discards    !! Discarded object 
         logical, dimension(:), intent(in)    :: lspill_list !! Logical array of bodies to spill into the discards
      end subroutine rmvs_discard_spill

      !> Steps the Swiftest nbody system forward in time one stepsize
      module subroutine rmvs_step_system(cb, pl, tp, config)
         implicit none
         class(rmvs_cb),            intent(inout) :: cb      !! WHM central body object  
         class(rmvs_pl),            intent(inout) :: pl      !! WHM central body object  
         class(rmvs_tp),            intent(inout) :: tp      !! WHM central body object  
         class(rmvs_configuration), intent(in)    :: config  !! Input collection of user-defined configuration parameters 
      end subroutine rmvs_step_system
   end interface

end module rmvs_classes
