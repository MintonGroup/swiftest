module rmvs_classes
   !! author: David A. Minton
   !!
   !! Definition of classes and methods specific to the Democratic Heliocentric Method
   !! Partially adapted from David E. Kaufmann's Swifter module: module_whm.f90
   use swiftest_globals
   use swiftest_classes
   use whm_classes
   implicit none

   !********************************************************************************************************************************
   ! rmvs_configuration class definitions and method interfaces
   !*******************************************************************************************************************************
   type, public, extends(whm_configuration) :: rmvs_configuration
   contains
   end type

   !********************************************************************************************************************************
   ! rmvs_cb class definitions and method interfaces
   !*******************************************************************************************************************************
   !> RMVS central body particle class
   type, public, extends(whm_cb) :: rmvs_cb
   contains
   end type rmvs_cb

   !********************************************************************************************************************************
   !                                    rmvs_pl class definitions and method interfaces
   !*******************************************************************************************************************************

   !> RMVS massive body particle class
   type, public, extends(whm_pl) :: rmvs_pl
   
      !! Note to developers: If you add componenets to this class, be sure to update methods and subroutines that traverse the
      !!    component list, such as rmvs_setup_pl and rmvs_discard_spill_pl
   contains
   end type rmvs_pl

   interface
   end interface

   !********************************************************************************************************************************
   !  rmvs_tp class definitions and method interfaces
   !*******************************************************************************************************************************

   !! RMVS test particle class
   type, public, extends(whm_tp) :: rmvs_tp
      !! Note to developers: If you add componenets to this class, be sure to update methods and subroutines that traverse the
      !!    component list, such as rmvs_setup_tp and rmvs_discard_spill_tp
   contains
   end type rmvs_tp

   interface
   end interface

   !********************************************************************************************************************************
   !  rmvs_nbody_system class definitions and method interfaces
   !********************************************************************************************************************************
   !> An abstract class for the RMVS integrator nbody system 
   type, public, extends(whm_nbody_system) :: rmvs_nbody_system
      !> In the RMVS integrator, only test particles are discarded
   contains
      private
      
   end type rmvs_nbody_system

   interface

   end interface

   !> Interfaces for all non-type bound whm methods that are implemented in separate submodules 
   interface
      module subroutine rmvs_chk_ind(xr, vr, dt, r2crit, iflag)
         implicit none
         real(DP), intent(in)      :: dt, r2crit
         real(DP), dimension(NDIM), intent(in) :: xr, vr
         integer(I4B), intent(out)      :: iflag
      end subroutine rmvs_chk_ind
   end interface

end module rmvs_classes
