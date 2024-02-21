! Copyight 2022 - David Minton, Carlisle Wishard, Jennifer Pouplin, Jake Elliott, & Dana Singh
! This file is part of Swiftest.
! Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
! Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
! You should have received a copy of the GNU General Public License along with Swiftest. 
! If not, see: https://www.gnu.org/licenses. 

module encounter
   !! author: The Purdue Swiftest Team - David A. Minton, Carlisle A. Wishard, Jennifer L.L. Pouplin, and Jacob R. Elliott
   !!
   !! Definition of classes and methods used to determine close encounters
   use globals
   use base
   use netcdf_io
   implicit none
   public

   character(len=*), parameter :: ENCOUNTER_OUTFILE = 'encounters.nc'  !! Name of NetCDF output file for encounter information
   real(DP), parameter :: RSWEEP_FACTOR = 1.1_DP

   type, abstract :: encounter_list
      integer(I8B)                              :: nenc = 0   !! Total number of encounters
      real(DP)                                  :: t          !! Time of encounter
      logical                                   :: lcollision !! Indicates if the encounter resulted in at least one collision
      real(DP),     dimension(:),   allocatable :: tcollision !! Time of collision
      logical,      dimension(:),   allocatable :: lclosest   !! indicates that thie pair of bodies is in currently at its closest approach point
      logical,      dimension(:),   allocatable :: lvdotr     !! relative vdotr flag
      integer(I4B), dimension(:),   allocatable :: status     !! status of the interaction
      integer(I4B), dimension(:),   allocatable :: index1     !! position of the first body in the encounter
      integer(I4B), dimension(:),   allocatable :: index2     !! position of the second body in the encounter
      integer(I4B), dimension(:),   allocatable :: id1        !! id of the first body in the encounter
      integer(I4B), dimension(:),   allocatable :: id2        !! id of the second body in the encounter
      real(DP),     dimension(:,:), allocatable :: r1         !! the position of body 1 in the encounter
      real(DP),     dimension(:,:), allocatable :: r2         !! the position of body 2 in the encounter
      real(DP),     dimension(:,:), allocatable :: v1         !! the velocity of body 1 in the encounter
      real(DP),     dimension(:,:), allocatable :: v2         !! the velocity of body 2 in the encounter
      integer(I4B), dimension(:),   allocatable :: level      !! Recursion level (used in SyMBA)
   contains
      procedure :: setup   => encounter_util_setup_list        !! A constructor that sets the number of encounters and allocates and initializes all arrays  
      procedure :: append  => encounter_util_append_list  !! Appends elements from one structure to another
      procedure :: copy    => encounter_util_copy_list    !! Copies elements from the source encounter list into self.
      procedure :: dealloc => encounter_util_dealloc_list !! Deallocates all allocatables
      procedure :: spill   => encounter_util_spill_list   !! "Spills" bodies from one object to another depending on the results of a mask (uses the PACK intrinsic)
      procedure :: resize  => encounter_util_resize_list  !! Checks the current size of the encounter list against the required size and extends it by a factor of 2 more than requested if it is too small.
   end type encounter_list
   

   type, extends(base_object) :: encounter_snapshot
      !! A simplified version of a SyMBA nbody system object for storing minimal snapshots of the system state during encounters
      class(base_object), allocatable :: pl    !! Massive body data structure
      class(base_object), allocatable :: tp    !! Test particle data structure
      real(DP)                        :: t     !! Simulation time when snapshot was taken
      integer(I8B)                    :: iloop !! Loop number at time of snapshot
   contains
      procedure :: write_frame => encounter_io_netcdf_write_frame_snapshot !! Writes a frame of encounter data to file 
      procedure :: dealloc     => encounter_util_dealloc_snapshot          !! Deallocates all allocatables
      procedure :: get_idvals  => encounter_util_get_idvalues_snapshot     !! Gets an array of all id values saved in this snapshot
      final     ::                encounter_final_snapshot
   end type encounter_snapshot

   !> NetCDF dimension and variable names for the enounter save object
   type, extends(netcdf_parameters) :: encounter_netcdf_parameters
   contains
      procedure :: initialize => encounter_io_netcdf_initialize_output !! Initialize a set of parameters used to identify a NetCDF output object
      procedure :: open       => encounter_io_netcdf_open              !! Open an encounter NetCDF file
      final     ::               encounter_final_netcdf_parameters     !! Finalizer will close the NetCDF file
   end type encounter_netcdf_parameters


   !> A class that that is used to store simulation history data between file output
   type, extends(base_storage) :: encounter_storage
      class(encounter_netcdf_parameters), allocatable :: nc             !! NetCDF object attached to this storage object
   contains
      procedure :: dump             => encounter_io_netcdf_dump        !! Dumps contents of encounter history to file
      procedure :: dealloc          => encounter_util_dealloc_storage  !! Deallocates all allocatables
      procedure :: get_index_values => encounter_util_get_vals_storage !! Gets the unique values of the indices of a storage object (i.e. body id or time value)
      procedure :: make_index_map   => encounter_util_index_map        !! Maps body id values to storage index values so we don't have to use unlimited dimensions for id
      procedure :: take_snapshot    => encounter_util_snapshot         !! Take a minimal snapshot of the system through an encounter
      final     ::                     encounter_final_storage
   end type encounter_storage


   type encounter_bounding_box_1D
      integer(I4B)                            :: n    !! Number of bodies with extents
      integer(I4B), dimension(:), allocatable :: ind  !! Sorted minimum/maximum extent indices (value > n indicates an ending index)
      integer(I8B), dimension(:), allocatable :: ibeg !! Beginning index for box
      integer(I8B), dimension(:), allocatable :: iend !! Ending index for box
   contains
      procedure :: sort    => encounter_check_sort_aabb_1D !! Sorts the bounding box extents along a single dimension prior to the sweep phase
      procedure :: dealloc => encounter_util_dealloc_aabb  !! Deallocates all allocatables
      final     ::            encounter_final_aabb    !! Finalize the axis-aligned bounding box (1D) - deallocates all allocatables
   end type


   type encounter_bounding_box
      type(encounter_bounding_box_1D) :: aabb
   contains
      procedure :: dealloc      => encounter_util_dealloc_bounding_box    !! Deallocates all allocatables
      procedure :: setup        => encounter_util_setup_aabb              !! Setup a new axis-aligned bounding box structure
      procedure :: sweep_single => encounter_check_sweep_aabb_single_list !! Sweeps the sorted bounding box extents and returns the encounter candidates
      procedure :: sweep_double => encounter_check_sweep_aabb_double_list !! Sweeps the sorted bounding box extents and returns the encounter candidates
      generic   :: sweep        => sweep_single, sweep_double
      final     ::                 encounter_final_bounding_box
   end type


   interface
      module subroutine encounter_check_all_plpl(param, npl, r, v, renc, dt, nenc, index1, index2, lvdotr)
         use base, only: base_parameters
         implicit none
         class(base_parameters),                  intent(inout) :: param  !! Current Swiftest run configuration parameter5s
         integer(I4B),                            intent(in)    :: npl    !! Total number of massive bodies
         real(DP),     dimension(:,:),            intent(in)    :: r      !! Position vectors of massive bodies
         real(DP),     dimension(:,:),            intent(in)    :: v      !! Velocity vectors of massive bodies
         real(DP),     dimension(:),              intent(in)    :: renc   !! Critical radii of massive bodies that defines an encounter 
         real(DP),                                intent(in)    :: dt     !! Step size
         logical,      dimension(:), allocatable, intent(out)   :: lvdotr !! Logical flag indicating the sign of v .dot. x
         integer(I4B), dimension(:), allocatable, intent(out)   :: index1 !! List of indices for body 1 in each encounter
         integer(I4B), dimension(:), allocatable, intent(out)   :: index2 !! List of indices for body 2 in each encounter
         integer(I8B),                            intent(out)   :: nenc   !! Total number of encounters
      end subroutine encounter_check_all_plpl

      module subroutine encounter_check_all_plplm(param, nplm, nplt, rplm, vplm, rplt, vplt, rencm, renct, dt, &
                                                  nenc, index1, index2, lvdotr)
         use base, only: base_parameters
         implicit none
         class(base_parameters),                  intent(inout) :: param  !! Current Swiftest run configuration parameter5s
         integer(I4B),                            intent(in)    :: nplm   !! Total number of fully interacting massive bodies 
         integer(I4B),                            intent(in)    :: nplt   !! Total number of partially interacting masive bodies (GM < GMTINY) 
         real(DP),     dimension(:,:),            intent(in)    :: rplm   !! Position vectors of fully interacting massive bodies
         real(DP),     dimension(:,:),            intent(in)    :: vplm   !! Velocity vectors of fully interacting massive bodies
         real(DP),     dimension(:,:),            intent(in)    :: rplt   !! Position vectors of partially interacting massive bodies
         real(DP),     dimension(:,:),            intent(in)    :: vplt   !! Velocity vectors of partially interacting massive bodies
         real(DP),     dimension(:),              intent(in)    :: rencm  !! Critical radii of fully interacting massive bodies that defines an encounter
         real(DP),     dimension(:),              intent(in)    :: renct  !! Critical radii of partially interacting massive bodies that defines an encounter
         real(DP),                                intent(in)    :: dt     !! Step size
         integer(I8B),                            intent(out)   :: nenc   !! Total number of encounters
         integer(I4B), dimension(:), allocatable, intent(out)   :: index1 !! List of indices for body 1 in each encounter
         integer(I4B), dimension(:), allocatable, intent(out)   :: index2 !! List of indices for body 2 in each encounter
         logical,      dimension(:), allocatable, intent(out)   :: lvdotr !! Logical flag indicating the sign of v .dot. x
      end subroutine encounter_check_all_plplm

      module subroutine encounter_check_all_pltp(param, npl, ntp, rpl, vpl, rtp, vtp, renc, dt, nenc, index1, index2, lvdotr)
         use base, only: base_parameters
         implicit none
         class(base_parameters),                  intent(inout) :: param  !! Current Swiftest run configuration parameter5s
         integer(I4B),                            intent(in)    :: npl    !! Total number of massive bodies 
         integer(I4B),                            intent(in)    :: ntp    !! Total number of test particles 
         real(DP),     dimension(:,:),            intent(in)    :: rpl    !! Position vectors of massive bodies
         real(DP),     dimension(:,:),            intent(in)    :: vpl    !! Velocity vectors of massive bodies
         real(DP),     dimension(:,:),            intent(in)    :: rtp    !! Position vectors of massive bodies
         real(DP),     dimension(:,:),            intent(in)    :: vtp    !! Velocity vectors of massive bodies
         real(DP),     dimension(:),              intent(in)    :: renc   !! Critical radii of massive bodies that defines an encounter
         real(DP),                                intent(in)    :: dt     !! Step size
         integer(I8B),                            intent(out)   :: nenc   !! Total number of encounters
         integer(I4B), dimension(:), allocatable, intent(out)   :: index1 !! List of indices for body 1 in each encounter
         integer(I4B), dimension(:), allocatable, intent(out)   :: index2 !! List of indices for body 2 in each encounter
         logical,      dimension(:), allocatable, intent(out)   :: lvdotr !! Logical flag indicating the sign of v .dot. x
      end subroutine encounter_check_all_pltp

      elemental module subroutine encounter_check_one(xr, yr, zr, vxr, vyr, vzr, renc, dt, lencounter, lvdotr)
         !$omp declare simd(encounter_check_one)
         implicit none
         real(DP), intent(in)  :: xr, yr, zr    !! Relative distance vector components
         real(DP), intent(in)  :: vxr, vyr, vzr !! Relative velocity vector components
         real(DP), intent(in)  :: renc          !! Critical encounter distance
         real(DP), intent(in)  :: dt            !! Step size
         logical,  intent(out) :: lencounter    !! Flag indicating that an encounter has occurred
         logical,  intent(out) :: lvdotr        !! Logical flag indicating the direction of the v .dot. r vector
      end subroutine encounter_check_one

      module subroutine encounter_check_collapse_ragged_list(ragged_list, n1, nenc, index1, index2, lvdotr)
         implicit none
         class(encounter_list), dimension(:),             intent(in)            :: ragged_list !! The ragged encounter list
         integer(I4B),                                    intent(in)            :: n1          !! Number of bodies 1
         integer(I8B),                                    intent(out)           :: nenc        !! Total number of encountersj 
         integer(I4B),         dimension(:), allocatable, intent(out)           :: index1      !! Array of indices for body 1
         integer(I4B),         dimension(:), allocatable, intent(out)           :: index2      !! Array of indices for body 1
         logical,              dimension(:), allocatable, intent(out), optional :: lvdotr      !! Array indicating which bodies are approaching
      end subroutine encounter_check_collapse_ragged_list

      pure module subroutine encounter_check_sort_aabb_1D(self, n, extent_arr)
         implicit none
         class(encounter_bounding_box_1D), intent(inout) :: self       !! Bounding box structure along a single dimension
         integer(I4B),                     intent(in)    :: n          !! Number of bodies with extents
         real(DP), dimension(:),           intent(in)    :: extent_arr !! Array of extents of size 2*n
      end subroutine encounter_check_sort_aabb_1D

      module subroutine encounter_check_sweep_aabb_double_list(self, n1, n2, r1, v1, r2, v2, renc1, renc2, dt, &
                                                               nenc, index1, index2, lvdotr)
         implicit none
         class(encounter_bounding_box),           intent(inout) :: self       !! Multi-dimensional bounding box structure
         integer(I4B),                            intent(in)    :: n1         !! Number of bodies 1
         integer(I4B),                            intent(in)    :: n2         !! Number of bodies 2
         real(DP),     dimension(:,:),            intent(in)    :: r1, v1     !! Array of indices of bodies 1
         real(DP),     dimension(:,:),            intent(in)    :: r2, v2     !! Array of indices of bodies 2
         real(DP),     dimension(:),              intent(in)    :: renc1      !! Radius of encounter regions of bodies 1
         real(DP),     dimension(:),              intent(in)    :: renc2      !! Radius of encounter regions of bodies 2
         real(DP),                                intent(in)    :: dt         !! Step size
         integer(I8B),                            intent(out)   :: nenc       !! Total number of encounter candidates
         integer(I4B), dimension(:), allocatable, intent(out)   :: index1     !! List of indices for body 1 in each encounter candidate pair
         integer(I4B), dimension(:), allocatable, intent(out)   :: index2     !! List of indices for body 2 in each encounter candidate pair
         logical,      dimension(:), allocatable, intent(out)   :: lvdotr     !! Logical array indicating which pairs are approaching
      end subroutine encounter_check_sweep_aabb_double_list

      module subroutine encounter_check_sweep_aabb_single_list(self, n, r, v, renc, dt, nenc, index1, index2, lvdotr)
         implicit none
         class(encounter_bounding_box),           intent(inout) :: self       !! Multi-dimensional bounding box structure
         integer(I4B),                            intent(in)    :: n          !! Number of bodies
         real(DP),     dimension(:,:),            intent(in)    :: r, v       !! Array of position and velocity vectors 
         real(DP),     dimension(:),              intent(in)    :: renc       !! Radius of encounter regions of bodies 1
         real(DP),                                intent(in)    :: dt         !! Step size
         integer(I8B),                            intent(out)   :: nenc       !! Total number of encounter candidates
         integer(I4B), dimension(:), allocatable, intent(out)   :: index1     !! List of indices for one body in each encounter candidate pair
         integer(I4B), dimension(:), allocatable, intent(out)   :: index2     !! List of indices for the other body in each encounter candidate pair
         logical,      dimension(:), allocatable, intent(out)   :: lvdotr     !! Logical array indicating which pairs are approaching
      end subroutine encounter_check_sweep_aabb_single_list

      module subroutine encounter_io_netcdf_dump(self, param)
         implicit none
         class(encounter_storage), intent(inout)        :: self   !! Encounter storage object
         class(base_parameters),   intent(inout)        :: param  !! Current run configuration parameters 
      end subroutine encounter_io_netcdf_dump

      module subroutine encounter_io_netcdf_initialize_output(self, param)
         implicit none
         class(encounter_netcdf_parameters), intent(inout) :: self    !! Parameters used to identify a particular NetCDF dataset
         class(base_parameters),             intent(in)    :: param   
      end subroutine encounter_io_netcdf_initialize_output

      module subroutine encounter_io_netcdf_open(self, param, readonly)
         implicit none
         class(encounter_netcdf_parameters), intent(inout) :: self     !! Parameters used to identify a particular NetCDF dataset
         class(base_parameters),             intent(in)    :: param    !! Current run configuration parameters
         logical, optional,                  intent(in)    :: readonly !! Logical flag indicating that this should be open read only
      end subroutine encounter_io_netcdf_open

      module subroutine encounter_io_netcdf_write_frame_snapshot(self, history, param)
         implicit none
         class(encounter_snapshot), intent(in)    :: self    !! Swiftest encounter structure
         class(encounter_storage),  intent(inout) :: history !! Encounter storage object
         class(base_parameters),    intent(inout) :: param   !! Current run configuration parameters
      end subroutine encounter_io_netcdf_write_frame_snapshot

      module subroutine encounter_util_setup_aabb(self, n, n_last)
         implicit none
         class(encounter_bounding_box), intent(inout) :: self   !! Swiftest encounter structure
         integer(I4B),                  intent(in)    :: n      !! Number of objects with bounding box extents
         integer(I4B),                  intent(in)    :: n_last !! Number of objects with bounding box extents the previous time this was called
      end subroutine encounter_util_setup_aabb

      module subroutine encounter_util_setup_list(self, n)
         implicit none
         class(encounter_list), intent(inout) :: self !! Swiftest encounter structure
         integer(I8B),          intent(in)    :: n    !! Number of encounters to allocate space for
      end subroutine encounter_util_setup_list

      module subroutine encounter_util_append_list(self, source, lsource_mask)
         implicit none
         class(encounter_list), intent(inout)        :: self         !! Swiftest encounter list object
         class(encounter_list), intent(in)           :: source       !! Source object to append
         logical, dimension(:), intent(in)           :: lsource_mask !! Logical mask indicating which elements to append to
      end subroutine encounter_util_append_list

      module subroutine encounter_util_copy_list(self, source)
         implicit none
         class(encounter_list), intent(inout) :: self   !! Encounter list 
         class(encounter_list), intent(in)    :: source !! Source object to copy into
      end subroutine encounter_util_copy_list

      module subroutine encounter_util_dealloc_aabb(self)
         implicit none
         class(encounter_bounding_box_1D), intent(inout) :: self !!Bounding box structure along a single dimension
      end subroutine encounter_util_dealloc_aabb

      module subroutine encounter_util_dealloc_bounding_box(self)
         implicit none
         class(encounter_bounding_box), intent(inout) :: self !! Bounding box structure
      end subroutine encounter_util_dealloc_bounding_box

      module subroutine encounter_util_dealloc_list(self)
         implicit none
         class(encounter_list), intent(inout) :: self !! Swiftest encounter list object
      end subroutine encounter_util_dealloc_list

      module subroutine encounter_util_dealloc_snapshot(self)
         implicit none
         class(encounter_snapshot), intent(inout) :: self !! Encounter shapshot object
      end subroutine encounter_util_dealloc_snapshot

      module subroutine encounter_util_dealloc_storage(self)
         implicit none
         class(encounter_storage), intent(inout) :: self !! Swiftest storage object
      end subroutine encounter_util_dealloc_storage

      module subroutine encounter_util_get_idvalues_snapshot(self, idvals)
         implicit none
         class(encounter_snapshot),               intent(in)  :: self   !! Encounter snapshot object
         integer(I4B), dimension(:), allocatable, intent(out) :: idvals !! Array of all id values saved in this snapshot
      end subroutine encounter_util_get_idvalues_snapshot

      module subroutine encounter_util_get_vals_storage(self, idvals, tvals)
         class(encounter_storage), intent(in)               :: self   !! Encounter storages object
         integer(I4B), dimension(:),  allocatable, intent(out) :: idvals !! Array of all id values in all snapshots
         real(DP),     dimension(:),  allocatable, intent(out) :: tvals  !! Array of all time values in all snapshots
      end subroutine encounter_util_get_vals_storage 

      module subroutine encounter_util_index_map(self)
         implicit none
         class(encounter_storage), intent(inout) :: self  !! Encounter storage object
      end subroutine encounter_util_index_map

      module subroutine encounter_util_resize_list(self, nnew)
         implicit none
         class(encounter_list), intent(inout) :: self !! Swiftest encounter list 
         integer(I8B),          intent(in)    :: nnew !! New size of list needed
      end subroutine encounter_util_resize_list

      module subroutine encounter_util_snapshot(self, param, nbody_system, t, arg)
         implicit none
         class(encounter_storage),  intent(inout)        :: self         !! Swiftest storage object
         class(base_parameters),    intent(inout)        :: param        !! Current run configuration parameters
         class(base_nbody_system),  intent(inout)        :: nbody_system !! Swiftest nbody system object to store
         real(DP),                  intent(in), optional :: t            !! Time of snapshot if different from system time
         character(*),              intent(in), optional :: arg          !! Optional argument (needed for extended storage type used in collision snapshots)
      end subroutine encounter_util_snapshot

      module subroutine encounter_util_spill_list(self, discards, lspill_list, ldestructive)
         implicit none
         class(encounter_list), intent(inout) :: self         !! Swiftest encounter list 
         class(encounter_list), intent(inout) :: discards     !! Discarded object 
         logical, dimension(:), intent(in)    :: lspill_list  !! Logical array of bodies to spill into the discards
         logical,               intent(in)    :: ldestructive !! Logical flag indicating whether or not this operation should alter body by removing the discard list
      end subroutine encounter_util_spill_list

   end interface

   contains

      subroutine encounter_final_aabb(self)
         !! author: David A. Minton
         !!
         !! Finalize the axis aligned bounding box (1D) - deallocates all allocatables
         implicit none
         ! Arguments
         type(encounter_bounding_box_1D), intent(inout) :: self

         call self%dealloc()

         return
      end subroutine encounter_final_aabb

      subroutine encounter_final_bounding_box(self)
         !! author: David A. Minton
         !!
         !! Finalize the bounding box object
         implicit none
         ! Arguments
         type(encounter_bounding_box), intent(inout) :: self

         call self%dealloc()

         return
      end subroutine encounter_final_bounding_box


      subroutine encounter_final_netcdf_parameters(self)
         !! author: David A. Minton
         !!
         !! Finalize the NetCDF by closing the file
         implicit none
         ! Arguments
         type(encounter_netcdf_parameters), intent(inout) :: self

         call self%close()

         return
      end subroutine encounter_final_netcdf_parameters


      subroutine encounter_final_snapshot(self)
         !! author: David A. Minton
         !!
         !! Deallocates allocatable arrays in an encounter snapshot
         implicit none
         ! Arguments
         type(encounter_snapshot),  intent(inout) :: self !! Encounter storage object

         if (allocated(self%pl)) deallocate(self%pl)
         if (allocated(self%tp)) deallocate(self%tp)
         self%t = 0.0_DP

         return
      end subroutine encounter_final_snapshot


      subroutine encounter_final_storage(self)
         !! author: David A. Minton
         !!
         !! Deallocates allocatable arrays in an encounter snapshot
         implicit none
         ! Arguments
         type(encounter_storage),  intent(inout) :: self !! Encounter storage object

         call self%dealloc()

         return
      end subroutine encounter_final_storage

end module encounter

