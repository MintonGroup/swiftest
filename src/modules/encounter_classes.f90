module encounter_classes
   !! author: The Purdue Swiftest Team - David A. Minton, Carlisle A. Wishard, Jennifer L.L. Pouplin, and Jacob R. Elliott
   !!
   !! Definition of classes and methods used to determine close encounters
   use swiftest_globals
   use swiftest_classes
   implicit none
   public

   integer(I4B), parameter :: SWEEPDIM = 3

   type :: encounter_list
      integer(I8B)                              :: nenc = 0  !! Total number of encounters
      logical,      dimension(:),   allocatable :: lvdotr !! relative vdotr flag
      integer(I4B), dimension(:),   allocatable :: status !! status of the interaction
      integer(I4B), dimension(:),   allocatable :: index1 !! position of the first body in the encounter
      integer(I4B), dimension(:),   allocatable :: index2 !! position of the second body in the encounter
      integer(I4B), dimension(:),   allocatable :: id1    !! id of the first body in the encounter
      integer(I4B), dimension(:),   allocatable :: id2    !! id of the second body in the encounter
      real(DP),     dimension(:,:), allocatable :: x1     !! the position of body 1 in the encounter
      real(DP),     dimension(:,:), allocatable :: x2     !! the position of body 2 in the encounter
      real(DP),     dimension(:,:), allocatable :: v1     !! the velocity of body 1 in the encounter
      real(DP),     dimension(:,:), allocatable :: v2     !! the velocity of body 2 in the encounter
      real(DP),     dimension(:),   allocatable :: t      !! Time of encounter
   contains
      procedure :: setup   => encounter_setup_list        !! A constructor that sets the number of encounters and allocates and initializes all arrays  
      procedure :: append  => encounter_util_append_list  !! Appends elements from one structure to another
      procedure :: copy    => encounter_util_copy_list    !! Copies elements from the source encounter list into self.
      procedure :: dealloc => encounter_util_dealloc_list !! Deallocates all allocatables
      procedure :: spill   => encounter_util_spill_list   !! "Spills" bodies from one object to another depending on the results of a mask (uses the PACK intrinsic)
      procedure :: resize  => encounter_util_resize_list  !! Checks the current size of the encounter list against the required size and extends it by a factor of 2 more than requested if it is too small.
      procedure :: write   => encounter_io_write_list     !! Write close encounter data to output binary file
      final     :: encounter_util_final_list            !! Finalize the encounter list - deallocates all allocatables
   end type encounter_list

   type encounter_bounding_box_1D
      integer(I4B)                            :: n    !! Number of bodies with extents
      integer(I4B), dimension(:), allocatable :: ind  !! Sorted minimum/maximum extent indices (value > n indicates an ending index)
      integer(I8B), dimension(:), allocatable :: ibeg !! Beginning index for box
      integer(I8B), dimension(:), allocatable :: iend !! Ending index for box
   contains
      procedure :: sort  => encounter_check_sort_aabb_1D !! Sorts the bounding box extents along a single dimension prior to the sweep phase
      procedure :: dealloc => encounter_util_dealloc_aabb !! Deallocates all allocatables
      final     :: encounter_util_final_aabb             !! Finalize the axis-aligned bounding box (1D) - deallocates all allocatables
   end type

   type encounter_bounding_box
      type(encounter_bounding_box_1D), dimension(SWEEPDIM) :: aabb
   contains
      procedure :: setup        => encounter_setup_aabb      !! Setup a new axis-aligned bounding box structure
      procedure :: sweep_single => encounter_check_sweep_aabb_single_list !! Sweeps the sorted bounding box extents and returns the encounter candidates
      procedure :: sweep_double => encounter_check_sweep_aabb_double_list !! Sweeps the sorted bounding box extents and returns the encounter candidates
      generic   :: sweep        => sweep_single, sweep_double
   end type

   interface
      module subroutine encounter_check_all_plpl(param, npl, x, v, renc, dt, nenc, index1, index2, lvdotr)
         use swiftest_classes, only: swiftest_parameters
         implicit none
         class(swiftest_parameters),              intent(inout) :: param  !! Current Swiftest run configuration parameter5s
         integer(I4B),                            intent(in)    :: npl    !! Total number of massive bodies
         real(DP),     dimension(:,:),            intent(in)    :: x      !! Position vectors of massive bodies
         real(DP),     dimension(:,:),            intent(in)    :: v      !! Velocity vectors of massive bodies
         real(DP),     dimension(:),              intent(in)    :: renc   !! Critical radii of massive bodies that defines an encounter 
         real(DP),                                intent(in)    :: dt     !! Step size
         logical,      dimension(:), allocatable, intent(out)   :: lvdotr !! Logical flag indicating the sign of v .dot. x
         integer(I4B), dimension(:), allocatable, intent(out)   :: index1 !! List of indices for body 1 in each encounter
         integer(I4B), dimension(:), allocatable, intent(out)   :: index2 !! List of indices for body 2 in each encounter
         integer(I8B),                            intent(out)   :: nenc   !! Total number of encounters
      end subroutine encounter_check_all_plpl

      module subroutine encounter_check_all_plplm(param, nplm, nplt, xplm, vplm, xplt, vplt, rencm, renct, dt, &
                                                  nenc, index1, index2, lvdotr)
         use swiftest_classes, only: swiftest_parameters
         implicit none
         class(swiftest_parameters),              intent(inout) :: param  !! Current Swiftest run configuration parameter5s
         integer(I4B),                            intent(in)    :: nplm   !! Total number of fully interacting massive bodies 
         integer(I4B),                            intent(in)    :: nplt   !! Total number of partially interacting masive bodies (GM < GMTINY) 
         real(DP),     dimension(:,:),            intent(in)    :: xplm   !! Position vectors of fully interacting massive bodies
         real(DP),     dimension(:,:),            intent(in)    :: vplm   !! Velocity vectors of fully interacting massive bodies
         real(DP),     dimension(:,:),            intent(in)    :: xplt   !! Position vectors of partially interacting massive bodies
         real(DP),     dimension(:,:),            intent(in)    :: vplt   !! Velocity vectors of partially interacting massive bodies
         real(DP),     dimension(:),              intent(in)    :: rencm  !! Critical radii of fully interacting massive bodies that defines an encounter
         real(DP),     dimension(:),              intent(in)    :: renct  !! Critical radii of partially interacting massive bodies that defines an encounter
         real(DP),                                intent(in)    :: dt     !! Step size
         integer(I8B),                            intent(out)   :: nenc   !! Total number of encounters
         integer(I4B), dimension(:), allocatable, intent(out)   :: index1 !! List of indices for body 1 in each encounter
         integer(I4B), dimension(:), allocatable, intent(out)   :: index2 !! List of indices for body 2 in each encounter
         logical,      dimension(:), allocatable, intent(out)   :: lvdotr !! Logical flag indicating the sign of v .dot. x
      end subroutine encounter_check_all_plplm

      module subroutine encounter_check_all_pltp(param, npl, ntp, xpl, vpl, xtp, vtp, renc, dt, nenc, index1, index2, lvdotr)
         use swiftest_classes, only: swiftest_parameters
         implicit none
         class(swiftest_parameters),              intent(inout) :: param  !! Current Swiftest run configuration parameter5s
         integer(I4B),                            intent(in)    :: npl    !! Total number of massive bodies 
         integer(I4B),                            intent(in)    :: ntp    !! Total number of test particles 
         real(DP),     dimension(:,:),            intent(in)    :: xpl    !! Position vectors of massive bodies
         real(DP),     dimension(:,:),            intent(in)    :: vpl    !! Velocity vectors of massive bodies
         real(DP),     dimension(:,:),            intent(in)    :: xtp    !! Position vectors of massive bodies
         real(DP),     dimension(:,:),            intent(in)    :: vtp    !! Velocity vectors of massive bodies
         real(DP),     dimension(:),              intent(in)    :: renc   !! Critical radii of massive bodies that defines an encounter
         real(DP),                                intent(in)    :: dt     !! Step size
         integer(I8B),                            intent(out)   :: nenc   !! Total number of encounters
         integer(I4B), dimension(:), allocatable, intent(out)   :: index1 !! List of indices for body 1 in each encounter
         integer(I4B), dimension(:), allocatable, intent(out)   :: index2 !! List of indices for body 2 in each encounter
         logical,      dimension(:), allocatable, intent(out)   :: lvdotr !! Logical flag indicating the sign of v .dot. x
      end subroutine encounter_check_all_pltp

      module elemental subroutine encounter_check_one(xr, yr, zr, vxr, vyr, vzr, renc, dt, lencounter, lvdotr)
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
         type(encounter_list), dimension(:),              intent(in)            :: ragged_list !! The ragged encounter list
         integer(I4B),                                    intent(in)            :: n1          !! Number of bodies 1
         integer(I8B),                                    intent(out)           :: nenc        !! Total number of encountersj 
         integer(I4B),         dimension(:), allocatable, intent(out)           :: index1      !! Array of indices for body 1
         integer(I4B),         dimension(:), allocatable, intent(out)           :: index2      !! Array of indices for body 1
         logical,              dimension(:), allocatable, intent(out), optional :: lvdotr      !! Array indicating which bodies are approaching
      end subroutine encounter_check_collapse_ragged_list

      module pure subroutine encounter_check_sort_aabb_1D(self, n, extent_arr)
         implicit none
         class(encounter_bounding_box_1D), intent(inout) :: self       !! Bounding box structure along a single dimension
         integer(I4B),                     intent(in)    :: n          !! Number of bodies with extents
         real(DP), dimension(:),           intent(in)    :: extent_arr !! Array of extents of size 2*n
      end subroutine encounter_check_sort_aabb_1D

      module subroutine encounter_check_sweep_aabb_double_list(self, n1, n2, x1, v1, x2, v2, renc1, renc2, dt, &
                                                               nenc, index1, index2, lvdotr)
         implicit none
         class(encounter_bounding_box),           intent(inout) :: self       !! Multi-dimensional bounding box structure
         integer(I4B),                            intent(in)    :: n1         !! Number of bodies 1
         integer(I4B),                            intent(in)    :: n2         !! Number of bodies 2
         real(DP),     dimension(:,:),            intent(in)    :: x1, v1     !! Array of indices of bodies 1
         real(DP),     dimension(:,:),            intent(in)    :: x2, v2     !! Array of indices of bodies 2
         real(DP),     dimension(:),              intent(in)    :: renc1      !! Radius of encounter regions of bodies 1
         real(DP),     dimension(:),              intent(in)    :: renc2      !! Radius of encounter regions of bodies 2
         real(DP),                                intent(in)    :: dt         !! Step size
         integer(I8B),                            intent(out)   :: nenc       !! Total number of encounter candidates
         integer(I4B), dimension(:), allocatable, intent(out)   :: index1     !! List of indices for body 1 in each encounter candidate pair
         integer(I4B), dimension(:), allocatable, intent(out)   :: index2     !! List of indices for body 2 in each encounter candidate pair
         logical,      dimension(:), allocatable, intent(out)   :: lvdotr     !! Logical array indicating which pairs are approaching
      end subroutine encounter_check_sweep_aabb_double_list

      module subroutine encounter_check_sweep_aabb_single_list(self, n, x, v, renc, dt, nenc, index1, index2, lvdotr)
         implicit none
         class(encounter_bounding_box),           intent(inout) :: self       !! Multi-dimensional bounding box structure
         integer(I4B),                            intent(in)    :: n          !! Number of bodies
         real(DP),     dimension(:,:),            intent(in)    :: x, v       !! Array of position and velocity vectors 
         real(DP),     dimension(:),              intent(in)    :: renc       !! Radius of encounter regions of bodies 1
         real(DP),                                intent(in)    :: dt         !! Step size
         integer(I8B),                            intent(out)   :: nenc       !! Total number of encounter candidates
         integer(I4B), dimension(:), allocatable, intent(out)   :: index1     !! List of indices for one body in each encounter candidate pair
         integer(I4B), dimension(:), allocatable, intent(out)   :: index2     !! List of indices for the other body in each encounter candidate pair
         logical,      dimension(:), allocatable, intent(out)   :: lvdotr     !! Logical array indicating which pairs are approaching
      end subroutine encounter_check_sweep_aabb_single_list

      module subroutine encounter_io_write_frame(iu, t, id1, id2, Gmass1, Gmass2, radius1, radius2, xh1, xh2, vh1, vh2)
         implicit none
         integer(I4B),           intent(in) :: iu               !! Open file unit number
         real(DP),               intent(in) :: t                !! Time of encounter
         integer(I4B),           intent(in) :: id1, id2         !! ids of the two encountering bodies
         real(DP),               intent(in) :: Gmass1, Gmass2   !! G*mass of the two encountering bodies
         real(DP),               intent(in) :: radius1, radius2 !! Radii of the two encountering bodies
         real(DP), dimension(:), intent(in) :: xh1, xh2         !! Swiftestcentric position vectors of the two encountering bodies 
         real(DP), dimension(:), intent(in) :: vh1, vh2         !! Swiftestcentric velocity vectors of the two encountering bodies 
      end subroutine encounter_io_write_frame

      module subroutine encounter_io_write_list(self, pl, encbody, param)
         use swiftest_classes, only : swiftest_pl, swiftest_body, swiftest_parameters
         implicit none
         class(encounter_list),      intent(in) :: self    !! Swiftest encounter list object
         class(swiftest_pl),         intent(in) :: pl      !! Swiftest massive body object
         class(swiftest_body),       intent(in) :: encbody !! Encountering body - Swiftest generic body object (pl or tp) 
         class(swiftest_parameters), intent(in) :: param   !! Current run configuration parameters 
      end subroutine encounter_io_write_list

      module subroutine encounter_setup_aabb(self, n, n_last)
         implicit none
         class(encounter_bounding_box), intent(inout) :: self   !! Swiftest encounter structure
         integer(I4B),                  intent(in)    :: n      !! Number of objects with bounding box extents
         integer(I4B),                  intent(in)    :: n_last !! Number of objects with bounding box extents the previous time this was called
      end subroutine encounter_setup_aabb

      module subroutine encounter_setup_list(self, n)
         implicit none
         class(encounter_list), intent(inout) :: self !! Swiftest encounter structure
         integer(I8B),          intent(in)    :: n    !! Number of encounters to allocate space for
      end subroutine encounter_setup_list

      module subroutine encounter_util_append_list(self, source, lsource_mask)
         implicit none
         class(encounter_list), intent(inout) :: self         !! Swiftest encounter list object
         class(encounter_list), intent(in)    :: source       !! Source object to append
         logical, dimension(:), intent(in)    :: lsource_mask !! Logical mask indicating which elements to append to
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

      module subroutine encounter_util_dealloc_list(self)
         implicit none
         class(encounter_list), intent(inout) :: self !! Swiftest encounter list object
      end subroutine encounter_util_dealloc_list

      module subroutine encounter_util_final_aabb(self)
         implicit none
         type(encounter_bounding_box_1D), intent(inout) :: self !!Bounding box structure along a single dimension
      end subroutine encounter_util_final_aabb

      module subroutine encounter_util_final_list(self)
         implicit none
         type(encounter_list), intent(inout) :: self !! Swiftest encounter list object
      end subroutine encounter_util_final_list

      module subroutine encounter_util_resize_list(self, nnew)
         implicit none
         class(encounter_list), intent(inout) :: self !! Swiftest encounter list 
         integer(I8B),          intent(in)    :: nnew !! New size of list needed
      end subroutine encounter_util_resize_list

      module subroutine encounter_util_spill_list(self, discards, lspill_list, ldestructive)
         implicit none
         class(encounter_list), intent(inout) :: self         !! Swiftest encounter list 
         class(encounter_list), intent(inout) :: discards     !! Discarded object 
         logical, dimension(:), intent(in)    :: lspill_list  !! Logical array of bodies to spill into the discards
         logical,               intent(in)    :: ldestructive !! Logical flag indicating whether or not this operation should alter body by removing the discard list
      end subroutine encounter_util_spill_list

   end interface

end module encounter_classes

