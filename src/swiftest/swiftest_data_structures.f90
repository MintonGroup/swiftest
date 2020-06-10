module swiftest_data_structures
   !! author: The Purdue Swiftest Team -  David A. Minton, Carlisle A. Wishard, Jennifer L.L. Pouplin, and Jacob R. Elliott
   !!
   !! Definition of data and structures generic to all integrators.
   !! Adapted from David E. Kaufmann's Swifter modules: module_swifter.f90
   use swiftest_globals
   use user
   implicit none

   type swiftest_pl
      integer(I4B)                                :: npl    !! Number of massive bodies
      integer(I4B), dimension(:),     allocatable :: name   !! External identifier (hash)
      integer(I4B), dimension(:),     allocatable :: status !! Status
      real(DP),     dimension(:),     allocatable :: mass   !! Mass
      real(DP),     dimension(:),     allocatable :: radius !! Radius
      real(DP),     dimension(:),     allocatable :: rhill  !! Hill's sphere radius
      real(DP),     dimension(:,:),   allocatable :: xh     !! Heliocentric position
      real(DP),     dimension(:,:),   allocatable :: vh     !! Heliocentric velocity
      real(DP),     dimension(:,:),   allocatable :: xb     !! Barycentric position
      real(DP),     dimension(:,:),   allocatable :: vb     !! Barycentric velocity
   contains
      procedure :: set_initial_conditions => swiftest_init_pl
      procedure :: read_from_file => swiftest_read_pl_in 
   end type swiftest_pl

   type swiftest_tp
      integer(I4B)                                :: ntp    !! Number of test particles
      integer(I4B), dimension(:),     allocatable :: name   !! External identifier (hash)
      integer(I4B), dimension(:),     allocatable :: status !! Status
      integer(I4B), dimension(:),     allocatable :: isperi !! Perihelion passage flag
      real(DP),     dimension(:),     allocatable :: peri   !! Perihelion distance
      real(DP),     dimension(:),     allocatable :: atp    !! Semimajor axis following perihelion passage
      real(DP),     dimension(:,:),   allocatable :: xh     !! Heliocentric position
      real(DP),     dimension(:,:),   allocatable :: vh     !! Heliocentric velocity
      real(DP),     dimension(:,:),   allocatable :: xb     !! Barycentric position
      real(DP),     dimension(:,:),   allocatable :: vb     !! Barycentric velocity
   end type swiftest_tp

   interface
      module subroutine swiftest_init_pl(swiftest_plA) 
         class(swiftest_pl), intent(inout)  :: swiftest_plA  !! Swiftest data structure to store massive body initial conditions
      end subroutine swiftest_init_pl

      module subroutine swiftest_read_pl_in(swiftest_plA,param) 
         class(swiftest_pl), intent(inout)  :: swiftest_plA  !! Swiftest data structure to store massive body initial conditions
         type(user_input_parameters),intent(in) :: param    !! Input collection of user-defined parameters
      end subroutine swiftest_read_pl_in
   end interface



end module swiftest_data_structures
