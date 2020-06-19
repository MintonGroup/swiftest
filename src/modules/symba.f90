module symba
   !! author: The Purdue Swiftest Team - David A. Minton, Carlisle A. Wishard, Jennifer L.L. Pouplin, and Jacob R. Elliott
   !!
   !! Definition of classes and methods specific to the Symplectic Massive Body Algorithm
   !!
   !! Adapted from David E. Kaufmann's Swifter modules: symba.f90
   use swiftest_globals
   use helio
   implicit none

   integer(I4B), public, parameter :: NENMAX = 32767
   integer(I4B), public, parameter :: NTENC = 3
   real(DP), public, parameter     :: RHSCALE = 6.5_DP
   real(DP), public, parameter     :: RSHELL = 0.48075_DP

   !********************************************************************************************************************************
   !                                    symba_tp class definitions and method interfaces
   !*******************************************************************************************************************************

   !! SyMBA test particle class
   type, public, extends(helio_pl) :: symba_tp
      integer(I4B), dimension(:),     allocatable :: nplenc  !! Number of encounters with massive bodies this time step
      integer(I4B), dimension(:),     allocatable :: levelg  !! Level at which this particle should be moved
      integer(I4B), dimension(:),     allocatable :: levelm  !! Deepest encounter level achieved this time step
   contains
      procedure, public :: alloc => symba_allocate_tp
      final :: symba_deallocate_tp
   end type symba_tp

   !********************************************************************************************************************************
   !                                    symba_pl class definitions and method interfaces
   !*******************************************************************************************************************************

   !! SyMBA massive body particle class
   type, public, extends(symba_tp) :: symba_pl
      real(DP)                                  :: eoffset      !! Energy offset (net energy lost in mergers)
      logical, dimension(:),        allocatable :: lmerged      !! Flag indicating whether body has merged with another this time step
      integer(I4B), dimension(:),   allocatable :: ntpenc       !! Number of encounters with test particles this time step
      integer(I4B), dimension(:),   allocatable :: nchild       !! Number of children in merger list
      integer(I4B), dimension(:),   allocatable :: index_parent !! Position of the parent of id
      integer(I4B), dimension(:,:), allocatable :: index_child  !! Position of the children of id
   contains
      procedure, public :: alloc => symba_allocate_pl
      final :: symba_deallocate_pl
   end type symba_pl

   !********************************************************************************************************************************
   !                                    symba_encounter class definitions and method interfaces
   !*******************************************************************************************************************************


   !! Generic abstract class structure for a SyMBA encounter class
   type, private, extends(swiftest_body) :: symba_encounter
      logical     , dimension(:),     allocatable :: lvdotr !! Relative vdotr flag
      integer(I4B), dimension(:),     allocatable :: level  !! Encounter recursion level
   contains
      procedure :: alloc => symba_allocate_encounter
      procedure :: set_from_file => symba_encounter_dummy_input
      final :: symba_deallocate_encounter
   end type symba_encounter

   !********************************************************************************************************************************
   !                                    symba_plplenc class definitions and method interfaces
   !*******************************************************************************************************************************

   !! Class structure for a massive body-massive body encounter
   type, public, extends(symba_encounter) :: symba_plplenc
      integer(I4B), dimension(:), allocatable :: index1       !! Position of the first massive body in encounter
      integer(I4B), dimension(:), allocatable :: index2       !! Position of the second massive body in encounter
      integer(I4B), dimension(:), allocatable :: enc_child    !! The child of the encounter
      integer(I4B), dimension(:), allocatable :: enc_parent   !! The child of the encounter
   contains
      procedure :: alloc => symba_allocate_plplenc
      final :: symba_deallocate_plplenc
   end type symba_plplenc

   !********************************************************************************************************************************
   !                                    symba_pltpenc class definitions and method interfaces
   !*******************************************************************************************************************************

   !! Class structure for a massive body-test particle encounter
   type, public, extends(symba_encounter) :: symba_pltpenc
      integer(I4B), dimension(:), allocatable :: indexpl    !! Index position within the main symba structure for the first massive body in an encounter
      integer(I4B), dimension(:), allocatable :: indextp    !! Index position within the main symba structure for the second massive body in an encounter
   contains
      procedure :: alloc => symba_pltpenc_allocate
      final :: symba_deallocate_pltpenc
   end type symba_pltpenc

   !********************************************************************************************************************************
   !                                    symba_merger class definitions and method interfaces
   !********************************************************************************************************************************

   !! Class structure for merger structure
   type, public, extends(swiftest_pl) :: symba_merger
      integer(I4B), dimension(:), allocatable :: index_ps !! Index position within the main symba structure for the body being merged
      integer(I4B), dimension(:), allocatable :: ncomp    !! Number of component bodies in this one during this merger
   contains
      procedure :: alloc => symba_allocate_merger
      final :: symba_deallocate_merger
   end type symba_merger 

!> Only the constructor and destructor method implementations are listed here. All other methods are implemented in the symba submodules.
interface
!! Interfaces for all helio particle methods that are implemented in separate submodules 

   module subroutine io_discard_write_symba(t, mtiny, npl, ntp, nsppl, nsptp, nmergeadd, nmergesub, symba_plA, & 
      discard_plA, discard_tpA, mergeadd_list, mergesub_list, fname, lbig_discard)
      implicit none
      logical , intent(in)       :: lbig_discard
      integer(I4B), intent(in)       :: npl, ntp, nsppl, nsptp, nmergeadd, nmergesub
      real(DP), intent(in)         :: t, mtiny
      character(*), intent(in)       :: fname
      type(symba_pl), intent(inout)      :: symba_plA
      type(swiftest_tp), intent(inout)      :: discard_tpA
      type(swiftest_pl), intent(inout)      :: discard_plA
      type(symba_merger), intent(inout)     :: mergeadd_list, mergesub_list
   end subroutine io_discard_write_symba

   module subroutine symba_casedisruption (t, dt, index_enc, nmergeadd, nmergesub, mergeadd_list, mergesub_list, eoffset, vbs, & 
      symba_plA, nplplenc, plplenc_list, fragmax, mres, rres, m1, m2, rad1, rad2, x1, x2, v1, v2, config)
      implicit none
      integer(I4B), intent(in)       :: index_enc
      integer(I4B), intent(inout)        :: nmergeadd, nmergesub, nplplenc, fragmax
      real(DP), intent(in)          :: t, dt
      real(DP), intent(inout)        :: eoffset, m1, m2, rad1, rad2
      real(DP), dimension(:), intent(inout)     :: mres, rres
      real(DP), dimension(:), intent(in)     :: vbs
      real(DP), dimension(:), intent(inout)    :: x1, x2, v1, v2
      type(symba_plplenc), intent(inout)      :: plplenc_list
      type(symba_merger), intent(inout)     :: mergeadd_list, mergesub_list
      type(symba_pl), intent(inout)      :: symba_plA
      type(swiftest_configuration)        :: config

   end subroutine symba_casedisruption

   module subroutine symba_casehitandrun (t, dt, index_enc, nmergeadd, nmergesub, mergeadd_list, mergesub_list, eoffset, vbs, & 
      symba_plA, nplplenc, plplenc_list, fragmax, mres, rres, m1, m2, rad1, rad2, x1, x2, v1, v2, config)
      implicit none
      integer(I4B), intent(in)       :: index_enc
      integer(I4B), intent(inout)        :: nmergeadd, nmergesub, nplplenc, fragmax
      real(DP), intent(in)          :: t, dt
      real(DP), intent(inout)        :: eoffset, m1, m2, rad1, rad2
      real(DP), dimension(:), intent(inout)     :: mres, rres
      real(DP), dimension(:), intent(in)     :: vbs
      real(DP), dimension(:), intent(inout)    :: x1, x2, v1, v2
      type(symba_plplenc), intent(inout)      :: plplenc_list
      type(symba_merger), intent(inout)     :: mergeadd_list, mergesub_list
      type(symba_pl), intent(inout)      :: symba_plA

   end subroutine symba_casehitandrun

   module subroutine symba_casemerge (t, index_enc, nmergeadd, nmergesub, mergeadd_list, mergesub_list, eoffset, vbs, npl, &
      symba_plA, nplplenc, plplenc_list, array_index1_child, array_index2_child, m1, m2, rad1, rad2, x1, x2, v1, v2)
      implicit none
      integer(I4B), intent(in)       :: index_enc
      integer(I4B), intent(inout)        :: npl, nmergeadd, nmergesub, nplplenc
      real(DP), intent(in)          :: t
      real(DP), intent(inout)        :: eoffset, m1, m2, rad1, rad2
      real(DP), dimension(:), intent(in)     :: vbs
      real(DP), dimension(:), intent(inout)    :: x1, x2, v1, v2
      type(symba_plplenc), intent(inout)      :: plplenc_list
      type(symba_merger), intent(inout)     :: mergeadd_list, mergesub_list
      type(symba_pl), intent(inout)      :: symba_plA
      integer(I4B), dimension(npl), intent(inout)   :: array_index1_child, array_index2_child
   end subroutine symba_casemerge

   module subroutine symba_caseresolve (t, dt, index_enc, nmergeadd, nmergesub, mergeadd_list, mergesub_list, &
      eoffset, vbs, & 
      npl, symba_plA, nplplenc, plplenc_list, regime, config%nplmax, config%ntpmax, fragmax, mres, rres, array_index1_child, &
      array_index2_child, m1, m2, rad1, rad2, x1, x2, v1, v2)
      implicit none
      integer(I4B), intent(in)       :: index_enc, config%nplmax, config%ntpmax
      integer(I4B), intent(inout)        :: npl, nmergeadd, nmergesub, nplplenc, fragmax
      real(DP), intent(in)          :: t, dt
      real(DP), intent(inout)        :: eoffset, m1, m2, rad1, rad2
      real(DP), dimension(:), intent(inout)     :: mres, rres
      real(DP), dimension(:), intent(in)     :: vbs
      real(DP), dimension(:), intent(inout)    :: x1, x2, v1, v2
      type(symba_plplenc), intent(inout)      :: plplenc_list
      type(symba_merger), intent(inout)     :: mergeadd_list, mergesub_list
      type(symba_pl), intent(inout)      :: symba_plA
      integer(I4B), intent(in)       :: regime
      integer(I4B), dimension(npl), intent(inout)   :: array_index1_child, array_index2_child
   end subroutine symba_caseresolve

   module subroutine symba_casesupercatastrophic (t, dt, index_enc, nmergeadd, nmergesub, mergeadd_list, mergesub_list, &
      eoffset, vbs, & 
      symba_plA, nplplenc, plplenc_list, config%nplmax, config%ntpmax, fragmax, mres, rres, m1, m2, rad1, rad2, x1, x2, v1, v2)
      implicit none
      integer(I4B), intent(in)       :: index_enc, config%nplmax, config%ntpmax
      integer(I4B), intent(inout)        :: nmergeadd, nmergesub, nplplenc, fragmax
      real(DP), intent(in)          :: t, dt
      real(DP), intent(inout)        :: eoffset, m1, m2, rad1, rad2
      real(DP), dimension(:), intent(inout)     :: mres, rres
      real(DP), dimension(:), intent(in)     :: vbs
      real(DP), dimension(:), intent(inout)    :: x1, x2, v1, v2
      type(symba_plplenc), intent(inout)      :: plplenc_list
      type(symba_merger), intent(inout)     :: mergeadd_list, mergesub_list
      TYPE(symba_pl), INTENT(INOUT)      :: symba_plA
   end subroutine symba_casesupercatastrophic

   module subroutine symba_chk(xr, vr, rhill1, rhill2, dt, irec, lencounter, lvdotr)
      implicit none
      logical , intent(out)   :: lencounter, lvdotr
      integer(I4B), intent(in)    :: irec
      real(DP), intent(in)      :: rhill1, rhill2, dt
      real(DP), dimension(:), intent(in) :: xr, vr
   end subroutine symba_chk

   module subroutine symba_chk_eucl(num_encounters, k_plpl, symba_plA, dt, lencounter, lvdotr, nplplenc)
      implicit none
      type(symba_pl), intent(in)      :: symba_plA
      integer(I4B), dimension(num_encounters), intent(out) :: lencounter, lvdotr
      integer(I4B), intent(in)    :: num_encounters
      integer(I4B), dimension(2,num_encounters),intent(in)   :: k_plpl
      real(DP), intent(in)      :: dt
      integer(I4B), intent(inout)   :: nplplenc
   end subroutine symba_chk_eucl

   module subroutine symba_chk_eucl_pltp(num_encounters, k_pltp, symba_plA, symba_tpA, dt, lencounter, lvdotr, npltpenc)
      implicit none
      type(symba_pl), intent(in)      :: symba_plA
      type(symba_tp), intent(in)      :: symba_tpA
      integer(I4B), dimension(num_encounters), intent(out) :: lencounter, lvdotr
      integer(I4B), intent(in)    :: num_encounters
      integer(I4B), dimension(2,num_encounters),intent(in)   :: k_pltp
      real(DP), intent(in)      :: dt
      integer(I4B), intent(inout)   :: npltpenc
   end subroutine symba_chk_eucl_pltp

   module subroutine symba_discard_merge_pl(t, npl, symba_plA, nplplenc, plplenc_list)
      implicit none
      integer(I4B), intent(in)        :: nplplenc
      integer(I4B), intent(inout)       :: npl
      real(DP), intent(in)        :: t
      type(symba_pl)         :: symba_plA
      type(symba_plplenc), intent(in)      :: plplenc_list
   end subroutine symba_discard_merge_pl

   module subroutine symba_discard_peri_pl(t, npl, symba_plA, msys, qmin, qmin_alo, qmin_ahi, qmin_coord, ldiscards)
      implicit none
      logical , intent(inout) :: ldiscards
      integer(I4B), intent(in)   :: npl
      real(DP), intent(in)   :: t, msys, qmin, qmin_alo, qmin_ahi
      character(*), intent(in)   :: qmin_coord
      type(symba_pl), intent(inout)   :: symba_plA
   end subroutine symba_discard_peri_pl

   module subroutine symba_discard_pl(t, npl, config%nplmax, nsp, symba_plA, rmin, rmax, config%rmaxu, qmin, qmin_coord,   &
      qmin_alo, qmin_ahi, config%j2rp2, config%j4rp4, eoffset)
      implicit none
      integer(I4B), intent(in)   :: config%nplmax
      integer(I4B), intent(inout) :: npl, nsp
      real(DP), intent(in)   :: t, rmin, rmax, config%rmaxu, qmin, qmin_alo, qmin_ahi, config%j2rp2, config%j4rp4
      real(DP), intent(inout)   :: eoffset
      character(*), intent(in)   :: qmin_coord
      type(symba_pl), intent(inout)   :: symba_plA
   end subroutine symba_discard_pl

   module subroutine symba_discard_sun_pl(t, npl, msys, swiftest_plA, rmin, rmax, config%rmaxu, ldiscards)
      implicit none
      logical , intent(inout) :: ldiscards
      integer(I4B), intent(in)   :: npl
      real(DP), intent(in)   :: t, msys, rmin, rmax, config%rmaxu
      type(swiftest_pl), intent(inout)   :: swiftest_plA
   end subroutine symba_discard_sun_pl

   module subroutine symba_discard_tp(t, npl, ntp, nsp, symba_plA, symba_tpA, dt, &
      rmin, rmax, config%rmaxu, qmin, qmin_coord, qmin_alo, qmin_ahi, lclose, lrhill_present)
      implicit none
      logical , intent(in)   :: lclose, lrhill_present
      integer(I4B), intent(in)   :: npl
      integer(I4B), intent(inout) :: ntp, nsp
      real(DP), intent(in)   :: t, dt, rmin, rmax, config%rmaxu, qmin, qmin_alo, qmin_ahi
      character(*), intent(in)   :: qmin_coord
      type(symba_pl), intent(inout)   :: symba_plA
      type(symba_tp), intent(inout)   :: symba_tpA
   end subroutine symba_discard_tp

   module subroutine symba_energy(npl, config%nplmax, swiftest_plA, config%j2rp2, config%j4rp4, ke, pe, te, htot)
      implicit none
      integer(I4B), intent(in)      :: npl, config%nplmax
      real(DP), intent(in)       :: config%j2rp2, config%j4rp4
      real(DP), intent(out)      :: ke, pe, te
      real(DP), dimension(NDIM), intent(out) :: htot
      type(swiftest_pl), intent(inout)      :: swiftest_plA
   end subroutine symba_energy

   module subroutine symba_fragmentation(t, dt, index_enc, nmergeadd, nmergesub, mergeadd_list, &
      mergesub_list, eoffset, vbs, encounter_file, out_type, npl, ntp, &
      symba_plA, symba_tpA, nplplenc, npltpenc, pltpenc_list, plplenc_list, &
      config%nplmax, config%ntpmax, fragmax)
      implicit none
      integer(I4B), intent(in)       :: index_enc, config%nplmax, config%ntpmax
      integer(I4B), intent(inout)        :: nmergeadd, nmergesub, nplplenc, npltpenc, fragmax
      integer(I4B), intent(inout)        :: npl, ntp
      real(DP), intent(in)          :: t, dt
      real(DP), intent(inout)        :: eoffset
      real(DP), dimension(NDIM), intent(in)     :: vbs
      character(*), intent(in)       :: encounter_file, out_type
      type(symba_plplenc), intent(inout)      :: plplenc_list
      type(symba_pltpenc), intent(inout)      :: pltpenc_list
      type(symba_merger), intent(inout)     :: mergeadd_list, mergesub_list
      type(symba_pl), intent(inout)      :: symba_plA
      type(symba_tp), intent(inout)      :: symba_tpA
   end subroutine symba_fragmentation

   module subroutine symba_getacch(lextra_force, t, npl, nplm, config%nplmax, symba_plA, config%j2rp2, config%j4rp4, nplplenc, &
      plplenc_list)
      implicit none
      logical , intent(in)        :: lextra_force
      integer(I4B), intent(in)        :: npl, nplm, config%nplmax, nplplenc
      real(DP), intent(in)        :: t, config%j2rp2, config%j4rp4
      type(symba_pl), intent(inout)      :: symba_plA
      type(symba_plplenc), intent(in)      :: plplenc_list
   end subroutine symba_getacch

   module subroutine symba_getacch_tp(lextra_force, t, npl, nplm, config%nplmax, ntp, config%ntpmax, symba_plA, symba_tpA, &
      xh, config%j2rp2, config%j4rp4,  &
      npltpenc, pltpenc_list)
      implicit none
      logical , intent(in)        :: lextra_force
      integer(I4B), intent(in)        :: npl, nplm, config%nplmax, ntp, config%ntpmax, npltpenc
      real(DP), intent(in)        :: t, config%j2rp2, config%j4rp4
      real(DP), dimension(NDIM, npl), intent(in)   :: xh
      type(symba_pl), intent(inout)      :: symba_plA
      type(symba_tp), intent(inout)      :: symba_tpA
      type(symba_pltpenc), intent(in)      :: pltpenc_list
   end subroutine symba_getacch_tp

   module subroutine symba_getacch_eucl(lextra_force, t, npl, nplm, config%nplmax, symba_plA, config%j2rp2, config%j4rp4, nplplenc, &
      plplenc_list, num_plpl_comparisons, k_plpl)
      implicit none
      logical , intent(in)        :: lextra_force
      integer(I4B), intent(in)        :: npl, nplm, config%nplmax, nplplenc, num_plpl_comparisons
      real(DP), intent(in)        :: t, config%j2rp2, config%j4rp4
      type(symba_pl), intent(inout)      :: symba_plA
      type(symba_plplenc), intent(in)      :: plplenc_list
      integer(I4B), dimension(num_plpl_comparisons,2),intent(in) :: k_plpl
   end subroutine symba_getacch_eucl

   module subroutine symba_getacch_tp_eucl(lextra_force, t, npl, nplm, config%nplmax, ntp, config%ntpmax, symba_plA, symba_tpA, &
      xh, config%j2rp2, config%j4rp4, npltpenc, pltpenc_list, num_pltp_comparisons, k_pltp)
      implicit none
      logical , intent(in)        :: lextra_force
      integer(I4B), intent(in)        :: npl, nplm, config%nplmax, ntp, config%ntpmax, npltpenc, num_pltp_comparisons
      real(DP), intent(in)        :: t, config%j2rp2, config%j4rp4
      real(DP), dimension(NDIM, npl), intent(in)   :: xh
      type(symba_pl), intent(inout)      :: symba_plA
      type(symba_tp), intent(inout)      :: symba_tpA
      type(symba_pltpenc), intent(in)      :: pltpenc_list
      integer(I4B), dimension(num_pltp_comparisons,2), intent(in) :: k_pltp
   end subroutine symba_getacch_tp_eucl

   module subroutine symba_helio_drift(irec, npl, symba_plA, dt)
      implicit none
      integer(I4B), intent(in)   :: irec, npl
      real(DP), intent(in)   :: dt
      type(symba_pl), intent(inout) :: symba_plA
   end subroutine symba_helio_drift

   module subroutine symba_helio_drift_tp(irec, ntp, symba_tpA, mu, dt)
      implicit none
      integer(I4B), intent(in)   :: irec, ntp
      real(DP), intent(in)   :: mu, dt
      type(symba_tp), intent(inout) :: symba_tpA
   end subroutine symba_helio_drift_tp

   module subroutine symba_helio_getacch(lflag, lextra_force, t, npl, nplm, config%nplmax, helio_plA, config%j2rp2, config%j4rp4)
      implicit none
      logical , intent(in)   :: lflag, lextra_force
      integer(I4B), intent(in)   :: npl, nplm, config%nplmax
      real(DP), intent(in)    :: t, config%j2rp2, config%j4rp4
      type(helio_pl), intent(inout)  :: helio_plA
   end subroutine symba_helio_getacch

   module subroutine symba_helio_getacch_int(npl, nplm, helio_plA)
      implicit none
      integer(I4B), intent(in)   :: npl, nplm
      type(helio_pl), intent(inout) :: helio_plA
   end subroutine symba_helio_getacch_int

   module subroutine symba_kick(irec, nplplenc, npltpenc, plplenc_list, pltpenc_list, dt, sgn, symba_plA, &
      symba_tpA)
      implicit none
      integer(I4B), intent(in)   :: irec, nplplenc, npltpenc
      real(DP), intent(in)     :: dt, sgn
      type(symba_plplenc), intent(in) :: plplenc_list
      type(symba_pltpenc), intent(in) :: pltpenc_list
      type(symba_pl), intent(inout)   :: symba_plA
      type(symba_tp), intent(inout)   :: symba_tpA
   end subroutine symba_kick

   module subroutine symba_merge_pl(t, dt, index_enc, nplplenc, plplenc_list, nmergeadd, nmergesub, &
      mergeadd_list, mergesub_list, eoffset, vbs, encounter_file, out_type, npl, symba_plA, &
      symba_tpA)
      implicit none
      integer(I4B), intent(in)       :: index_enc, nplplenc
      integer(I4B), intent(inout)        :: nmergeadd, nmergesub, npl
      real(DP), intent(in)          :: t, dt
      real(DP), intent(inout)        :: eoffset
      real(DP), dimension(NDIM), intent(in)     :: vbs
      character(*), intent(in)       :: encounter_file, out_type
      type(symba_plplenc), intent(inout)      :: plplenc_list
      type(symba_merger),  intent(inout)      :: mergeadd_list, mergesub_list
      type(symba_pl), intent(inout)      :: symba_plA
      type(symba_tp), intent(inout)      :: symba_tpA
   end subroutine symba_merge_pl

   module subroutine symba_merge_tp(t, dt, index_enc, npltpenc, pltpenc_list, vbs, encounter_file, out_type, symba_plA, symba_tpA)
      implicit none
      integer(I4B), intent(in)      :: index_enc, npltpenc
      real(DP), intent(in)       :: t, dt
      real(DP), dimension(NDIM), intent(in)  :: vbs
      character(*), intent(in)      :: encounter_file, out_type
      type(symba_pltpenc), intent(inout)   :: pltpenc_list
      type(symba_pl), intent(inout)      :: symba_plA
      type(symba_tp), intent(inout)      :: symba_tpA
   end subroutine symba_merge_tp

   module subroutine symba_peri(lfirst, npl, symba_plA, msys, qmin_coord)
      implicit none
      logical , intent(in)   :: lfirst
      integer(I4B), intent(in)   :: npl
      real(DP), intent(in)    :: msys
      character(*), intent(in)   :: qmin_coord
      type(symba_pl), intent(inout)  :: symba_plA
   end subroutine symba_peri

   module subroutine symba_rearray(t, npl, ntp, nsppl, nsptp, symba_plA, symba_tpA, nmergeadd, &
      mergeadd_list, discard_plA, discard_tpA, config%nplmax, config%j2rp2, config%j4rp4)
      implicit none
      integer(I4B), intent(inout)        :: npl, ntp, nsppl, nsptp, nmergeadd, config%nplmax !change to fragadd
      real(DP), intent(in)          :: t, config%j2rp2, config%j4rp4
      type(symba_pl), intent(inout)      :: symba_plA
      type(symba_tp), intent(inout)      :: symba_tpA
      type(swiftest_tp), intent(inout)      :: discard_tpA
      type(swiftest_pl), intent(inout)      :: discard_plA
      type(symba_merger), intent(inout)     :: mergeadd_list !change to fragadd_list

   end subroutine symba_rearray

   module subroutine symba_reorder_pl(npl, symba_plA)
      implicit none
      integer(I4B), intent(in) :: npl
      type(symba_pl), intent(inout)  :: symba_plA
      integer(I4B)         :: i
      integer(I4B), dimension(:), allocatable   :: index
      real(DP), dimension(:), allocatable   :: mass
      real(DP), dimension(:,:), allocatable   :: symba_plwkspa
      integer(I4B), dimension(:,:), allocatable :: symba_plwkspa_id_status
   end subroutine symba_reorder_pl

   module subroutine symba_setup(npl, ntp, symba_plA, symba_tpA, symba_pl1p, symba_tp1p, swiftest_pl1p, &
      swiftest_tp1p)
      implicit none
      integer(I4B), intent(in)         :: npl, ntp
      type(swiftest_pl), pointer         :: swiftest_pl1p
      type(swiftest_tp), pointer         :: swiftest_tp1p
      type(symba_pl), intent(inout) :: symba_plA
      type(symba_tp), intent(inout) :: symba_tpA
      type(symba_pl), pointer          :: symba_pl1p
      type(symba_tp), pointer          :: symba_tp1p
   end subroutine symba_setup

   !> Method to remove the inactive symba test particles and spill them to a discard object
   module subroutine symba_spill_tp(self,discard)
      implicit none
      class(symba_tp), intent(inout) :: self    !! Swiftest test particle object to input
      class(symba_tp), intent(inout) :: discard !! Discarded body list
   end subroutine symba_spill_tp

   !> Method to remove the inactive symba massive bodies and spill them to a discard object
   module subroutine symba_spill_pl(self,discard)
      implicit none
      class(symba_pl), intent(inout) :: self    !! Swiftest test particle object to input
      class(symba_pl), intent(inout) :: discard !! Discarded body list
   end subroutine symba_spill_pl

   module subroutine symba_step_eucl(lfirst, lextra_force, lclose, t, npl, config%nplmax, ntp, config%ntpmax, symba_plA, symba_tpA, config%j2rp2, config%j4rp4,&
      dt,nplplenc, npltpenc, plplenc_list, pltpenc_list, nmergeadd, nmergesub, mergeadd_list, mergesub_list, eoffset,&
      mtiny,encounter_file, out_type, num_plpl_comparisons, k_plpl, num_pltp_comparisons, k_pltp)
      implicit none
      logical , intent(in)    :: lextra_force, lclose
      logical , intent(inout)   :: lfirst
      integer(I4B), intent(in)    :: npl, config%nplmax, ntp, config%ntpmax
      integer(I4B), intent(inout)   :: nplplenc, npltpenc, nmergeadd, nmergesub
      real(DP), intent(in)      :: t, config%j2rp2, config%j4rp4, dt, mtiny
      real(DP), intent(inout)     :: eoffset
      character(*), intent(in)    :: encounter_file, out_type
      type(symba_pl), intent(inout)   :: symba_plA
      type(symba_tp), intent(inout)   :: symba_tpA
      type(symba_plplenc), intent(inout) :: plplenc_list
      type(symba_pltpenc), intent(inout) :: pltpenc_list
      type(symba_merger), intent(inout)  :: mergeadd_list, mergesub_list
      integer(I4B), intent(in)    :: num_plpl_comparisons, num_pltp_comparisons
      integer(I4B), dimension(2,num_plpl_comparisons),intent(in) :: k_plpl
      integer(I4B), dimension(2,num_pltp_comparisons),intent(in) :: k_pltp
   end subroutine symba_step_eucl

   module subroutine symba_step(lfirst, lextra_force, lclose, t, npl, config%nplmax, ntp, config%ntpmax, symba_plA, &
      symba_tpA, config%j2rp2, config%j4rp4, dt, nplplenc, npltpenc, plplenc_list, pltpenc_list, nmergeadd, &
      nmergesub, mergeadd_list, mergesub_list, eoffset, mtiny, encounter_file, out_type, &
      fragmax)
      implicit none
      logical , intent(in)    :: lextra_force, lclose
      logical , intent(inout)   :: lfirst
      integer(I4B), intent(in)    :: npl, config%nplmax, ntp, config%ntpmax
      integer(I4B), intent(inout)   :: nplplenc, npltpenc, nmergeadd, nmergesub, fragmax
      real(DP), intent(in)      :: t, config%j2rp2, config%j4rp4, dt, mtiny
      real(DP), intent(inout)     :: eoffset
      character(*), intent(in)    :: encounter_file, out_type
      type(symba_pl), intent(inout)   :: symba_plA
      type(symba_tp), intent(inout)   :: symba_tpA
      type(symba_plplenc), intent(inout) :: plplenc_list
      type(symba_pltpenc), intent(inout) :: pltpenc_list
      type(symba_merger), intent(inout)  :: mergeadd_list, mergesub_list
   end subroutine symba_step

   ! for testing purposes only _ use with symba_step_test
   module subroutine symba_step_helio(lfirst, lextra_force, t, npl, nplm, config%nplmax, ntp, config%ntpmax, helio_plA, helio_tpA, config%j2rp2,   &
      config%j4rp4, dt)
      implicit none
      logical , intent(in)   :: lextra_force
      logical , intent(inout)   :: lfirst
      integer(I4B), intent(in)   :: npl, nplm, config%nplmax, ntp, config%ntpmax
      real(DP), intent(in)   :: t, config%j2rp2, config%j4rp4, dt
      type(helio_pl), intent(inout) :: helio_plA
      type(helio_tp), intent(inout) :: helio_tpA
   end subroutine symba_step_helio

   module subroutine symba_step_helio_pl(lfirst, lextra_force, t, npl, nplm, config%nplmax, helio_plA, config%j2rp2, config%j4rp4, dt, xbeg, xend,   &
      ptb, pte)
      implicit none
      logical , intent(in)       :: lextra_force
      logical , intent(inout)      :: lfirst
      integer(I4B), intent(in)       :: npl, nplm, config%nplmax
      real(DP), intent(in)       :: t, config%j2rp2, config%j4rp4, dt
      real(DP), dimension(NDIM, nplm), intent(out) :: xbeg, xend
      real(DP), dimension(NDIM), intent(out)   :: ptb, pte
      type(helio_pl), intent(inout)     :: helio_plA
   end subroutine symba_step_helio_pl

   module subroutine symba_step_interp_eucl(lextra_force, lclose, t, npl, nplm, config%nplmax, ntp, config%ntpmax, symba_plA, symba_tpA, config%j2rp2,&
      config%j4rp4, dt, eoffset, mtiny, nplplenc, npltpenc, plplenc_list, pltpenc_list, nmergeadd, nmergesub, mergeadd_list,&
      mergesub_list, encounter_file, out_type, num_plpl_comparisons, k_plpl, num_pltp_comparisons, k_pltp)
      implicit none
      logical , intent(in)    :: lextra_force, lclose
      integer(I4B), intent(in)    :: npl, nplm, config%nplmax, ntp, config%ntpmax, nplplenc, npltpenc, num_pltp_comparisons
      integer(I4B), intent(inout)   :: nmergeadd, nmergesub
      real(DP), intent(in)      :: t, config%j2rp2, config%j4rp4, dt, mtiny
      real(DP), intent(inout)     :: eoffset
      character(*), intent(in)    :: encounter_file, out_type
      type(symba_pl), intent(inout)   :: symba_plA
      type(symba_tp), intent(inout)   :: symba_tpA
      type(symba_plplenc), intent(inout) :: plplenc_list
      type(symba_pltpenc), intent(inout) :: pltpenc_list
      type(symba_merger), intent(inout)  :: mergeadd_list, mergesub_list
      integer(I4B), intent(in)       :: num_plpl_comparisons
      integer(I4B), dimension(num_plpl_comparisons,2),intent(in) :: k_plpl
      integer(I4B), dimension(2,num_pltp_comparisons),intent(in) :: k_pltp
   end subroutine symba_step_interp_eucl

   module subroutine symba_step_interp(lextra_force, lclose, t, npl, nplm, config%nplmax, ntp, config%ntpmax, symba_plA, symba_tpA, config%j2rp2,   &
      config%j4rp4, dt, eoffset, mtiny, nplplenc, npltpenc, plplenc_list, pltpenc_list, nmergeadd, nmergesub, mergeadd_list,   &
      mergesub_list, encounter_file, out_type, fragmax)
      implicit none
      logical , intent(in)    :: lextra_force, lclose
      integer(I4B), intent(in)    :: npl, nplm, config%nplmax, ntp, config%ntpmax, nplplenc, npltpenc
      integer(I4B), intent(inout)   :: nmergeadd, nmergesub, fragmax
      real(DP), intent(in)      :: t, config%j2rp2, config%j4rp4, dt, mtiny
      real(DP), intent(inout)     :: eoffset
      character(*), intent(in)    :: encounter_file, out_type
      type(symba_pl), intent(inout)   :: symba_plA
      type(symba_tp), intent(inout)   :: symba_tpA
      type(symba_plplenc), intent(inout) :: plplenc_list
      type(symba_pltpenc), intent(inout) :: pltpenc_list
      type(symba_merger), intent(inout)  :: mergeadd_list, mergesub_list
   end subroutine symba_step_interp

   recursive module subroutine symba_step_recur(lclose, t, ireci, npl, nplm, ntp, symba_plA, symba_tpA, dt0, eoffset, nplplenc, &
      npltpenc, plplenc_list, pltpenc_list, nmergeadd, nmergesub, mergeadd_list, mergesub_list, encounter_file, & 
      out_type, config%nplmax, config%ntpmax, fragmax)
      implicit none
      logical , intent(in)    :: lclose
      integer(I4B), intent(in)    :: ireci, npl, nplm, ntp, nplplenc, npltpenc, config%nplmax, config%ntpmax, fragmax
      integer(I4B), intent(inout)   :: nmergeadd, nmergesub
      real(DP), intent(in)      :: t, dt0
      real(DP), intent(inout)     :: eoffset
      character(*), intent(in)    :: encounter_file, out_type
      type(symba_pl), intent(inout)   :: symba_plA
      type(symba_tp), intent(inout)   :: symba_tpA
      type(symba_plplenc), intent(inout) :: plplenc_list
      type(symba_pltpenc), intent(inout) :: pltpenc_list
      type(symba_merger), intent(inout)  :: mergeadd_list, mergesub_list
   end subroutine symba_step_recur

   module subroutine symba_user_getacch(t, npl, symba_plA)
      implicit none
      integer(I4B), intent(in)   :: npl
      real(DP), intent(in)    :: t
      type(symba_pl), intent(inout):: symba_plA
   end subroutine symba_user_getacch

   module subroutine symba_user_getacch_tp(t, ntp, symba_tpA)
      implicit none
      integer(I4B), intent(in)   :: ntp
      real(DP), intent(in)    :: t
      type(symba_tp), intent(inout)  :: symba_tpA
   end subroutine symba_user_getacch_tp

end interface

contains
   !! SyMBA constructor and desctructor methods
   subroutine symba_allocate_tp(self,n)
      !! SyMBA test particle constructor method
      implicit none

      class(symba_tp), intent(inout) :: self !! Symba test particle object
      integer, intent(in)            :: n    !! Number of test particles to allocate

      call self%helio_tp%alloc(n)

      if (self%is_allocated) then
         write(*,*) 'Symba test particle structure already alllocated'
         return
      end if
      write(*,*) 'Allocating the Swiftest test particle structure'
      
      if (n <= 0) return
      allocate(self%nplenc(n))
      allocate(self%levelg(n))
      allocate(self%levelm(n))

      self%nplenc(:) = 0
      self%levelg(:) = 0
      self%levelm(:) = 0
      return
   end subroutine symba_allocate_tp

   subroutine symba_deallocate_tp(self)
      !! SyMBA test particle destructor/finalizer
      implicit none

      type(symba_tp), intent(inout)    :: self

      if (self%is_allocated) then
         deallocate(self%nplenc)
         deallocate(self%levelg)
         deallocate(self%levelm)
      end if
      return
   end subroutine symba_deallocate_tp

   subroutine symba_allocate_pl(self,n)
      !! SyMBA massive body constructor method
      implicit none

      class(symba_pl), intent(inout) :: self !! SyMBA massive body particle object
      integer, intent(in)            :: n    !! Number of massive body particles to allocate

      call self%helio_pl%alloc(n)

      if (self%is_allocated) then
         write(*,*) 'Symba massive body structure already alllocated'
         return
      end if
      if (n <= 0) return
      allocate(self%lmerged(n))
      allocate(self%nplenc(n))
      allocate(self%ntpenc(n))
      allocate(self%levelg(n))
      allocate(self%levelm(n))
      allocate(self%nchild(n))
      allocate(self%index_parent(n))
      allocate(self%index_child(n,n))

      self%lmerged(:) = .false.
      self%nplenc(:) = 0
      self%ntpenc(:) = 0
      self%levelg(:) = 0
      self%levelm(:) = 0
      self%nchild(:) = 0
      self%index_parent(:) = 1
      self%index_child(:,:) = 1

      return
   end subroutine symba_allocate_pl

   subroutine symba_deallocate_pl(self)
      !! SyMBA massive body destructor/finalizer
      implicit none

      type(symba_pl), intent(inout)    :: self
      if (self%is_allocated) then
         deallocate(self%lmerged)
         deallocate(self%nplenc)
         deallocate(self%ntpenc)
         deallocate(self%levelg)
         deallocate(self%levelm)
         deallocate(self%nchild)
         deallocate(self%index_parent)
         deallocate(self%index_child)
      end if
      return
   end subroutine symba_deallocate_pl

   subroutine symba_allocate_encounter(self,n)
      !! Basic Symba encounter structure constructor method
      implicit none

      class(symba_encounter), intent(inout) :: self !! SyMBA encounter super class
      integer, intent(in)                   :: n    !! Number of test particles to allocate
     
      call self%swiftest_body%alloc(n)
      if (n <= 0) return

      if (self%is_allocated) then
         write(*,*) 'SyMBA encounter structure already alllocated'
         return
      end if
      write(*,*) 'Allocating the Symba encounter superclass'

      allocate(self%lvdotr(n))
      allocate(self%level(n))

      self%lvdotr(:) = .false.
      self%level(:) = 0

      return
   end subroutine symba_allocate_encounter
   
   subroutine symba_deallocate_encounter(self)
      !! SyMBA encounter superclass destructor/finalizer
      implicit none

      type(symba_encounter), intent(inout)    :: self

      if (self%is_allocated) then
         deallocate(self%lvdotr)
         deallocate(self%level)
      end if
      return
   end subroutine symba_deallocate_encounter

   subroutine symba_encounter_dummy_input(self,config) 
      !! This method is needed in order to extend the abstract type swiftest_body. It does nothing
      implicit none
      class(symba_encounter), intent(inout)  :: self  !! SyMBA encounter data structure 
      type(swiftest_configuration),intent(in) :: config !! Input collection of user-defined parameters
      return
   end subroutine symba_encounter_dummy_input

   subroutine symba_pltpenc_allocate(self,n)
      !! SyMBA massive body-test particle encounter structure constructor method
      implicit none

      class(symba_pltpenc), intent(inout) :: self !! SyMBA massive body-test particle encounter class
      integer, intent(in)                 :: n    !! Number of encounter slots to allocate

      call self%symba_encounter%alloc(n)
      if (n <= 0) return

      if (self%is_allocated) then
         write(*,*) 'SyMBA pl-tp encounter structure already alllocated'
         return
      end if
      write(*,*) 'Allocating the Symba pl-tp encounter class'

      allocate(self%indexpl(n))
      allocate(self%indextp(n))

      self%indexpl(:) = 1
      self%indextp(:) = 1
   end subroutine symba_pltpenc_allocate

   subroutine symba_deallocate_pltpenc(self)
      !! SyMBA massive body-test particle encounter destructor/finalizer
      implicit none

      type(symba_pltpenc), intent(inout)    :: self

      if (self%is_allocated) then
         deallocate(self%indexpl)
         deallocate(self%indextp)
      end if
      return
   end subroutine symba_deallocate_pltpenc

   subroutine symba_allocate_plplenc(self,n)
      !! SyMBA massive body-massive body particle encounter structure constructor method
      implicit none

      class(symba_plplenc), intent(inout) :: self !! SyMBA massive body-massive body encounter class
      integer, intent(in)                 :: n    !! Number of encounter slots to allocate

      call self%symba_encounter%alloc(n)
      if (n <= 0) return

      if (self%is_allocated) then
         write(*,*) 'SyMBA pl-pl encounter structure already alllocated'
         return
      end if
      write(*,*) 'Allocating the Symba pl-pl encounter class'
    
      allocate(self%index1(n))
      allocate(self%index2(n))
      allocate(self%enc_child(n))
      allocate(self%enc_parent(n))
      
      self%index1(:) = 1
      self%index2(:) = 1
      self%enc_child(:) = 1
      self%enc_parent(:) = 1

      return
   end subroutine symba_allocate_plplenc

   subroutine symba_deallocate_plplenc(self)
      !! SyMBA massive body-massive body encounter destructor/finalizer
      implicit none

      type(symba_plplenc), intent(inout)    :: self

      if (self%is_allocated) then
         deallocate(self%index1)
         deallocate(self%index2)
         deallocate(self%enc_child)
         deallocate(self%enc_parent)
      end if
      return
   end subroutine symba_deallocate_plplenc

   subroutine symba_allocate_merger(self,n)
      !! SyMBA merger encounter structure constructor method
      implicit none

      class(symba_merger), intent(inout) :: self !! SyMBA merger class
      integer, intent(in)                 :: n    !! Number of encounter slots to allocate

      call self%swiftest_pl%alloc(n)
      if (n <= 0) return

      if (self%is_allocated) then
         write(*,*) 'SyMBA merger structure already alllocated'
         return
      end if
      write(*,*) 'Allocating the SyMBA merger  class'
    
      allocate(self%index_ps(n))
      allocate(self%ncomp(n))
      
      self%index_ps(:) = 1
      self%ncomp(:) = 0

   end subroutine symba_allocate_merger

   subroutine symba_deallocate_merger(self)
      !! SyMBA merger destructor/finalizer
      implicit none

      type(symba_merger), intent(inout)    :: self

      if (self%is_allocated) then
         deallocate(self%index_ps)
         deallocate(self%ncomp)
      end if
      return
   end subroutine symba_deallocate_merger

end module symba
