submodule (swiftest_classes) s_setup
   use swiftest
contains

   module subroutine setup_construct_system(system, param)
      !! author: David A. Minton
      !!
      !! Constructor for a Swiftest nbody system. Creates the nbody system object based on the user-input integrator
      !! 
      implicit none
      ! Arguments
      class(swiftest_nbody_system),  allocatable,  intent(inout) :: system     !! Swiftest system object
      class(swiftest_parameters),                  intent(in)    :: param     !! Swiftest parameters

      select case(param%integrator)
      case (BS)
         write(*,*) 'Bulirsch-Stoer integrator not yet enabled'
      case (HELIO)
         allocate(helio_nbody_system :: system)
         select type(system)
         class is (helio_nbody_system)
            allocate(helio_cb :: system%cb)
            allocate(helio_pl :: system%pl)
            allocate(helio_tp :: system%tp)
            allocate(helio_tp :: system%tp_discards)
         end select
      case (RA15)
         write(*,*) 'Radau integrator not yet enabled'
      case (TU4)
         write(*,*) 'TU4 integrator not yet enabled'
      case (WHM)
         allocate(whm_nbody_system :: system)
         select type(system)
         class is (whm_nbody_system)
            allocate(whm_cb :: system%cb)
            allocate(whm_pl :: system%pl)
            allocate(whm_tp :: system%tp)
            allocate(whm_tp :: system%tp_discards)
         end select
      case (RMVS)
         allocate(rmvs_nbody_system :: system)
         select type(system)
         class is (rmvs_nbody_system)
            allocate(rmvs_cb :: system%cb)
            allocate(rmvs_pl :: system%pl)
            allocate(rmvs_tp :: system%tp)
            allocate(rmvs_tp :: system%tp_discards)
         end select
      case (SYMBA)
         allocate(symba_nbody_system :: system)
         select type(system)
         class is (symba_nbody_system)
            allocate(symba_cb :: system%cb)
            allocate(symba_pl :: system%pl)
            allocate(symba_tp :: system%tp)
            allocate(symba_pl :: system%pl_discards)
            allocate(symba_tp :: system%tp_discards)
            allocate(symba_pl :: system%mergeadd_list)
            allocate(symba_pl :: system%mergesub_list)
            allocate(symba_plplenc :: system%plplenc_list)
            allocate(symba_pltpenc :: system%pltpenc_list)
         end select
      case (RINGMOONS)
         write(*,*) 'RINGMOONS-SyMBA integrator not yet enabled'
      case default
         write(*,*) 'Unkown integrator',param%integrator
         call util_exit(FAILURE)
      end select

      return
   end subroutine setup_construct_system


   module subroutine setup_initialize_system(self, param)
      !! author: David A. Minton
      !!
      !! Wrapper method to initialize a basic Swiftest nbody system from files
      !!
      implicit none
      ! Arguments
      class(swiftest_nbody_system), intent(inout) :: self   !! Swiftest system object
      class(swiftest_parameters),   intent(inout) :: param  !! Current run configuration parameters
  
      call self%cb%initialize(param)
      call self%pl%initialize(param)
      call self%tp%initialize(param)
      call util_valid(self%pl, self%tp)
      call self%set_msys()
      call self%pl%set_mu(self%cb) 
      call self%tp%set_mu(self%cb) 
      call self%pl%eucl_index()
      if (.not.param%lrhill_present) call self%pl%set_rhill(self%cb)
      return
   end subroutine setup_initialize_system


   module subroutine setup_body(self, n, param)
      !! author: David A. Minton
      !!
      !! Constructor for base Swiftest particle class. Allocates space for all particles and
      !! initializes all components with a value.
      !! Note: Timing tests indicate that (NDIM, n) is more efficient than (NDIM, n) 
      implicit none
      ! Arguments
      class(swiftest_body),       intent(inout) :: self  !! Swiftest generic body object
      integer(I4B),               intent(in)    :: n     !! Number of particles to allocate space for
      class(swiftest_parameters), intent(in)    :: param !! Current run configuration parameter

      self%nbody = n
      if (n <= 0) return
      self%lfirst = .true.

      allocate(self%id(n))
      allocate(self%name(n))
      allocate(self%status(n))
      allocate(self%ldiscard(n))
      allocate(self%xh(NDIM, n))
      allocate(self%vh(NDIM, n))
      allocate(self%xb(NDIM, n))
      allocate(self%vb(NDIM, n))
      allocate(self%ah(NDIM, n))
      allocate(self%ir3h(n))
      allocate(self%mu(n))
      allocate(self%lmask(n))

      self%id(:)   = 0
      self%name(:) = "UNNAMED"
      self%status(:) = INACTIVE
      self%lmask(:)  = .false.
      self%ldiscard(:) = .false.
      self%xh(:,:)   = 0.0_DP
      self%vh(:,:)   = 0.0_DP
      self%xb(:,:)   = 0.0_DP
      self%vb(:,:)   = 0.0_DP
      self%ah(:,:)   = 0.0_DP
      self%ir3h(:)   = 0.0_DP
      self%mu(:)     = 0.0_DP

      if (param%loblatecb) then
         allocate(self%aobl(NDIM, n))
         self%aobl(:,:) = 0.0_DP
      end if
      if (param%ltides) then
         allocate(self%atide(NDIM, n))
         self%atide(:,:) = 0.0_DP
      end if
      if (param%lgr) then
         allocate(self%agr(NDIM, n))
         self%agr(:,:) = 0.0_DP
      end if

      return
   end subroutine setup_body


   module subroutine setup_pl(self, n, param)
      !! author: David A. Minton
      !!
      !! Constructor for base Swiftest massive body class. Allocates space for all particles and
      !! initializes all components with a value. 
      implicit none
      ! Arguments
      class(swiftest_pl),         intent(inout) :: self  !! Swiftest massive body object
      integer(I4B),               intent(in)    :: n     !! Number of particles to allocate space for
      class(swiftest_parameters), intent(in)    :: param !! Current run configuration parameter

      !> Call allocation method for parent class
      !> The parent class here is the abstract swiftest_body class, so we can't use the type-bound procedure
      call setup_body(self, n, param)
      if (n <= 0) return 

      allocate(self%mass(n))
      allocate(self%Gmass(n))
      allocate(self%rhill(n))

      self%mass(:) = 0.0_DP
      self%Gmass(:) = 0.0_DP
      self%rhill(:) = 0.0_DP

      self%nplpl = 0   

      if (param%lclose) then
         allocate(self%radius(n))
         allocate(self%density(n))
         self%radius(:) = 0.0_DP
         self%density(:) = 1.0_DP
      end if

      if (param%lrotation) then
         allocate(self%rot(NDIM, n))
         allocate(self%Ip(NDIM, n))
         self%rot(:,:) = 0.0_DP
         self%Ip(:,:) = 0.0_DP
      end if

      if (param%ltides) then
         allocate(self%k2(n))
         allocate(self%Q(n))
         allocate(self%tlag(n))
         self%k2(:) = 0.0_DP
         self%Q(:) = 0.0_DP
         self%tlag(:) = 0.0_DP
      end if
      
      return
   end subroutine setup_pl
   

   module subroutine setup_tp(self, n, param)
      !! author: David A. Minton
      !!
      !! Constructor for base Swiftest test particle particle class. Allocates space for 
      !! all particles and initializes all components with a value. 
      implicit none
      ! Arguments
      class(swiftest_tp),         intent(inout) :: self  !! Swiftest test particle object
      integer(I4B),               intent(in)    :: n     !! Number of particles to allocate space for
      class(swiftest_parameters), intent(in)    :: param !! Current run configuration parameter

      !> Call allocation method for parent class
      !> The parent class here is the abstract swiftest_body class, so we can't use the type-bound procedure
      call setup_body(self, n, param)
      if (n <= 0) return

      allocate(self%isperi(n))
      allocate(self%peri(n))
      allocate(self%atp(n))

      self%isperi(:) = 0
      self%peri(:)   = 0.0_DP
      self%atp(:)    = 0.0_DP

      return
   end subroutine setup_tp

end submodule s_setup
