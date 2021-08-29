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
      class(swiftest_parameters),                  intent(inout) :: param     !! Swiftest parameters

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
            allocate(symba_tp :: system%tp_discards)
            allocate(symba_merger :: system%pl_adds)
            allocate(symba_merger :: system%pl_discards)
            allocate(symba_pltpenc :: system%pltpenc_list)
            allocate(symba_plplenc :: system%plplenc_list)
            allocate(symba_plplenc :: system%plplcollision_list)
         end select
      case (RINGMOONS)
         write(*,*) 'RINGMOONS-SyMBA integrator not yet enabled'
      case default
         write(*,*) 'Unkown integrator',param%integrator
         call util_exit(FAILURE)
      end select

      return
   end subroutine setup_construct_system


   module subroutine setup_encounter(self, n)
      !! author: David A. Minton
      !!
      !! A constructor that sets the number of encounters and allocates and initializes all arrays  
      !!
      implicit none
      ! Arguments
      class(swiftest_encounter), intent(inout) :: self !! Swiftest encounter structure
      integer(I4B),              intent(in)    :: n    !! Number of encounters to allocate space for

      self%nenc = n
      if (n == 0) return

      if (allocated(self%lvdotr)) deallocate(self%lvdotr)
      if (allocated(self%status)) deallocate(self%status)
      if (allocated(self%kidx)) deallocate(self%kidx)
      if (allocated(self%index1)) deallocate(self%index1)
      if (allocated(self%index2)) deallocate(self%index2)
      if (allocated(self%id1)) deallocate(self%id1)
      if (allocated(self%id2)) deallocate(self%id2)
      if (allocated(self%x1)) deallocate(self%x1)
      if (allocated(self%x2)) deallocate(self%x2)
      if (allocated(self%v1)) deallocate(self%v1)
      if (allocated(self%v2)) deallocate(self%v2)
      if (allocated(self%t)) deallocate(self%t)

      allocate(self%lvdotr(n))
      allocate(self%status(n))
      allocate(self%kidx(n))
      allocate(self%index1(n))
      allocate(self%index2(n))
      allocate(self%id1(n))
      allocate(self%id2(n))
      allocate(self%x1(NDIM,n))
      allocate(self%x2(NDIM,n))
      allocate(self%v1(NDIM,n))
      allocate(self%v2(NDIM,n))
      allocate(self%t(n))

      self%lvdotr(:) = .false.
      self%status(:) = INACTIVE
      self%kidx(:) = 0_I8B
      self%index1(:) = 0
      self%index2(:) = 0
      self%id1(:) = 0
      self%id2(:) = 0
      self%x1(:,:) = 0.0_DP
      self%x2(:,:) = 0.0_DP
      self%v1(:,:) = 0.0_DP
      self%v2(:,:) = 0.0_DP
      self%t(:) = 0.0_DP

      return
   end subroutine setup_encounter


   module subroutine setup_initialize_particle_info_system(self, param)
      !! author: David A. Minton
      !!
      !! Setup up particle information metadata from initial conditions
      !
      implicit none
      ! Arguments
      class(swiftest_nbody_system), intent(inout) :: self  !! Swiftest nbody system object
      class(swiftest_parameters),   intent(inout) :: param !! Current run configuration parameters
      ! Internals
      logical :: ltrack_origin
      integer(I4B) :: i

      associate(cb => self%cb, pl => self%pl, npl => self%pl%nbody, tp => self%tp, ntp => self%tp%nbody)

         if (param%nciu%ltrack_origin) then
            cb%info%origin_type = "Initial conditions"
            cb%info%origin_time = param%t0
            cb%info%origin_xh(:) = 0.0_DP
            cb%info%origin_vh(:) = 0.0_DP
            do i = 1, self%pl%nbody
               pl%info(i)%origin_type = "Initial conditions"
               pl%info(i)%origin_time = param%t0
               pl%info(i)%origin_xh(:) = self%pl%xh(:,i)
               pl%info(i)%origin_vh(:) = self%pl%vh(:,i)
            end do
            do i = 1, self%tp%nbody
               tp%info(i)%origin_type = "Initial conditions"
               tp%info(i)%origin_time = param%t0
               tp%info(i)%origin_xh(:) = self%tp%xh(:,i)
               tp%info(i)%origin_vh(:) = self%tp%vh(:,i)
            end do
         end if

         cb%info%particle_type = CB_TYPE_NAME
         if (npl > 0) pl%info(1:npl)%particle_type = PL_TYPE_NAME
         if (ntp > 0) tp%info(1:ntp)%particle_type = TP_TYPE_NAME

      end associate

      return
   end subroutine setup_initialize_particle_info_system


   module subroutine setup_initialize_system(self, param)
      !! author: David A. Minton
      !!
      !! Wrapper method to initialize a basic Swiftest nbody system from files
      !!
      implicit none
      ! Arguments
      class(swiftest_nbody_system), intent(inout) :: self   !! Swiftest system object
      class(swiftest_parameters),   intent(inout) :: param  !! Current run configuration parameters
 
      call self%cb%read_in(param)
      call self%pl%read_in(param)
      call self%tp%read_in(param)
      call self%validate_ids(param)
      call self%set_msys()
      call self%pl%set_mu(self%cb) 
      call self%tp%set_mu(self%cb) 
      if (param%in_form == EL) then
         call self%pl%el2xv(self%cb)
         call self%tp%el2xv(self%cb)
      end if
      call self%pl%index(param)
      if (.not.param%lrhill_present) call self%pl%set_rhill(self%cb)
      self%pl%lfirst = param%lfirstkick
      self%tp%lfirst = param%lfirstkick

      if (param%lrestart) then
         call self%read_particle_info(param)
      else
         call self%init_particle_info(param)
      end if
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
      ! Internals
      integer(I4B) :: i
      self%nbody = n
      if (n <= 0) return
      self%lfirst = .true.

      if (allocated(self%id)) deallocate(self%id)
      if (allocated(self%info)) deallocate(self%info)
      if (allocated(self%status)) deallocate(self%status)
      if (allocated(self%ldiscard)) deallocate(self%ldiscard)
      if (allocated(self%xh)) deallocate(self%xh)
      if (allocated(self%vh)) deallocate(self%vh)
      if (allocated(self%xb)) deallocate(self%xb)
      if (allocated(self%vb)) deallocate(self%vb)
      if (allocated(self%ah)) deallocate(self%ah)
      if (allocated(self%ir3h)) deallocate(self%ir3h)
      if (allocated(self%mu)) deallocate(self%mu)
      if (allocated(self%lmask)) deallocate(self%lmask)

      allocate(self%id(n))
      allocate(self%info(n))
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
      self%info(:)%name = "UNNAMED"
      self%info(:)%particle_type = "UKNOWN"
      self%info(:)%origin_type = "UNKNOWN"
      self%info(:)%origin_time = -1.0_DP
      do i = 1, n
         self%info(i)%origin_xh(:) = 0.0_DP
         self%info(i)%origin_vh(:) = 0.0_DP
      end do
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
         if (allocated(self%aobl)) deallocate(self%aobl)
         allocate(self%aobl(NDIM, n))
         self%aobl(:,:) = 0.0_DP
      end if
      if (param%ltides) then
         if (allocated(self%atide)) deallocate(self%lmask)
         allocate(self%atide(NDIM, n))
         self%atide(:,:) = 0.0_DP
      end if
      if (param%lgr) then
         if (allocated(self%agr)) deallocate(self%lmask)
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

      if (allocated(self%info)) deallocate(self%info)
      if (allocated(self%mass)) deallocate(self%mass)
      if (allocated(self%Gmass)) deallocate(self%Gmass)
      if (allocated(self%rhill)) deallocate(self%rhill)

      allocate(swiftest_particle_info :: self%info(n))
      allocate(self%mass(n))
      allocate(self%Gmass(n))
      allocate(self%rhill(n))

      self%mass(:) = 0.0_DP
      self%Gmass(:) = 0.0_DP
      self%rhill(:) = 0.0_DP

      self%nplpl = 0   

      if (param%lclose) then
         if (allocated(self%radius)) deallocate(self%radius)
         if (allocated(self%density)) deallocate(self%density)
         allocate(self%radius(n))
         allocate(self%density(n))
         self%radius(:) = 0.0_DP
         self%density(:) = 1.0_DP
      end if

      if (param%lrotation) then
         if (allocated(self%rot)) deallocate(self%rot)
         if (allocated(self%Ip)) deallocate(self%Ip)
         allocate(self%rot(NDIM, n))
         allocate(self%Ip(NDIM, n))
         self%rot(:,:) = 0.0_DP
         self%Ip(:,:) = 0.0_DP
      end if

      if (param%ltides) then
         if (allocated(self%k2)) deallocate(self%k2)
         if (allocated(self%Q)) deallocate(self%Q)
         if (allocated(self%tlag)) deallocate(self%tlag)
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

      if (allocated(self%info)) deallocate(self%info)
      if (allocated(self%isperi)) deallocate(self%isperi)
      if (allocated(self%peri)) deallocate(self%peri)
      if (allocated(self%atp)) deallocate(self%atp)

      allocate(swiftest_particle_info :: self%info(n))
      allocate(self%isperi(n))
      allocate(self%peri(n))
      allocate(self%atp(n))

      self%info(:)%particle_type = TP_TYPE_NAME
      self%isperi(:) = 0
      self%peri(:)   = 0.0_DP
      self%atp(:)    = 0.0_DP

      return
   end subroutine setup_tp

end submodule s_setup
