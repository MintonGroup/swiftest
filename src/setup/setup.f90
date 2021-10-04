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


   module subroutine setup_finalize_system(self, param)
      !! author: David A. Minton
      !!
      !! Runs any finalization subroutines when ending the simulation.
      !!
      implicit none
      ! Arguments
      class(swiftest_nbody_system), intent(inout) :: self   !! Swiftest system object
      class(swiftest_parameters),   intent(inout) :: param  !! Current run configuration parameters
      integer(I4B)                                :: ierr

      associate(system => self)
         if ((param%out_type == NETCDF_FLOAT_TYPE) .or. (param%out_type == NETCDF_DOUBLE_TYPE)) then
            call param%nciu%close()
         end if
      end associate

      return
   end subroutine setup_finalize_system


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
      integer(I4B) :: i

      associate(cb => self%cb, pl => self%pl, npl => self%pl%nbody, tp => self%tp, ntp => self%tp%nbody)

         call cb%info%set_value(particle_type=CB_TYPE_NAME, status="ACTIVE", origin_type="Initial conditions", origin_time=param%t0, origin_xh=[0.0_DP, 0.0_DP, 0.0_DP], origin_vh=[0.0_DP, 0.0_DP, 0.0_DP], collision_id=0)
         do i = 1, self%pl%nbody
            call pl%info(i)%set_value(particle_type=PL_TYPE_NAME, status="ACTIVE", origin_type="Initial conditions", origin_time=param%t0, origin_xh=self%pl%xh(:,i), origin_vh=self%pl%vh(:,i), collision_id=0)
         end do
         do i = 1, self%tp%nbody
            call tp%info(i)%set_value(particle_type=TP_TYPE_NAME, status="ACTIVE", origin_type="Initial conditions", origin_time=param%t0, origin_xh=self%tp%xh(:,i), origin_vh=self%tp%vh(:,i), collision_id=0)
         end do

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
      integer(I4B)                                :: ierr

      associate(system => self, cb => self%cb, pl => self%pl, tp => self%tp)

         call system%read_in(param)
         call system%validate_ids(param)
         call system%set_msys()
         call pl%set_mu(cb) 
         call tp%set_mu(cb) 
         if (param%in_form == EL) then
            call pl%el2xv(cb)
            call tp%el2xv(cb)
         end if
         call pl%flatten(param)
         if (.not.param%lrhill_present) call pl%set_rhill(cb)
         pl%lfirst = param%lfirstkick
         tp%lfirst = param%lfirstkick

         if (param%lrestart) then
            call system%read_particle_info(param)
         else
            call system%init_particle_info(param)
         end if
      end associate

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

      if (n < 0) return

      self%lfirst = .true.

      if (allocated(self%info)) deallocate(self%info)
      if (allocated(self%id)) deallocate(self%id)
      if (allocated(self%status)) deallocate(self%status)
      if (allocated(self%ldiscard)) deallocate(self%ldiscard)
      if (allocated(self%lmask)) deallocate(self%lmask)
      if (allocated(self%mu)) deallocate(self%mu)
      if (allocated(self%xh)) deallocate(self%xh)
      if (allocated(self%vh)) deallocate(self%vh)
      if (allocated(self%xb)) deallocate(self%xb)
      if (allocated(self%vb)) deallocate(self%vb)
      if (allocated(self%ah)) deallocate(self%ah)
      if (allocated(self%aobl)) deallocate(self%aobl)
      if (allocated(self%agr)) deallocate(self%lmask)
      if (allocated(self%atide)) deallocate(self%lmask)
      if (allocated(self%ir3h)) deallocate(self%ir3h)
      if (allocated(self%a)) deallocate(self%a)
      if (allocated(self%e)) deallocate(self%e)
      if (allocated(self%e)) deallocate(self%e)
      if (allocated(self%inc)) deallocate(self%inc)
      if (allocated(self%capom)) deallocate(self%capom)
      if (allocated(self%omega)) deallocate(self%omega)
      if (allocated(self%capm)) deallocate(self%capm)

      self%nbody = n
      if (n == 0) return

      allocate(self%info(n))
      allocate(self%id(n))
      allocate(self%status(n))
      allocate(self%ldiscard(n))
      allocate(self%lmask(n))
      allocate(self%mu(n))
      allocate(self%xh(NDIM, n))
      allocate(self%vh(NDIM, n))
      allocate(self%xb(NDIM, n))
      allocate(self%vb(NDIM, n))
      allocate(self%ah(NDIM, n))
      allocate(self%ir3h(n))

      self%id(:) = 0
      do i = 1, n
         call self%info(i)%set_value(&
            name = "UNNAMED", &
            particle_type = "UNKNOWN", &
            status = "INACTIVE", & 
            origin_type = "UNKNOWN", &
            collision_id = 0, &
            origin_time = -huge(1.0_DP), & 
            origin_xh = [0.0_DP, 0.0_DP, 0.0_DP], &
            origin_vh = [0.0_DP, 0.0_DP, 0.0_DP], &
            discard_time = -huge(1.0_DP), & 
            discard_xh = [0.0_DP, 0.0_DP, 0.0_DP], &
            discard_vh = [0.0_DP, 0.0_DP, 0.0_DP], &
            discard_body_id = -1  &
         )
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
      if (n < 0) return 

      if (allocated(self%mass)) deallocate(self%mass)
      if (allocated(self%Gmass)) deallocate(self%Gmass)
      if (allocated(self%rhill)) deallocate(self%rhill)
      if (allocated(self%renc)) deallocate(self%renc)
      if (allocated(self%radius)) deallocate(self%radius)
      if (allocated(self%density)) deallocate(self%density)
      if (allocated(self%rot)) deallocate(self%rot)
      if (allocated(self%Ip)) deallocate(self%Ip)
      if (allocated(self%k2)) deallocate(self%k2)
      if (allocated(self%Q)) deallocate(self%Q)
      if (allocated(self%tlag)) deallocate(self%tlag)

      if (n == 0) return

      allocate(self%mass(n))
      allocate(self%Gmass(n))
      allocate(self%rhill(n))
      allocate(self%renc(n))

      self%mass(:) = 0.0_DP
      self%Gmass(:) = 0.0_DP
      self%rhill(:) = 0.0_DP
      self%renc(:) = 0.0_DP

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
      if (n < 0) return

      if (allocated(self%isperi)) deallocate(self%isperi)
      if (allocated(self%peri)) deallocate(self%peri)
      if (allocated(self%atp)) deallocate(self%atp)

      if (n == 0) return

      allocate(self%isperi(n))
      allocate(self%peri(n))
      allocate(self%atp(n))

      self%isperi(:) = 0
      self%peri(:)   = 0.0_DP
      self%atp(:)    = 0.0_DP

      return
   end subroutine setup_tp

end submodule s_setup
