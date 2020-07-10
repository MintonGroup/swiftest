submodule (swiftest_classes) setup_implementations

contains
   module procedure setup_construct_system
      !! author: David A. Minton
      !!
      !! Constructor for a Swiftest nbody system. Creates the nbody system object based on the user-input integrator
      !! 
      use swiftest
      implicit none

      select case(integrator)
      case (BS)
         write(*,*) 'Bulirsch-Stoer integrator not yet enabled'
      case (HELIO)
         write(*,*) 'Democratic Heliocentric integrator not yet enabled'
      case (RA15)
         write(*,*) 'Radau integrator not yet enabled'
      case (TU4)
         write(*,*) 'TU4 integrator not yet enabled'
      case (WHM)
         allocate(whm_nbody_system :: system)
         allocate(whm_configuration :: config)
         select type(system)
         class is (whm_nbody_system)
            allocate(whm_cb :: system%cb)
            allocate(whm_pl :: system%pl)
            allocate(whm_tp :: system%tp)
            allocate(whm_tp :: system%tp_discards)
         end select
      case (RMVS)
         write(*,*) 'RMVS integrator not yet enabled'
      case (SYMBA)
         write(*,*) 'SyMBA integrator not yet enabled'
      case (RINGMOONS)
         write(*,*) 'RINGMOONS-SyMBA integrator not yet enabled'
      case default
         write(*,*) 'Unkown integrator',integrator
         call util_exit(FAILURE)
      end select

      return
   end procedure setup_construct_system

   module procedure setup_body
      !! author: David A. Minton
      !!
      !! Constructor for base Swiftest particle class. Allocates space for all particles and
      !! initializes all components with a value.
      !! Note: Timing tests indicate that (NDIM, n) is more efficient than (NDIM, n) 
      use swiftest
      implicit none

      self%nbody = n
      if (n <= 0) return

      !write(*,*) 'Allocating the basic Swiftest particle'
      allocate(self%name(n))
      allocate(self%status(n))
      allocate(self%ldiscard(n))
      allocate(self%xh(NDIM, n))
      allocate(self%vh(NDIM, n))
      allocate(self%xb(NDIM, n))
      allocate(self%vb(NDIM, n))
      allocate(self%ah(NDIM, n))
      allocate(self%aobl(NDIM, n))
      allocate(self%a(n))
      allocate(self%e(n))
      allocate(self%inc(n))
      allocate(self%capom(n))
      allocate(self%omega(n))
      allocate(self%capm(n))
      allocate(self%mu(n))

      self%name(:)   = 0
      self%status(:) = INACTIVE
      self%ldiscard(:) = .false.
      self%xh(:,:)   = 0.0_DP
      self%vh(:,:)   = 0.0_DP
      self%xb(:,:)   = 0.0_DP
      self%vb(:,:)   = 0.0_DP
      self%ah(:,:)   = 0.0_DP
      self%aobl(:,:)  = 0.0_DP
      self%a(:)      = 0.0_DP
      self%e(:)      = 0.0_DP
      self%inc(:)    = 0.0_DP
      self%capom(:)  = 0.0_DP
      self%omega(:)  = 0.0_DP
      self%capm(:)   = 0.0_DP
      self%a(:)      = 0.0_DP
      self%mu(:)     = 0.0_DP

      return
   end procedure setup_body

   module procedure setup_pl
      !! author: David A. Minton
      !!
      !! Constructor for base Swiftest massive body class. Allocates space for all particles and
      !! initializes all components with a value. 
      use swiftest
      implicit none

      !> Call allocation method for parent class
      !> The parent class here is the abstract swiftest_body class, so we can't use the type-bound procedure
      call setup_body(self, n)
      if (n <= 0) return 

      allocate(self%mass(n))
      allocate(self%Gmass(n))
      allocate(self%rhill(n))
      allocate(self%radius(n))
      allocate(self%density(n))
      allocate(self%Ip(NDIM, n))
      allocate(self%rot(NDIM, n))
      allocate(self%k2(n))
      allocate(self%Q(n))

      self%mass(:) = 0.0_DP
      self%Gmass(:) = 0.0_DP
      self%rhill(:) = 0.0_DP
      self%radius(:) = 0.0_DP
      self%density(:) = 0.0_DP
      self%Ip(:,:) = 0.0_DP
      return
   end procedure setup_pl

   module procedure setup_tp
      !! author: David A. Minton
      !!
      !! Constructor for base Swiftest test particle particle class. Allocates space for 
      !! all particles and initializes all components with a value. 
      use swiftest
      implicit none

      !> Call allocation method for parent class
      !> The parent class here is the abstract swiftest_body class, so we can't use the type-bound procedure
      call setup_body(self, n)
      if (n <= 0) return

      allocate(self%isperi(n))
      allocate(self%peri(n))
      allocate(self%atp(n))

      self%isperi(:) = 0
      self%peri(:)   = 0.0_DP
      self%atp(:)    = 0.0_DP

      return
   end procedure setup_tp

   module procedure setup_set_msys
      !! author: David A. Minton
      !!
      !! Sets the value of msys and the vector mass quantities based on the total mass of the system
      self%msys = self%cb%mass + sum(self%pl%mass(1:self%pl%nbody))

      return
   end procedure setup_set_msys

   module procedure setup_set_mu_pl
      !! author: David A. Minton
      !!
      !! Computes G * (M + m) for each massive body
      use swiftest
      implicit none

      if (self%nbody > 0) self%mu(:) = cb%Gmass + self%Gmass(:)

      return
   end procedure setup_set_mu_pl

   module procedure setup_set_mu_tp
      !! author: David A. Minton
      !!
      !! Converts certain scalar values to arrays so that they can be used in elemental functions
      use swiftest
      implicit none

      if (self%nbody > 0) self%mu(:) = cb%Gmass

      return
   end procedure setup_set_mu_tp

end submodule setup_implementations
