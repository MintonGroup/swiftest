submodule (swiftest_classes) s_nbody_allocate

contains

   module procedure nbody_allocate_body
      !! author: David A. Minton
      !!
      !! Constructor for base Swiftest particle class. Allocates space for all particles and
      !! initializes all components with a value. Also sets the is_allocated flag to true.
      use swiftest
      implicit none

      self%nbody = n
      if (n <= 0) return

      !write(*,*) 'Allocating the basic Swiftest particle'
      allocate(self%name(n))
      allocate(self%status(n))
      allocate(self%xh(NDIM,n))
      allocate(self%vh(NDM,n))
      allocate(self%xb(NDIM,n))
      allocate(self%vb(NDIM,n))
      allocate(self%mu_vec(n))
      allocate(self%dt_vec(n))
      allocate(self%a(n))
      allocate(self%e(n))
      allocate(self%inc(n))
      allocate(self%capom(n))
      allocate(self%omega(n))
      allocate(self%capm(n))

      self%name(:)   = 0
      self%status(:) = INACTIVE
      self%xh(:,:)   = 0.0_DP
      self%vh(:,:)   = 0.0_DP
      self%xb(:,:)   = 0.0_DP
      self%vb(:,:)   = 0.0_DP
      self%a(:)      = 0.0_DP
      self%e(:)      = 0.0_DP
      self%inc(:)    = 0.0_DP
      self%capom(:)  = 0.0_DP
      self%omega(:)  = 0.0_DP
      self%capm(:)   = 0.0_DP
      self%a(:)      = 0.0_DP
      self%mu_vec(:) = 0.0_DP
      self%dt_vec(:) = 0.0_DP

      return
   end procedure nbody_allocate_body

   module procedure nbody_allocate_pl
      !! author: David A. Minton
      !!
      !! Constructor for base Swiftest massive body class. Allocates space for all particles and
      !! initializes all components with a value. 
      use swiftest
      implicit none

      if (self%is_allocated) then
         !write(*,*) 'Swiftest massive body structure already alllocated'
         return
      end if

      !> Call allocation method for parent class
      call self%swiftest_body%alloc(n)
      if (n <= 0) return 

      allocate(self%mass(n))
      allocate(self%Gmass(n))
      allocate(self%rhill(n))
      allocate(self%radius(n))
      allocate(self%density(n))
      allocate(self%Ip(NDIM,n))
      allocate(self%rot(NDIM,n))
      allocate(self%k2(n))
      allocate(self%Q(n))

      self%mass(:) = 0.0_DP
      self%Gmass(:) = 0.0_DP
      self%rhill(:) = 0.0_DP
      self%radius(:) = 0.0_DP
      self%density(:) = 0.0_DP
      self%Ip(:,):) = 0.0_DP
      return
   end procedure nbody_allocate_pl

   module procedure nbody_allocate_tp
      !! author: David A. Minton
      !!
      !! Constructor for base Swiftest test particle particle class. Allocates space for 
      !! all particles and initializes all components with a value. 
      use swiftest
      implicit none

      !> Call allocation method for parent class
      call self%swiftest_body%alloc(n)
      if (n <= 0) return

      return
   end procedure nbody_allocate_tp

   ! module procedure nbody_deallocate_body
   !    !! author: David A. Minton
   !    !!
   !    !! Finalizer for base Swiftest particle class. Deallocates all components and sets 
   !    !! is_allocated flag to false. Mostly this is redundant, so this serves as a placeholder
   !    !! in case future updates include pointers as part of the class.
   !    use swiftest
   !    implicit none
      
   !    if (self%is_allocated) then
   !       deallocate(self%name)
   !       deallocate(self%status)
   !       if (allocated(self%lspill_list)) deallocate(self%lspill_list)
   !       self%is_allocated = .false.
   !    end if
   !    return
   ! end procedure nbody_deallocate_body

   ! module procedure nbody_allocate_tp
   !    !! author: David A. Minton
   !    !!
   !    !! Constructor for base Swiftest test particle particle class. Allocates space for 
   !    !! all particles and initializes all components with a value. 
   !    use swiftest
   !    implicit none

   !    if (self%is_allocated) then
   !       !write(*,*) 'Swiftest test particle structure already alllocated'
   !       return
   !    end if
   !    !write(*,*) 'Allocating the Swiftest test particle'

   !    !> Call allocation method for parent class
   !    call self%swiftest_body%alloc(n)
   !    if (n <= 0) return

   !    allocate(self%peri(n))
   !    allocate(self%atp(n))
   !    allocate(self%isperi(n))
   !    allocate(self%xh(NDIM,n))
   !    allocate(self%vh(NDIM,n))
   !    allocate(self%xb(NDIM,n))
   !    allocate(self%vb(NDIM,n))

   !    self%peri(:) = 0.0_DP
   !    self%atp(:) = 0.0_DP
   !    self%isperi(:) = 0.0_DP
   !    self%xh(:,:) = 0.0_DP
   !    self%vh(:,:) = 0.0_DP
   !    self%xb(:,:) = 0.0_DP
   !    self%vb(:,:) = 0.0_DP

   !    return
   ! end procedure nbody_allocate_tp

   ! module procedure nbody_deallocate_tp
   !    !! author: David A. Minton
   !    !!
   !    !! Finalizer for base Swiftest particle class.
   !    !! Basic Swiftest test particle destructor/finalizer
   !    use swiftest
   !    implicit none

   !    if (self%is_allocated) then
   !       deallocate(self%isperi)
   !       deallocate(self%peri)
   !       deallocate(self%atp)
   !       deallocate(self%xh)
   !       deallocate(self%vh)
   !       deallocate(self%xb)
   !       deallocate(self%vb)
   !    end if
   !    return
   ! end procedure nbody_deallocate_tp

end submodule s_nbody_allocate
