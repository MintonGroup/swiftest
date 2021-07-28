submodule(symba_classes) s_symba_setup
   use swiftest
contains
   module subroutine symba_setup_pl(self, n)
      !! author: David A. Minton
      !!
      !! Allocate SyMBA test particle structure
      !!
      !! Equivalent in functionality to David E. Kaufmann's Swifter routine symba_setup.f90
      implicit none
      ! Arguments
      class(symba_pl), intent(inout) :: self !! SyMBA test particle object
      integer(I4B),    intent(in)    :: n    !! Number of massive bodies to allocate
      ! Internals
      integer(I4B)                   :: i

      !> Call allocation method for parent class
      call setup_pl(self, n) 
      if (n <= 0) return
      allocate(self%lcollision(n))
      allocate(self%lencounter(n))
      allocate(self%nplenc(n))
      allocate(self%ntpenc(n))
      allocate(self%levelg(n))
      allocate(self%levelm(n))
      allocate(self%isperi(n))
      allocate(self%peri(n))
      allocate(self%atp(n))
      allocate(self%kin(n))
      allocate(self%info(n))

      self%lcollision(:) = .false.
      self%lencounter(:) = .false.
      self%nplenc(:) = 0
      self%ntpenc(:) = 0
      self%levelg(:) = -1
      self%levelm(:) = -1
      self%isperi(:) = 0
      self%peri(:) = 0.0_DP
      self%atp(:) = 0.0_DP
      self%kin(:)%nchild = 0
      self%kin(:)%parent = [(i, i=1, n)]
      return
   end subroutine symba_setup_pl

   module subroutine symba_setup_pltpenc(self, n)
      !! author: David A. Minton
      !!
      !! A constructor that sets the number of encounters and allocates and initializes all arrays  
      !!
      implicit none
      ! Arguments
      class(symba_pltpenc), intent(inout) :: self !! Symba pl-tp encounter structure
      integer,              intent(in)    :: n    !! Number of encounters to allocate space for

      self%nenc = n
      if (n == 0) return
      if (allocated(self%lvdotr)) deallocate(self%lvdotr)
      if (allocated(self%status)) deallocate(self%status)
      if (allocated(self%level)) deallocate(self%level)
      if (allocated(self%index1)) deallocate(self%index1)
      if (allocated(self%index2)) deallocate(self%index2)
      allocate(self%lvdotr(n))
      allocate(self%status(n))
      allocate(self%level(n))
      allocate(self%index1(n))
      allocate(self%index2(n))
      self%lvdotr(:) = .false.
      self%status(:) = INACTIVE
      self%level(:) = -1
      self%index1(:) = 0
      self%index2(:) = 0
      return
   end subroutine symba_setup_pltpenc

   module subroutine symba_setup_plplenc(self,n)
      !! author: David A. Minton
      !!
      !! A constructor that sets the number of encounters and allocates and initializes all arrays  
      !
      implicit none
      ! Arguments
      class(symba_plplenc), intent(inout) :: self !! Symba pl-tp encounter structure
      integer,              intent(in)    :: n    !! Number of encounters to allocate space for

      call symba_setup_pltpenc(self, n)
      if (n == 0) return
      if (allocated(self%xh1)) deallocate(self%xh1)
      if (allocated(self%xh2)) deallocate(self%xh2)
      if (allocated(self%vb1)) deallocate(self%vb1)
      if (allocated(self%vb2)) deallocate(self%vb2)
      allocate(self%xh1(NDIM,n))
      allocate(self%xh2(NDIM,n))
      allocate(self%vb1(NDIM,n))
      allocate(self%vb2(NDIM,n))
      self%xh1(:,:) = 0.0_DP
      self%xh2(:,:) = 0.0_DP
      self%vb1(:,:) = 0.0_DP
      self%vb2(:,:) = 0.0_DP
      return
   end subroutine symba_setup_plplenc

   module subroutine symba_setup_initialize_system(self, param)
      !! author: David A. Minton
      !!
      !! Initialize an SyMBA nbody system from files and sets up the planetocentric structures.
      !! This subroutine will also sort the massive bodies in descending order by mass
      !! 
      implicit none
      ! Arguments
      class(symba_nbody_system),  intent(inout) :: self    !! SyMBA system object
      class(swiftest_parameters), intent(inout) :: param  !! Current run configuration parameters 
      ! Internals
      integer(I4B) :: i, j

      ! Call parent method
      associate(system => self)
         call whm_setup_initialize_system(system, param)
         call system%mergeadd_list%setup(1)
         call system%mergesub_list%setup(1)
         call system%pltpenc_list%setup(1)
         call system%plplenc_list%setup(1)
         select type(pl => system%pl)
         class is (symba_pl)
            call pl%sort("mass", ascending=.false.)
            select type(param)
            class is (symba_parameters)
               pl%lmtiny(:) = pl%Gmass(:) > param%MTINY
               pl%nplm = count(pl%lmtiny(:))
            end select
         end select
      end associate
      return
   end subroutine symba_setup_initialize_system

   module subroutine symba_setup_tp(self,n)
      !! author: David A. Minton
      !!
      !! Allocate WHM test particle structure
      !!
      !! Equivalent in functionality to David E. Kaufmann's Swifter routine whm_setup.f90
      implicit none
      ! Arguments
      class(symba_tp), intent(inout) :: self !! SyMBA test particle object
      integer,         intent(in)    :: n    !! Number of test particles to allocate

      !> Call allocation method for parent class
      call setup_tp(self, n) 
      if (n <= 0) return
      allocate(self%nplenc(n))
      allocate(self%levelg(n))
      allocate(self%levelm(n))
      self%nplenc(:) = 0
      self%levelg(:) = -1
      self%levelm(:) = -1
      return
   end subroutine symba_setup_tp

end submodule s_symba_setup
