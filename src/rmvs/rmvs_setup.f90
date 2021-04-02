submodule(rmvs_classes) s_rmvs_setup
   use swiftest
contains
   module subroutine rmvs_setup_pl(self,n)
      !! author: David A. Minton
      !!
      !! Allocate RMVS test particle structure
      !!
      !! Equivalent in functionality to David E. Kaufmann's Swifter routine rmvs_setup.f90
      implicit none
      ! Arguments
      class(rmvs_pl),                intent(inout) :: self !! RMVS test particle object
      integer(I4B),                  intent(in)    :: n    !! Number of test particles to allocate
      ! Internals
      integer(I4B)                                 :: i,j

      !> Call allocation method for parent class
      call whm_setup_pl(self, n) 
      if (n <= 0) return

      allocate(self%nenc(n))
      allocate(self%xout(NDIM, n, 0:NTENC))
      allocate(self%vout(NDIM, n, 0:NTENC))
      allocate(self%xin(NDIM, n, 0:NTPHENC))
      allocate(self%vin(NDIM, n, 0:NTPHENC))
      allocate(self%aoblin(NDIM, n, 0:NTPHENC))
      allocate(self%plind(n,n))
      allocate(self%tpenc(n))
      allocate(self%plenc(n))
      allocate(self%cbenc(n))
      self%plenc(:)%nbody = n

      self%nenc          = 0
      self%xout(:,:,:)   = 0.0_DP
      self%vout(:,:,:)   = 0.0_DP
      self%xin(:,:,:)    = 0.0_DP
      self%vin(:,:,:)    = 0.0_DP
      self%aoblin(:,:,:) = 0.0_DP

      do j = 1, n
         self%plind(j,:) = [(i,i=1,n)] 
         self%plind(j,2:n) = pack(self%plind(j,1:n), self%plind(j,1:n) /= j)
      end do

      return
   end subroutine rmvs_setup_pl 

   module subroutine rmvs_setup_tp(self,n)
      !! author: David A. Minton
      !!
      !! Allocate WHM test particle structure
      !!
      !! Equivalent in functionality to David E. Kaufmann's Swifter routine whm_setup.f90
      implicit none
      ! Arguments
      class(rmvs_tp),              intent(inout)   :: self !! RMVS test particle object
      integer,                     intent(in)      :: n    !! Number of test particles to allocate

      !> Call allocation method for parent class
      call whm_setup_tp(self, n) 
      if (n <= 0) return

      allocate(self%lperi(n))
      allocate(self%plperP(n))
      allocate(self%plencP(n))

      self%lperi(:)  = .false.
      self%plperP(:) = 0
      self%plencP(:) = 0

      return
   end subroutine rmvs_setup_tp

   module subroutine rmvs_setup_system(self, config)
      !! author: David A. Minton
      !!
      !! Wrapper method to initialize a basic Swiftest nbody system from files
      !!
      implicit none
      ! Arguments
      class(rmvs_nbody_system),      intent(inout) :: self    !! RMVS system object
      class(swiftest_configuration), intent(inout) :: config  !! Input collection of  configuration parameters 
      ! Internals
      integer(I4B) :: i
      ! Call parent method
      call whm_setup_system(self, config)

      ! Set up the tp-planet encounter structures
      select type(pl => self%pl)
      class is(rmvs_pl)
         select type(cb => self%cb)
         class is (rmvs_cb)
            select type (tp => self%tp)
            class is (rmvs_tp)
               associate(npl => pl%nbody)
                  tp%cb = cb
                  do i = 1, npl
                     allocate(pl%plenc(i)%Gmass(npl))
                     allocate(pl%plenc(i)%xh(NDIM,npl))
                     allocate(pl%plenc(i)%vh(NDIM,npl))
                     pl%cbenc(i)          = cb
                     pl%cbenc(i)%Gmass    = pl%Gmass(i)
                     pl%plenc(i)%Gmass(1) = cb%Gmass
                  end do
               end associate
            end select
         end select
      end select

   end subroutine rmvs_setup_system
   
   module subroutine rmvs_setup_set_beg_end(self, xbeg, xend, vbeg)
      !! author: David A. Minton
      !! 
      !! Sets one or more of the values of xbeg, xend, and vbeg
      implicit none
      ! Arguments
      class(rmvs_tp),                intent(inout) :: self !! RMVS test particle object
      real(DP), dimension(:,:),           optional :: xbeg, xend, vbeg

      if (present(xbeg)) then
         if (allocated(self%xbeg)) deallocate(self%xbeg)
         allocate(self%xbeg, source=xbeg)
      end if
      if (present(xend)) then
         if (allocated(self%xend)) deallocate(self%xend)
         allocate(self%xend, source=xend)
      end if
      if (present(vbeg)) then
         if (allocated(self%vbeg)) deallocate(self%vbeg)
         allocate(self%vbeg, source=vbeg)
      end if

      return

   end subroutine rmvs_setup_set_beg_end

end submodule s_rmvs_setup