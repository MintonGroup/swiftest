submodule(helio_classes) s_helio_setup
contains
   module procedure helio_setup_pl
      !! author: David A. Minton & Carlisle A. Wishard
      !!
      !! Allocate Helio planet structure
      !!
      !! Equivalent in functionality to David E. Kaufmann's Swifter routine helio_setup.f90
      use swiftest
      implicit none

      !> Call allocation method for parent class
      call setup_pl(self, n) 
      if (n <= 0) return

      allocate(self%eta(n))
      allocate(self%muj(n))
      allocate(self%xj(NDIM, n))
      allocate(self%vj(NDIM, n))
      allocate(self%ir3j(n))

      self%eta(:)   = 0.0_DP
      self%muj(:)   = 0.0_DP
      self%xj(:,:)  = 0.0_DP
      self%vj(:,:)  = 0.0_DP
      self%ir3j(:) = 0.0_DP

      return
   end procedure helio_setup_pl 

   module procedure helio_setup_tp
      !! author: David A. Minton & Carlisle A. Wishard
      !!
      !! Allocate Helio test particle structure
      !!
      !! Equivalent in functionality to David E. Kaufmann's Swifter routine helio_setup.f90
      use swiftest
      implicit none

      !> Call allocation method for parent class
      call setup_tp(self, n) 
      if (n <= 0) return

      return
   end procedure helio_setup_tp

   module procedure helio_setup_system
      !! author: David A. Minton & Carlisle A. Wishard
      !!
      !! Wrapper method to initialize a basic Swiftest nbody system from files
      !!
      use swiftest
      implicit none

      call io_read_initialize_system(self, config)
      ! Make sure that the discard list gets allocated initially
      call self%tp_discards%setup(self%tp%nbody)

      if (self%pl%nbody > 0) then
         select type(pl => self%pl)
         class is (helio_pl)
            call pl%set_mu(self%cb)
            if (config%lgr) call pl%gr_vh2pv(config)
            call pl%eucl_index()
         end select
      end if

      if (self%tp%nbody > 0) then
         select type(tp => self%tp)
         class is (helio_tp)
            call tp%set_mu(self%cb)
            if (config%lgr) call tp%gr_vh2pv(config)
         end select
      end if

   end procedure helio_setup_system

   module procedure helio_setup_set_beg_end
      !! author: David A. Minton & Carlisle A. Wishard
      !! 
      !! Sets one or more of the values of xbeg and xend
      use swiftest
      implicit none

      if (present(xbeg)) then
         if (allocated(self%xbeg)) deallocate(self%xbeg)
         allocate(self%xbeg, source=xbeg)
      end if
      if (present(xend)) then
         if (allocated(self%xend)) deallocate(self%xend)
         allocate(self%xend, source=xend)
      end if

      return

   end procedure helio_setup_set_beg_end

end submodule s_helio_setup