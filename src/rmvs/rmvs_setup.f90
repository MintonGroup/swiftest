submodule(rmvs_classes) s_rmvs_setup
contains
   module procedure rmvs_setup_pl
      !! author: David A. Minton
      !!
      !! Allocate RMVS test particle structure
      !!
      !! Equivalent in functionality to David E. Kaufmann's Swifter routine rmvs_setup.f90
      use swiftest
      implicit none

      !> Call allocation method for parent class
      call whm_setup_pl(self, n) 
      if (n <= 0) return

      allocate(self%nenc(n))
      allocate(self%tpenc1P(n))
      allocate(self%xout(NDIM, n, 0:NTENC))
      allocate(self%vout(NDIM, n, 0:NTENC))
      allocate(self%xin(NDIM, n, 0:NTPHENC))
      allocate(self%vin(NDIM, n, 0:NTPHENC))
      allocate(self%xpc(NDIM, n, 0:NTPHENC))
      allocate(self%aoblin(NDIM, n, 0:NTPHENC))

      self%nenc       = 0
      self%tpenc1P(:) = 0
      self%xout(:,:,:)  = 0.0_DP
      self%vout(:,:,:)  = 0.0_DP
      self%xin(:,:,:)  = 0.0_DP
      self%vin(:,:,:)  = 0.0_DP
      self%xpc(:,:,:)   = 0.0_DP
      self%aoblin(:,:,:)   = 0.0_DP

      return
   end procedure rmvs_setup_pl 

   module procedure rmvs_setup_tp
      !! author: David A. Minton
      !!
      !! Allocate WHM test particle structure
      !!
      !! Equivalent in functionality to David E. Kaufmann's Swifter routine whm_setup.f90
      use swiftest
      implicit none

      !> Call allocation method for parent class
      call whm_setup_tp(self, n) 
      if (n <= 0) return

      allocate(self%lperi(n))
      allocate(self%xpc(NDIM, n))
      allocate(self%vpc(NDIM, n))
      allocate(self%apc(NDIM, n))
      allocate(self%plperP(n))
      allocate(self%plencP(n))
      allocate(self%tpencP(n))

      self%lperi(:)  = .false.
      self%xpc(:,:)  = 0.0_DP
      self%vpc(:,:)  = 0.0_DP
      self%apc(:,:)  = 0.0_DP
      self%plperP(:) = 0
      self%plencP(:) = 0
      self%tpencP(:) = 0

      return
   end procedure rmvs_setup_tp

   module procedure rmvs_setup_system
      !! author: David A. Minton
      !!
      !! Wrapper method to initialize a basic Swiftest nbody system from files
      !!
      implicit none

      ! Call parent method
      call whm_setup_system(self, config)

   end procedure rmvs_setup_system


   module procedure rmvs_setup_encounter
      !! author: David A. Minton
      !!
      !! When encounters are detected, this method will call the interpolation methods for the planets and 
      !! creates a Swiftest test particle structure for each planet's encountering test particles to simplify the 
      !! planetocentric calculations. This subroutine is not based on an existing one from Swift and Swifter
      !!
      implicit none
      integer(I4B) :: i, j, k, link, nenc

      associate(npl => self%nbody)
         allocate(self%tpenc(npl))
         allocate(self%plenc(npl, 0:NTPHENC))
         allocate(self%cbenc(npl))
         do i = 1, npl
            nenc = self%nenc(i) 
            if (nenc > 0) then
               ! There are inner encounters with this planet...first make the planet a central body
               self%cbenc(i)%Gmass = self%Gmass(i)

               ! Next create an encountering test particle structure
               call self%tpenc(i)%setup(nenc)  
               link = self%tpenc1P(i)
               do j = 1, nenc
                  self%tpenc(i)%name(j) = tp%name(link)
                  self%tpenc(i)%status(j) = tp%status(link)
                  link = tp%tpencP(link)
               end do

               call self%tpenc(i)%set_mu(self%cbenc(i)) ! Make sure that the test particles get the proper value of mu 
               
               ! Now create a planetocentric "planet" structure containing the *other* planets (plus the Sun) in it
               do j = 0, NTPHENC
                  call self%plenc(i, j)%setup(npl)
                  self%plenc(i, j)%status(:) = self%status(:)
               end do
               do k = 0, NTPHENC
                  do j = 1, npl
                     if (j == i) then ! We will substitute the Sun in the array location occupied by the encountering planet
                        self%plenc(i, k)%name(j) = 0
                        self%plenc(i, k)%Gmass(j) = cb%Gmass
                        self%plenc(i, k)%mass(j) = cb%mass
                        self%plenc(i, k)%radius(j) = cb%radius
                     else
                        self%plenc(i, k)%name(j) = self%name(i)
                        self%plenc(i, k)%Gmass(j) = self%Gmass(i)
                        self%plenc(i, k)%mass(j) = self%mass(i)
                     end if
                  end do
               end do
               
            end if

         end do
      end associate


   end procedure rmvs_setup_encounter

   module procedure rmvs_destruct_encounter
      !! author: David A. Minton
      !!
      !! Deallocates all of the encountering particle data structures for next time
      !!
      use swiftest
      implicit none


      deallocate(self%tpenc)
      deallocate(self%cbenc)
      deallocate(self%plenc)

      return

   end procedure rmvs_destruct_encounter


   module procedure rmvs_setup_set_beg_end
      !! author: David A. Minton
      !! 
      !! Sets one or more of the values of xbeg, xend, and vbeg
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
      if (present(vbeg)) then
         if (allocated(self%vbeg)) deallocate(self%vbeg)
         allocate(self%vbeg, source=vbeg)
      end if

      return

   end procedure rmvs_setup_set_beg_end


end submodule s_rmvs_setup
