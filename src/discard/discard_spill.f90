submodule (swiftest_classes) s_discard_spill
!! This submodule contains the methods spill, spill_body, spill_pl, and spill_tp for basic swiftest particles
contains
   module procedure discard_spill_body
      !! author: David A. Minton
      !!
      !! Move spilled (discarded) Swiftest particle structure from active list to discard list
      use swiftest    
      implicit none
      integer(I4B) :: i

      associate(np => self%nbody)
         !call symba_spill_body(self, discard)

         ! Pack the discarded bodies into the discard object and do the inverse with the original
         discards%name(:)   = pack(self%name(1:np),         self%ldiscard(1:np))
         self%name(:)       = pack(self%name(1:np),   .not. self%ldiscard(1:np))
         discards%status(:) = pack(self%status(1:np), self%ldiscard(1:np))
         do i = 1, NDIM
            discards%xh(i,:) = pack(self%xh(i,1:np), self%ldiscard(1:np))
            discards%vh(i,:) = pack(self%vh(i,1:np), self%ldiscard(1:np))
            discards%xb(i,:) = pack(self%xb(i,1:np), self%ldiscard(1:np))
            discards%vb(i,:) = pack(self%vb(i,1:np), self%ldiscard(1:np))
         end do
         discards%mu_vec(:) = pack(self%mu_vec(1:np), self%ldiscard(1:np))
         discards%dt_vec(:) = pack(self%dt_vec(1:np), self%ldiscard(1:np))
         discards%a(:) = pack(self%a(1:np), self%ldiscard(1:np))
         discards%e(:) = pack(self%e(1:np), self%ldiscard(1:np))
         discards%inc(:) = pack(self%inc(1:np), self%ldiscard(1:np))
         discards%capom(:) = pack(self%capom(1:np), self%ldiscard(1:np))
         discards%omega(:) = pack(self%omega(1:np), self%ldiscard(1:np))
         discards%capm(:) = pack(self%capm(1:np), self%ldiscard(1:np))
      
         ! Pack the kept bodies back into the original object
         self%name(:)   = pack(self%name(1:np),   .not. self%ldiscard(1:np))
         self%status(:) = pack(self%status(1:np), .not. self%ldiscard(1:np))
         self%mu_vec(:) = pack(self%mu_vec(1:np), .not. self%ldiscard(1:np))
         self%dt_vec(:) = pack(self%dt_vec(1:np), .not. self%ldiscard(1:np))
      end associate
      ! This is the base class, so will be the last to be called. Therefore we must reset spill flag for the next discard operation.
      self%nbody = self%nbody - discards%nbody
      
   end procedure discard_spill_body

   !module procedure discard_spill_body
      !! author: David A. Minton
      !!
      !! Move spilled (discarded) Swiftest particle structure from active list to discard list
      !use swiftest    
      !implicit none
   !
   !   associate(np => self%nbody)
   !
   !
   !   end associate
   !   return 
!   end procedure discard_spill_body

   module procedure discard_spill_pl
      !! author: David A. Minton
      !!
      !! Move spilled (discarded) Swiftest massive body particle structure from active list to discard list
      use swiftest
      implicit none
      
      associate(npl => self%nbody, nspill => self%nspill)
         if (.not. self%lspill) then
            call discard%alloc(nspill) ! Create the discard object for this type
            self%lspill = .true.
         end if

         ! Pack the discarded bodies into the discard object
         discard%mass(:)   = pack(self%mass(1:npl),   self%ldiscard(1:npl))
         discard%radius(:) = pack(self%radius(1:npl), self%ldiscard(1:npl))
         discard%rhill(:)  = pack(self%rhill(1:npl),  self%ldiscard(1:npl))

         ! Pack the kept bodies back into the original object
         self%mass(:)   = pack(self%mass(1:npl),   .not. self%ldiscard(1:npl))
         self%radius(:) = pack(self%radius(1:npl), .not. self%ldiscard(1:npl))
         self%rhill(:)  = pack(self%rhill(1:npl),  .not. self%ldiscard(1:npl))

         ! Call the spill method for the parent class 
         call discard_spill_tp(self,discard)
      end associate
      
      return
   end procedure discard_spill_pl

   module procedure discard_spill_tp
      !! author: The Purdue Swiftest Team -  David A. Minton, Carlisle A. Wishard, Jennifer L.L. Pouplin, and Jacob R. Elliott
      !!
      !! Move spilled (discarded) Swiftest test particle structure from active list to discard list
      use swiftest
      implicit none

      integer(I4B) :: i

      associate(ntp => self%nbody, nspill => self%nspill)
         if (.not. self%lspill) then
            call discard%alloc(nspill) ! Create the discard object for this type
            self%lspill = .true.
         end if

         ! Pack the discarded bodies into the discard object
         discard%peri(:)   = pack(self%peri(1:ntp),   self%ldiscard(1:ntp))
         discard%atp(:)    = pack(self%atp(1:ntp),    self%ldiscard(1:ntp))
         discard%isperi(:) = pack(self%isperi(1:ntp), self%ldiscard(1:ntp))

         ! Pack the kept bodies back into the original object
         self%peri(:)   = pack(self%peri(1:ntp),   .not. self%ldiscard(1:ntp))
         self%atp(:)    = pack(self%atp(1:ntp),    .not. self%ldiscard(1:ntp))
         self%isperi(:) = pack(self%isperi(1:ntp), .not. self%ldiscard(1:ntp))

         do concurrent (i = 1:NDIM)
            discard%xh(i,:) = pack(self%xh(i,1:ntp), self%ldiscard(1:ntp))
            discard%vh(i,:) = pack(self%vh(i,1:ntp), self%ldiscard(1:ntp))
            discard%xb(i,:) = pack(self%xb(i,1:ntp), self%ldiscard(1:ntp))
            discard%vb(i,:) = pack(self%vb(i,1:ntp), self%ldiscard(1:ntp))
            self%xh(i,:)    = pack(self%xh(i,1:ntp), .not. self%ldiscard(1:ntp))
            self%vh(i,:)    = pack(self%vh(i,1:ntp), .not. self%ldiscard(1:ntp))
            self%xb(i,:)    = pack(self%xb(i,1:ntp), .not. self%ldiscard(1:ntp))
            self%vb(i,:)    = pack(self%vb(i,1:ntp), .not. self%ldiscard(1:ntp))
         end do

      end associate
      return
   end procedure discard_spill_tp
end submodule s_discard_spill



      


