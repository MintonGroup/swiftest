!! Copyright 2022 - David Minton, Carlisle Wishard, Jennifer Pouplin, Jake Elliott, & Dana Singh
!! This file is part of Swiftest.
!! Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
!! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!! Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
!! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!! You should have received a copy of the GNU General Public License along with Swiftest. 
!! If not, see: https://www.gnu.org/licenses. 

submodule(symba) s_symba_util
   use swiftest
   use fraggle
contains

   module subroutine symba_util_append_pl(self, source, lsource_mask)
      !! author: David A. Minton
      !!
      !! Append components from one massive body object to another. 
      !! This method will automatically resize the destination body if it is too small
      implicit none
      !! Arguments
      class(symba_pl),                 intent(inout) :: self         !! SyMBA massive body object
      class(swiftest_body),            intent(in)    :: source       !! Source object to append
      logical, dimension(:),           intent(in)    :: lsource_mask !! Logical mask indicating which elements to append to

      select type(source)
      class is (symba_pl)
         call util_append(self%levelg, source%levelg, lsource_mask=lsource_mask)
         call util_append(self%levelm, source%levelm, lsource_mask=lsource_mask)

         call swiftest_util_append_pl(self, source, lsource_mask) ! Note: helio_pl does not have its own append method, so we skip back to the base class
      class default
         write(*,*) "Invalid object passed to the append method. Source must be of class symba_pl or its descendents!"
         call base_util_exit(FAILURE)
      end select

      return
   end subroutine symba_util_append_pl


   module subroutine symba_util_append_tp(self, source, lsource_mask)
      !! author: David A. Minton
      !!
      !! Append components from test particle object to another. 
      !! This method will automatically resize the destination body if it is too small
      implicit none
      !! Arguments
      class(symba_tp),                 intent(inout) :: self         !! SyMBA test particle object
      class(swiftest_body),            intent(in)    :: source       !! Source object to append
      logical, dimension(:),           intent(in)    :: lsource_mask !! Logical mask indicating which elements to append to

      select type(source)
      class is (symba_tp)
         call util_append(self%levelg, source%levelg, lsource_mask=lsource_mask)
         call util_append(self%levelm, source%levelm, lsource_mask=lsource_mask)

         call swiftest_util_append_tp(self, source, lsource_mask) ! Note: helio_tp does not have its own append method, so we skip back to the base class
      class default
         write(*,*) "Invalid object passed to the append method. Source must be of class symba_tp or its descendents!"
         call base_util_exit(FAILURE)
      end select

      return
   end subroutine symba_util_append_tp


   module subroutine symba_util_dealloc_pl(self)
      !! author: David A. Minton
      !!
      !! Deallocates all allocatabale arrays
      implicit none
      ! Arguments
      class(symba_pl),  intent(inout) :: self !! SyMBA massive body object

      if (allocated(self%levelg)) deallocate(self%levelg)
      if (allocated(self%levelm)) deallocate(self%levelm)

      call self%helio_pl%dealloc()

      return
   end subroutine symba_util_dealloc_pl


   module subroutine symba_util_dealloc_system(self)
      !! author: David A. Minton
      !!
      !! Deallocates all allocatables and resets all values to defaults. Acts as a base for a finalizer
      implicit none
      ! Arguments
      class(symba_nbody_system), intent(inout) :: self

      self%irec = -1
      call self%helio_nbody_system%dealloc()

      return
   end subroutine symba_util_dealloc_system


   module subroutine symba_util_dealloc_tp(self)
      !! author: David A. Minton
      !!
      !! Deallocates all allocatabale arrays
      implicit none
      ! Arguments
      class(symba_tp),  intent(inout) :: self !! SyMBA test particle object

      if (allocated(self%levelg)) deallocate(self%levelg)
      if (allocated(self%levelm)) deallocate(self%levelm)

      call self%helio_tp%dealloc() 

      return
   end subroutine symba_util_dealloc_tp


   module subroutine symba_util_fill_pl(self, inserts, lfill_list)
      !! author: David A. Minton
      !!
      !! Insert new SyMBA test particle structure into an old one. 
      !! This is the inverse of a fill operation.
      !! 
      implicit none
      ! Arguments
      class(symba_pl),       intent(inout) :: self       !! SyMBA masive body object
      class(swiftest_body),  intent(in)    :: inserts    !! Inserted object 
      logical, dimension(:), intent(in)    :: lfill_list !! Logical array of bodies to merge into the keeps

      associate(keeps => self)
         select type(inserts)
         class is (symba_pl)
            call util_fill(keeps%levelg, inserts%levelg, lfill_list)
            call util_fill(keeps%levelm, inserts%levelm, lfill_list)

            call swiftest_util_fill_pl(keeps, inserts, lfill_list)  ! Note: helio_pl does not have its own fill method, so we skip back to the base class
         class default
            write(*,*) "Invalid object passed to the fill method. Source must be of class symba_pl or its descendents!"
            call base_util_exit(FAILURE)
         end select
      end associate

      return
   end subroutine symba_util_fill_pl


   module subroutine symba_util_fill_tp(self, inserts, lfill_list)
      !! author: David A. Minton
      !!
      !! Insert new SyMBA test particle structure into an old one. 
      !! This is the inverse of a fill operation.
      !! 
      implicit none
      ! Arguments
      class(symba_tp),       intent(inout) :: self       !! SyMBA test particle object
      class(swiftest_body),  intent(in)    :: inserts    !! Inserted object 
      logical, dimension(:), intent(in)    :: lfill_list !! Logical array of bodies to merge into the keeps

      associate(keeps => self)
         select type(inserts)
         class is (symba_tp)
            call util_fill(keeps%nplenc, inserts%nplenc, lfill_list)
            call util_fill(keeps%levelg, inserts%levelg, lfill_list)
            call util_fill(keeps%levelm, inserts%levelm, lfill_list)
            
            call swiftest_util_fill_tp(keeps, inserts, lfill_list) ! Note: helio_tp does not have its own fill method, so we skip back to the base class
         class default
            write(*,*) "Invalid object passed to the fill method. Source must be of class symba_tp or its descendents!"
            call base_util_exit(FAILURE)
         end select
      end associate

      return
   end subroutine symba_util_fill_tp


   module subroutine symba_util_flatten_eucl_plpl(self, param)
      !! author: Jacob R. Elliott and David A. Minton
      !!
      !! Turns i,j indices into k index for use in the Euclidean distance matrix. This also sets the lmtiny flag and computes the
      !! number of interactions that excludes semi-interacting bodies with each other (Gmass < GMTINY).
      !! This method will also sort the bodies in descending order by Mass
      !!
      !! Reference:
      !!
      !!    Mélodie Angeletti, Jean-Marie Bonny, Jonas Koko. Parallel Euclidean distance matrix computation on big datasets *. 
      !!       2019. hal-0204751
      implicit none
      ! Arguments
      class(symba_pl),            intent(inout) :: self  !! SyMBA massive body object
      class(swiftest_parameters), intent(inout) :: param !! Current run configuration parameters
      ! Internals
      integer(I8B) :: npl, nplm

      associate(pl => self, nplplm => self%nplplm)
         npl = int(self%nbody, kind=I8B)
         if (param%lmtiny_pl) then 
            pl%lmtiny(1:npl) = pl%Gmass(1:npl) < param%GMTINY 
            nplm = count(.not. pl%lmtiny(1:npl))
         else
            nplm = npl
         end if
         pl%nplm = int(nplm, kind=I4B)
         nplplm = nplm * npl - nplm * (nplm + 1_I8B) / 2_I8B ! number of entries in a strict lower triangle, npl x npl, minus first column including only mutually interacting bodies

         call swiftest_util_flatten_eucl_plpl(pl, param)
      end associate

      return
   end subroutine symba_util_flatten_eucl_plpl


   module subroutine symba_util_resize_pl(self, nnew)
      !! author: David A. Minton
      !!
      !! Checks the current size of a SyMBA massive body object against the requested size and resizes it if it is too small.
      implicit none
      ! Arguments
      class(symba_pl), intent(inout) :: self  !! SyMBA massive body object
      integer(I4B),    intent(in)    :: nnew  !! New size neded

      call util_resize(self%levelg, nnew)
      call util_resize(self%levelm, nnew)

      call swiftest_util_resize_pl(self, nnew)

      return
   end subroutine symba_util_resize_pl

   module subroutine symba_util_resize_tp(self, nnew)
      !! author: David A. Minton
      !!
      !! Checks the current size of a test particle object against the requested size and resizes it if it is too small.
      implicit none
      ! Arguments
      class(symba_tp), intent(inout) :: self  !! SyMBA test particle object
      integer(I4B),    intent(in)    :: nnew  !! New size neded

      call util_resize(self%levelg, nnew)
      call util_resize(self%levelm, nnew)

      call swiftest_util_resize_tp(self, nnew)

      return
   end subroutine symba_util_resize_tp

   module subroutine symba_util_set_renc(self, scale)
      !! author: David A. Minton
      !!
      !! Sets the critical radius for encounter given an input recursion depth
      !!
      implicit none
      ! Arguments
      class(symba_pl), intent(inout) :: self !! SyMBA massive body object
      integer(I4B),    intent(in)    :: scale !! Current recursion depth
      ! Internals
      integer(I4B) :: i
      real(DP)     :: rshell_irec

      associate(pl => self, npl => self%nbody)
         rshell_irec = 1._DP
         do i = 1, scale
            rshell_irec = rshell_irec * RSHELL
         end do
         pl%renc(1:npl) = pl%rhill(1:npl) * RHSCALE * rshell_irec
      end associate

      return
   end subroutine symba_util_set_renc


   module subroutine symba_util_setup_initialize_system(self, system_history, param)
      !! author: David A. Minton
      !!
      !! Initialize an SyMBA nbody system from files and sets up the planetocentric structures.
      !! This subroutine will also sort the massive bodies in descending order by mass
      !! 
      implicit none
      ! Arguments
      class(symba_nbody_system),               intent(inout) :: self           !! SyMBA nbody_system object
      class(swiftest_storage),    allocatable, intent(inout) :: system_history !! Stores the system history between output dumps
      class(swiftest_parameters),              intent(inout) :: param          !! Current run configuration parameters 
      ! Internals
      type(encounter_storage)  :: encounter_history
      type(collision_storage)  :: collision_history

      call encounter_history%setup(4096)
      call collision_history%setup(4096)
      ! Call parent method
      associate(nbody_system => self)
         call helio_util_setup_initialize_system(nbody_system, system_history, param)
         call nbody_system%pltp_encounter%setup(0_I8B)
         call nbody_system%plpl_encounter%setup(0_I8B)
         call nbody_system%plpl_collision%setup(0_I8B)

         if (param%lenc_save_trajectory .or. param%lenc_save_closest) then
            allocate(encounter_netcdf_parameters :: encounter_history%nc)
            select type(nc => encounter_history%nc)
            class is (encounter_netcdf_parameters)
               nc%file_name = ENCOUNTER_OUTFILE
               if (.not.param%lrestart) then
                  call nc%initialize(param)
                  call nc%close()
               end if
            end select
            allocate(nbody_system%encounter_history, source=encounter_history)
         end if
        
         allocate(collision_netcdf_parameters :: collision_history%nc)
         select type(nc => collision_history%nc)
         class is (collision_netcdf_parameters)
            nc%file_name = COLLISION_OUTFILE
            if (param%lrestart) then
               call nc%open(param) ! This will find the nc%max_idslot variable
            else
               call nc%initialize(param)
            end if
            call nc%close()
         end select
         allocate(nbody_system%collision_history, source=collision_history)

         select case(param%collision_model)
         case("MERGE")
            allocate(collision_basic :: nbody_system%collider)
         case("BOUNCE")
            allocate(collision_bounce :: nbody_system%collider)
         case("FRAGGLE")
            allocate(collision_fraggle :: nbody_system%collider)
         end select
         call nbody_system%collider%setup(nbody_system)

         nbody_system%collider%max_rot = MAX_ROT_SI * param%TU2S
         select type(nc => collision_history%nc)
         class is (collision_netcdf_parameters)
            nbody_system%collider%maxid_collision = nc%max_idslot
         end select

      end associate

      return
   end subroutine symba_util_setup_initialize_system


   module subroutine symba_util_setup_pl(self, n, param)
      !! author: David A. Minton
      !!
      !! Allocate SyMBA test particle structure
      !!
      !! Equivalent in functionality to David E. Kaufmann's Swifter routine symba_util_setup.f90
      implicit none
      ! Arguments
      class(symba_pl),            intent(inout) :: self  !! SyMBA massive body object
      integer(I4B),               intent(in)    :: n     !! Number of particles to allocate space for
      class(swiftest_parameters), intent(in)    :: param !! Current run configuration parameter

      !> Call allocation method for parent class. 
      call self%helio_pl%setup(n, param) 
      if (n == 0) return

      allocate(self%levelg(n))
      allocate(self%levelm(n))

      self%levelg(:) = -1
      self%levelm(:) = -1
      return
   end subroutine symba_util_setup_pl


   module subroutine symba_util_setup_tp(self, n, param)
      !! author: David A. Minton
      !!
      !! Allocate WHM test particle structure
      !!
      !! Equivalent in functionality to David E. Kaufmann's Swifter routine whm_util_setup.f90
      implicit none
      ! Arguments
      class(symba_tp),            intent(inout) :: self  !! SyMBA test particle object
      integer(I4B),               intent(in)    :: n     !! Number of particles to allocate space for
      class(swiftest_parameters), intent(in)    :: param !! Current run configuration parameter

      !> Call allocation method for parent class. 
      call self%helio_tp%setup(n, param) 
      if (n == 0) return

      allocate(self%levelg(n))
      allocate(self%levelm(n))

      self%levelg(:) = -1
      self%levelm(:) = -1
      
      return
   end subroutine symba_util_setup_tp


   module subroutine symba_util_sort_pl(self, sortby, ascending)
      !! author: David A. Minton
      !!
      !! Sort a SyMBA massive body object in-place. 
      !! sortby is a string indicating which array component to sort.
      implicit none
      ! Arguments
      class(symba_pl), intent(inout) :: self      !! SyMBA massive body object
      character(*),    intent(in)    :: sortby    !! Sorting attribute
      logical,         intent(in)    :: ascending !! Logical flag indicating whether or not the sorting should be in ascending or descending order
      ! Internals
      integer(I4B), dimension(:), allocatable :: ind
      integer(I4B)                        :: direction

      if (self%nbody == 0) return

      if (ascending) then
         direction = 1
      else
         direction = -1
      end if

      associate(pl => self, npl => self%nbody)
         select case(sortby)
         case("levelg")
            call util_sort(direction * pl%levelg(1:npl), ind)
         case("levelm")
            call util_sort(direction * pl%levelm(1:npl), ind)

         case default ! Look for components in the parent class
            call swiftest_util_sort_pl(pl, sortby, ascending)
            return
         end select

         call pl%rearrange(ind)

      end associate
      return
   end subroutine symba_util_sort_pl


   module subroutine symba_util_sort_tp(self, sortby, ascending)
      !! author: David A. Minton
      !!
      !! Sort a SyMBA test particle object in-place. 
      !! sortby is a string indicating which array component to sort.
      implicit none
      ! Arguments
      class(symba_tp), intent(inout) :: self      !! SyMBA test particle object
      character(*),    intent(in)    :: sortby    !! Sorting attribute
      logical,         intent(in)    :: ascending !! Logical flag indicating whether or not the sorting should be in ascending or descending order
      ! Internals
      integer(I4B), dimension(:), allocatable :: ind
      integer(I4B)                        :: direction

      if (self%nbody == 0) return

      if (ascending) then
         direction = 1
      else
         direction = -1
      end if

      associate(tp => self, ntp => self%nbody)
         select case(sortby)
         case("nplenc")
            call util_sort(direction * tp%nplenc(1:ntp), ind)
         case("levelg")
            call util_sort(direction * tp%levelg(1:ntp), ind)
         case("levelm")
            call util_sort(direction * tp%levelm(1:ntp), ind)
         case default ! Look for components in the parent class
            call swiftest_util_sort_tp(tp, sortby, ascending)
            return
         end select

         call tp%rearrange(ind)
      end associate

      return
   end subroutine symba_util_sort_tp


   module subroutine symba_util_sort_rearrange_pl(self, ind)
      !! author: David A. Minton
      !!
      !! Rearrange SyMBA massive body structure in-place from an index list.
      !! This is a helper utility used to make polymorphic sorting work on Swiftest structures.
      implicit none
      ! Arguments
      class(symba_pl),               intent(inout) :: self !! SyMBA massive body object
      integer(I4B),    dimension(:), intent(in)    :: ind  !! Index array used to restructure the body (should contain all 1:n index values in the desired order)

      associate(pl => self, npl => self%nbody)
         call util_sort_rearrange(pl%levelg,     ind, npl)
         call util_sort_rearrange(pl%levelm,     ind, npl)
         call swiftest_util_sort_rearrange_pl(pl,ind)
      end associate

      return
   end subroutine symba_util_sort_rearrange_pl


   module subroutine symba_util_sort_rearrange_tp(self, ind)
      !! author: David A. Minton
      !!
      !! Rearrange SyMBA test particle object in-place from an index list.
      !! This is a helper utility used to make polymorphic sorting work on Swiftest structures.
      implicit none
      ! Arguments
      class(symba_tp),               intent(inout) :: self !! SyMBA test particle object
      integer(I4B),    dimension(:), intent(in)    :: ind  !! Index array used to restructure the body (should contain all 1:n index values in the desired order)

      associate(tp => self, ntp => self%nbody)
         call util_sort_rearrange(tp%nplenc, ind, ntp)
         call util_sort_rearrange(tp%levelg, ind, ntp)
         call util_sort_rearrange(tp%levelm, ind, ntp)

         call swiftest_util_sort_rearrange_tp(tp,ind)
      end associate
      
      return
   end subroutine symba_util_sort_rearrange_tp


   module subroutine symba_util_spill_pl(self, discards, lspill_list, ldestructive)
      !! author: David A. Minton
      !!
      !! Move spilled (discarded) SyMBA massive body particle structure from active list to discard list
      !! Adapted from David E. Kaufmann's Swifter routine whm_discard_spill.f90
      implicit none
      ! Arguments
      class(symba_pl),       intent(inout) :: self        !! SyMBA massive body object
      class(swiftest_body),  intent(inout) :: discards    !! Discarded object 
      logical, dimension(:), intent(in)    :: lspill_list !! Logical array of bodies to spill into the discards
      logical,               intent(in)    :: ldestructive !! Logical flag indicating whether or not this operation should alter body by removing the discard list

      ! For each component, pack the discarded bodies into the discard object and do the inverse with the keeps
      !> Spill all the common components
      associate(keeps => self)
         select type(discards)
         class is (symba_pl)
            call util_spill(keeps%levelg, discards%levelg, lspill_list, ldestructive)
            call util_spill(keeps%levelm, discards%levelm, lspill_list, ldestructive)

            call swiftest_util_spill_pl(keeps, discards, lspill_list, ldestructive)
         class default
            write(*,*) "Invalid object passed to the spill method. Source must be of class symba_pl or its descendents!"
            call base_util_exit(FAILURE)
         end select
      end associate
     
      return
   end subroutine symba_util_spill_pl


   module subroutine symba_util_spill_tp(self, discards, lspill_list, ldestructive)
      !! author: David A. Minton
      !!
      !! Move spilled (discarded) SyMBA test particle structure from active list to discard list
      !! Adapted from David E. Kaufmann's Swifter routine whm_discard_spill.f90
      implicit none
      ! Arguments
      class(symba_tp),       intent(inout) :: self         !! SyMBA test particle object
      class(swiftest_body),  intent(inout) :: discards     !! Discarded object 
      logical, dimension(:), intent(in)    :: lspill_list  !! Logical array of bodies to spill into the discards
      logical,               intent(in)    :: ldestructive !! Logical flag indicating whether or not this operation should alter body by removing the discard list

      ! For each component, pack the discarded bodies into the discard object and do the inverse with the keeps
      !> Spill all the common components
      associate(keeps => self)
         select type(discards)
         class is (symba_tp)
            call util_spill(keeps%nplenc, discards%nplenc, lspill_list, ldestructive)
            call util_spill(keeps%levelg, discards%levelg, lspill_list, ldestructive)
            call util_spill(keeps%levelm, discards%levelm, lspill_list, ldestructive)

            call swiftest_util_spill_tp(keeps, discards, lspill_list, ldestructive)
         class default
            write(*,*) "Invalid object passed to the spill method. Source must be of class symba_tp or its descendents!"
            call base_util_exit(FAILURE)
         end select
      end associate
     
      return
   end subroutine symba_util_spill_tp



end submodule s_symba_util
