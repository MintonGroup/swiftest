!! Copyright 2022 - David Minton, Carlisle Wishard, Jennifer Pouplin, Jake Elliott, & Dana Singh
!! This file is part of Swiftest.
!! Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
!! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!! Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
!! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!! You should have received a copy of the GNU General Public License along with Swiftest. 
!! If not, see: https://www.gnu.org/licenses. 

submodule (netcdf_io) s_netcdf_io_implementations
   use netcdf
contains

   module subroutine netcdf_io_check(status, call_identifier)
      !! author: Carlisle A. Wishard, Dana Singh, and David A. Minton
      !!
      !! Checks the status of all NetCDF operations to catch errors
      implicit none
      ! Arguments
      integer, intent (in) :: status !! The status code returned by a NetCDF function
      character(len=*), intent(in), optional :: call_identifier !! String that indicates which calling function caused the error for diagnostic purposes

      if(status /= NF90_NOERR) then
         if (present(call_identifier)) write(*,*) "NetCDF error in ",trim(call_identifier)
         write(*,*) trim(nf90_strerror(status))
         call base_util_exit(FAILURE)
      end if

      return
   end subroutine netcdf_io_check


   module subroutine netcdf_io_close(self)
      !! author: Carlisle A. Wishard, Dana Singh, and David A. Minton
      !!
      !! Closes a NetCDF file
      implicit none
      ! Arguments
      class(netcdf_parameters),   intent(inout) :: self   !! Parameters used to identify a particular NetCDF dataset
      character(namelen) :: message

      if (self%lfile_is_open) then
#ifdef COARRAY
         write(message,*) this_image()
         message = "netcdf_io_close on image " // trim(adjustl(message))
#else
         message = "netcdf_io_close"
#endif
         call netcdf_io_check( nf90_close(self%id), message)
         self%lfile_is_open = .false.
      end if

      return
   end subroutine netcdf_io_close


   module subroutine netcdf_io_find_tslot(self, t, tslot)
      !! author: David A. Minton
      !! 
      !! Given an open NetCDF file and a value of time t, finds the index of the time value (aka the time slot) to place a new set of data.
      !! The returned value of tslot will correspond to the first index value where the value of t is greater than or equal to the saved time value.
      use, intrinsic :: ieee_exceptions
      use, intrinsic :: ieee_arithmetic
      implicit none
      ! Arguments
      class(netcdf_parameters), intent(inout) :: self  !! Parameters used to identify a particular NetCDF dataset
      real(DP),                 intent(in)    :: t     !! The value of time to search for
      integer(I4B),             intent(out)   :: tslot !! The index of the time slot where this data belongs
      ! Internals
      real(DP), dimension(:), allocatable :: tvals
      integer(I4B) :: i
      logical, dimension(size(IEEE_ALL))      :: fpe_halting_modes

      call ieee_get_halting_mode(IEEE_ALL,fpe_halting_modes)  ! Save the current halting modes so we can turn them off temporarily
      call ieee_set_halting_mode(IEEE_ALL,.false.)

      if (.not.self%lfile_is_open) return
      tslot = 0

      call netcdf_io_check( nf90_inquire_dimension(self%id, self%time_dimid, self%time_dimname, len=self%max_tslot), "netcdf_io_find_tslot nf90_inquire_dimension max_tslot"  )
      if (self%max_tslot > 0) then
         allocate(tvals(self%max_tslot))
         call netcdf_io_check( nf90_get_var(self%id, self%time_varid, tvals(:), start=[1]), "netcdf_io_find_tslot get_var"  )
         where(.not.ieee_is_normal(tvals(:))) tvals(:) = huge(1.0_DP)
      else
         allocate(tvals(1))
         tvals(1) = huge(1.0_DP)
         self%max_tslot = 1
      end if

      tslot = 1
      do i = 1, self%max_tslot
         if (t <= tvals(tslot)) exit
         tslot = tslot + 1
      end do
      self%max_tslot = max(self%max_tslot, tslot)
      self%tslot = tslot

      call ieee_set_halting_mode(IEEE_ALL,fpe_halting_modes)

      return
   end subroutine netcdf_io_find_tslot


   module subroutine netcdf_io_find_idslot(self, id, idslot)
      !! author: David A. Minton
      !! 
      !! Given an open NetCDF file and a value of id, finds the index of the id value (aka the id slot) to place a new set of data.
      !! The returned value of idslot will correspond to the first index value where the value of id is greater than or equal to the saved id value.
      implicit none
      ! Arguments
      class(netcdf_parameters), intent(inout) :: self   !! Parameters used to identify a particular NetCDF dataset
      integer(I4B),             intent(in)    :: id     !! The value of id to search for
      integer(I4B),             intent(out)   :: idslot !! The index of the id slot where this data belongs
      ! Internals
      integer(I4B), dimension(:), allocatable :: idvals

      if (.not.self%lfile_is_open) return

      if (.not.allocated(self%idvals)) call self%get_idvals()
      self%max_idslot = size(self%idvals)
      idslot = id + 1 
      if (idslot > self%max_idslot) then

         ! Update the idvals array
         allocate(idvals(idslot))
         idvals(:) = NF90_FILL_INT 
         idvals(1:self%max_idslot) = self%idvals(1:self%max_idslot)
         idvals(idslot) = id
         call move_alloc(idvals, self%idvals) 
         self%max_idslot = idslot
      end if

      self%idslot = idslot

      return
   end subroutine netcdf_io_find_idslot


   module subroutine netcdf_io_get_idvals(self)
      !! author: David A. Minton
      !! 
      !! Gets a full list of id values 
      implicit none
      ! Arguments
      class(netcdf_parameters),                            intent(inout) :: self   !! Parameters used to identify a particular NetCDF dataset
      ! Internals

      if (.not.self%lfile_is_open) return

      if (allocated(self%idvals)) deallocate(self%idvals)
      call netcdf_io_check( nf90_inquire_dimension(self%id, self%name_dimid, self%name_dimname, len=self%max_idslot), "netcdf_io_find_tslot nf90_inquire_dimension max_tslot"  )
      if (self%max_idslot > 0) then
         allocate(self%idvals(self%max_idslot))
         call netcdf_io_check( nf90_get_var(self%id, self%id_varid, self%idvals(:), start=[1]), "netcdf_io_find_idslot get_var"  )
      else
         allocate(self%idvals(1))
         self%idvals(1) = 0
      end if

      return
   end subroutine netcdf_io_get_idvals


   module subroutine netcdf_io_sync(self)
      !! author: David A. Minton
      !!
      !! Syncrhonize the disk and memory buffer of the NetCDF file (e.g. commit the frame files stored in memory to disk) 
      !!    
      implicit none
      ! Arguments
      class(netcdf_parameters), intent(inout) :: self !! Parameters used to identify a particular NetCDF dataset

      call netcdf_io_check( nf90_sync(self%id), "netcdf_io_sync nf90_sync"  )

      return
   end subroutine netcdf_io_sync



end submodule s_netcdf_io_implementations
