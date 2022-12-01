!! Copyright 2022 - David Minton, Carlisle Wishard, Jennifer Pouplin, Jake Elliott, & Dana Singh
!! This file is part of Swiftest.
!! Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
!! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!! Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
!! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!! You should have received a copy of the GNU General Public License along with Swiftest. 
!! If not, see: https://www.gnu.org/licenses. 

submodule(swiftest_classes) s_util_copy
   use swiftest
contains

   module subroutine util_copy_particle_info(self, source)
      !! author: David A. Minton
      !!
      !! Copies one set of information object components into another, component-by-component
      implicit none
      class(swiftest_particle_info),  intent(inout) :: self
      class(swiftest_particle_info),  intent(in)    :: source

      call self%set_value(&
         name = source%name, &
         particle_type = source%particle_type, &
         status = source%status, & 
         origin_type = source%origin_type, &
         origin_time = source%origin_time, & 
         collision_id = source%collision_id, &
         origin_xh = source%origin_xh(:), &
         origin_vh = source%origin_vh(:), &
         discard_time = source%discard_time, & 
         discard_xh = source%discard_xh(:), &
         discard_vh = source%discard_vh(:), &
         discard_body_id = source%discard_body_id &
      )

      return
   end subroutine util_copy_particle_info


   module subroutine util_copy_particle_info_arr(source, dest, idx)
      !! author: David A. Minton
      !!
      !! Copies contents from an array of one particle information objects to another.
      implicit none
      class(swiftest_particle_info), dimension(:), intent(in)             :: source !! Source object to copy into
      class(swiftest_particle_info), dimension(:), intent(inout)          :: dest   !! Swiftest body object with particle metadata information object
      integer(I4B),                  dimension(:), intent(in),   optional :: idx    !! Optional array of indices to draw the source object
      ! Internals
      integer(I4B) :: i, j, n, nsource, ndest

      if (size(source) == 0) return

      if (present(idx)) then
         n = size(idx)
      else
         n = size(source)
      end if

      nsource = size(source)
      ndest = size(dest)

      if ((n == 0) .or. (n > ndest) .or. (n > nsource)) then
         write(*,*) 'Particle info copy operation failed. n, nsource, ndest: ',n, nsource, ndest
         return
      end if

      do i = 1, n
         if (present(idx)) then
            j = idx(i)
         else
            j = i
         end if 
         call dest(i)%copy(source(j))
      end do
      
      return
   end subroutine util_copy_particle_info_arr


   module subroutine util_copy_store_system(self, system)
      !! author: David A. Minton
      !!
      !! Stores a snapshot of the nbody system so that later it can be retrieved for saving to file.
      implicit none
      class(storage_frame),         intent(inout) :: self   !! Swiftest storage frame object
      class(swiftest_nbody_system), intent(in)    :: system !! Swiftest n-body system object

      if (allocated(self%system)) deallocate(self%system)
      allocate(self%system, source=system)
      return

   end subroutine util_copy_store_system

end submodule s_util_copy
