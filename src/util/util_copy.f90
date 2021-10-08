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



end submodule s_util_copy
