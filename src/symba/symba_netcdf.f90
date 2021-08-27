submodule (symba_classes) s_symba_netcdf
   use swiftest
   use netcdf
contains

   subroutine check(status)
      !! author: Carlisle A. Wishard, Dana Singh, and David A. Minton
      !!
      !! Checks the status of all NetCDF operations to catch errors
      implicit none
      ! Arguments
      integer, intent (in) :: status

      if(status /= nf90_noerr) then
         write(*,*) trim(nf90_strerror(status))
         call util_exit(FAILURE)
      end if

      return
   end subroutine check

   module subroutine symba_netcdf_initialize_output(self, param)
      !! author: Carlisle A. Wishard, Dana Singh, and David A. Minton
      !!
      !! Initialize a NetCDF file system and defines all variables.
      implicit none
      ! Arguments
      class(symba_netcdf_parameters), intent(inout) :: self  !! Parameters used to identify a particular NetCDF dataset
      class(swiftest_parameters),     intent(in)    :: param !! Current run configuration parameters 


      call netcdf_initialize_output(self, param)

      ! Define the variables
      call check( nf90_def_var(self%ncid, ORIGIN_TYPE_VARNAME, NF90_CHAR, [self%str_dimid, self%id_dimid], self%origin_type_varid) )
      call check( nf90_def_var(self%ncid, ORIGIN_TIME_VARNAME, self%out_type, self%id_dimid, self%origin_time_varid) )
      call check( nf90_def_var(self%ncid, ORIGIN_XHX_VARNAME, self%out_type, self%id_dimid, self%origin_xhx_varid) )
      call check( nf90_def_var(self%ncid, ORIGIN_XHY_VARNAME, self%out_type, self%id_dimid, self%origin_xhy_varid) )
      call check( nf90_def_var(self%ncid, ORIGIN_XHZ_VARNAME, self%out_type, self%id_dimid, self%origin_xhz_varid) )
      call check( nf90_def_var(self%ncid, ORIGIN_VHX_VARNAME, self%out_type, self%id_dimid, self%origin_vhx_varid) )
      call check( nf90_def_var(self%ncid, ORIGIN_VHY_VARNAME, self%out_type, self%id_dimid, self%origin_vhy_varid) )
      call check( nf90_def_var(self%ncid, ORIGIN_VHZ_VARNAME, self%out_type, self%id_dimid,  self%origin_vhz_varid) )

      return

   end subroutine symba_netcdf_initialize_output


   module subroutine symba_netcdf_open(self, param)
      !! author: Carlisle A. Wishard, Dana Singh, and David A. Minton
      !!
      !! Opens a NetCDF file and does the variable inquiries to activate variable ids
      implicit none
      ! Arguments
      class(symba_netcdf_parameters),   intent(inout) :: self   !! Parameters used to identify a particular NetCDF dataset
      class(swiftest_parameters), intent(in)    :: param  !! Current run configuration parameters
      ! Internals
      integer(I4B) :: old_mode

      call netcdf_open(self, param)

      call check( nf90_inq_varid(self%ncid, ORIGIN_TYPE_VARNAME, self%origin_type_varid))
      call check( nf90_inq_varid(self%ncid, ORIGIN_TIME_VARNAME, self%origin_type_varid))
      call check( nf90_inq_varid(self%ncid, ORIGIN_XHX_VARNAME, self%origin_xhx_varid))
      call check( nf90_inq_varid(self%ncid, ORIGIN_XHY_VARNAME, self%origin_xhy_varid))
      call check( nf90_inq_varid(self%ncid, ORIGIN_XHZ_VARNAME, self%origin_xhz_varid))
      call check( nf90_inq_varid(self%ncid, ORIGIN_VHX_VARNAME, self%origin_vhx_varid))
      call check( nf90_inq_varid(self%ncid, ORIGIN_VHY_VARNAME, self%origin_vhy_varid))
      call check( nf90_inq_varid(self%ncid, ORIGIN_VHZ_VARNAME, self%origin_vhz_varid))

      return
   end subroutine symba_netcdf_open

   
   module subroutine symba_netcdf_write_frame_pl(self, iu, param)
      !! author: Carlisle A. Wishard, Dana Singh, and David A. Minton
      !!
      !! Write a frame of output of a SyMBA massive body data to the binary output file
      !!    Note: If outputting to orbital elements, but sure that the conversion is done prior to calling this method
      implicit none
      ! Arguments
      class(symba_pl),            intent(in)    :: self   !! SyMBA massive body object
      class(netcdf_parameters),   intent(inout) :: iu     !! Parameters used to identify a particular NetCDF dataset
      class(swiftest_parameters), intent(in)    :: param  !! Current run configuration parameters
      ! Internals
      integer(I4B)                              :: i, j, tslot, strlen, idslot
      integer(I4B), dimension(:), allocatable   :: ind
      character(len=:), allocatable             :: charstring

      call netcdf_write_frame_base(self, iu, param)
      tslot = int(param%ioutput, kind=I4B) + 1
      select type(iu)
      class is (symba_netcdf_parameters)
         associate(npl => self%nbody)
            if (npl == 0) return
            allocate(ind(npl))
            call util_sort(self%id(1:npl), ind(1:npl))
            select type(info => self%info)
            class is (symba_particle_info)
               do i = 1, npl
                  j = ind(i)
                  idslot = self%id(j) + 1
      
                  charstring = trim(adjustl(info(j)%origin_type))
                  strlen = len(charstring)
                  call check( nf90_put_var(iu%ncid, iu%origin_type_varid, charstring, start=[1, idslot], count=[strlen, 1]) )
               end do
            end select
         end associate

      end select

      return
   end subroutine symba_netcdf_write_frame_pl
end submodule s_symba_netcdf