submodule (swiftest_classes) s_netcdf
   use swiftest
   use netcdf
contains

   module subroutine netcdf_write_frame_base(self, iu, param)
      !! author: Carlisle A. Wishard, Dana Singh, and David A. Minton
      !!
      !! Write a frame of output of either test particle or massive body data to the binary output file
      !!    Note: If outputting to orbital elements, but sure that the conversion is done prior to calling this method
      !!
      !! Adapted from David E. Kaufmann's Swifter routine  io_write_frame.f90
      !! Adapted from Hal Levison's Swift routine io_write_frame.F
      implicit none
      ! Arguments
      class(swiftest_base),       intent(in)    :: self   !! Swiftest particle object
      class(netcdf_parameters),   intent(inout) :: iu     !! Parameters used to identify a particular NetCDF dataset
      class(swiftest_parameters), intent(in)    :: param  !! Current run configuration parameters
      ! Internals
      integer(I4B)                              :: i, j, id, ioutput
      integer(I4B), dimension(:), allocatable   :: ind

      !! Open the netCDF file
      call check( nf90_open(param%outfile, nf90_write, iu%ncid) )
      !call check( nf90_set_fill(iu%ncid, nf90_nofill, oldMode) )


      ! Calculate the output number that we are currently on
      ioutput = (param%t / param%dt) / param%istep_out

      select type(self)
      class is (swiftest_body)
         associate(n => self%nbody)
            if (n == 0) return
            allocate(ind(n))
            call util_sort(self%id(1:n), ind(1:n))

            do i = 1, n
               j = ind(i)
               id = self%id(j)
               if ((param%out_form == XV) .or. (param%out_form == XVEL)) then
                  call check( nf90_inq_varid(iu%ncid, XHX_VARNAME, iu%xhx_varid))
                  call check( nf90_inq_varid(iu%ncid, XHY_VARNAME, iu%xhy_varid))
                  call check( nf90_inq_varid(iu%ncid, XHZ_VARNAME, iu%xhz_varid))
                  call check( nf90_inq_varid(iu%ncid, VHX_VARNAME, iu%vhx_varid))
                  call check( nf90_inq_varid(iu%ncid, VHY_VARNAME, iu%vhy_varid))
                  call check( nf90_inq_varid(iu%ncid, VHZ_VARNAME, iu%vhz_varid))
      
                  call check( nf90_put_var(iu%ncid, iu%xhx_varid, self%xh(1, j), start=[ioutput + 1, id]) )
                  call check( nf90_put_var(iu%ncid, iu%xhy_varid, self%xh(2, j), start=[ioutput + 1, id]) )
                  call check( nf90_put_var(iu%ncid, iu%xhz_varid, self%xh(3, j), start=[ioutput + 1, id]) )
                  call check( nf90_put_var(iu%ncid, iu%vhx_varid, self%vh(1, j), start=[ioutput + 1, id]) )
                  call check( nf90_put_var(iu%ncid, iu%vhy_varid, self%vh(2, j), start=[ioutput + 1, id]) )
                  call check( nf90_put_var(iu%ncid, iu%vhz_varid, self%vh(3, j), start=[ioutput + 1, id]) )
               end if
               if ((param%out_form == EL) .or. (param%out_form == XVEL)) then
                  call check( nf90_inq_varid(iu%ncid, A_VARNAME, iu%a_varid))
                  call check( nf90_inq_varid(iu%ncid, E_VARNAME, iu%e_varid))
                  call check( nf90_inq_varid(iu%ncid, INC_VARNAME, iu%inc_varid))
                  call check( nf90_inq_varid(iu%ncid, CAPOM_VARNAME, iu%capom_varid))
                  call check( nf90_inq_varid(iu%ncid, OMEGA_VARNAME, iu%omega_varid))
                  call check( nf90_inq_varid(iu%ncid, CAPM_VARNAME, iu%capm_varid))

                  call check( nf90_put_var(iu%ncid, iu%a_varid, self%a(j), start=[ioutput + 1, id]) )
                  call check( nf90_put_var(iu%ncid, iu%e_varid, self%e(j), start=[ioutput + 1, id]) )
                  call check( nf90_put_var(iu%ncid, iu%inc_varid, self%inc(j), start=[ioutput + 1, id]) )
                  call check( nf90_put_var(iu%ncid, iu%capom_varid, self%capom(j), start=[ioutput + 1, id]) )
                  call check( nf90_put_var(iu%ncid, iu%omega_varid, self%omega(j), start=[ioutput + 1, id]) )
                  call check( nf90_put_var(iu%ncid, iu%capm_varid, self%capm(j), start=[ioutput + 1, id]) ) 
               end if
               select type(pl => self)  
               class is (swiftest_pl)  ! Additional output if the passed polymorphic object is a massive body
                  call check( nf90_inq_varid(iu%ncid, GMASS_VARNAME, iu%Gmass_varid))
                  call check( nf90_put_var(iu%ncid, iu%Gmass_varid, pl%Gmass(j), start=[ioutput + 1, id]) )
                  if (param%lrhill_present) then 
                     call check( nf90_inq_varid(iu%ncid, RHILL_VARNAME, iu%rhill_varid))
                     call check( nf90_put_var(iu%ncid, iu%rhill_varid, pl%rhill(j), start=[ioutput + 1, id]) )
                  end if
                  if (param%lclose) then
                     call check( nf90_inq_varid(iu%ncid, RADIUS_VARNAME, iu%radius_varid))
                     call check( nf90_put_var(iu%ncid, iu%radius_varid, pl%radius(j), start=[ioutput + 1, id]) )
                  end if
                  if (param%lrotation) then
                     call check( nf90_inq_varid(iu%ncid, IP1_VARNAME, iu%Ip1_varid))
                     call check( nf90_inq_varid(iu%ncid, IP2_VARNAME, iu%Ip2_varid))
                     call check( nf90_inq_varid(iu%ncid, IP3_VARNAME, iu%Ip3_varid))
                     call check( nf90_inq_varid(iu%ncid, ROTX_VARNAME, iu%rotx_varid))
                     call check( nf90_inq_varid(iu%ncid, ROTY_VARNAME, iu%roty_varid))
                     call check( nf90_inq_varid(iu%ncid, ROTZ_VARNAME, iu%rotz_varid))

                     call check( nf90_put_var(iu%ncid, iu%Ip1_varid, pl%Ip(1, j), start=[ioutput + 1, id]) )
                     call check( nf90_put_var(iu%ncid, iu%Ip2_varid, pl%Ip(2, j), start=[ioutput + 1, id]) )
                     call check( nf90_put_var(iu%ncid, iu%Ip3_varid, pl%Ip(3, j), start=[ioutput + 1, id]) )
                     call check( nf90_put_var(iu%ncid, iu%rotx_varid, pl%rot(1, j), start=[ioutput + 1, id]) )
                     call check( nf90_put_var(iu%ncid, iu%roty_varid, pl%rot(2, j), start=[ioutput + 1, id]) )
                     call check( nf90_put_var(iu%ncid, iu%rotz_varid, pl%rot(3, j), start=[ioutput + 1, id]) )
                  end if
                  if (param%ltides) then
                     call check( nf90_inq_varid(iu%ncid, K2_VARNAME, iu%k2_varid))
                     call check( nf90_inq_varid(iu%ncid, Q_VARNAME, iu%Q_varid))

                     call check( nf90_put_var(iu%ncid, iu%k2_varid, pl%k2(j), start=[ioutput + 1, id]) )
                     call check( nf90_put_var(iu%ncid, iu%Q_varid, pl%Q(j), start=[ioutput + 1, id]) )
                  end if
               end select
            end do
         end associate
      class is (swiftest_cb)
         id = self%id
         call check( nf90_inq_varid(iu%ncid, GMASS_VARNAME, iu%Gmass_varid))
         call check( nf90_put_var(iu%ncid, iu%Gmass_varid, self%Gmass, start=[ioutput + 1, id]) )
         call check( nf90_inq_varid(iu%ncid, RADIUS_VARNAME, iu%radius_varid))
         call check( nf90_put_var(iu%ncid, iu%radius_varid, self%radius, start=[ioutput + 1, id]) )
         if (param%lrotation) then
            call check( nf90_inq_varid(iu%ncid, IP1_VARNAME, iu%Ip1_varid))
            call check( nf90_inq_varid(iu%ncid, IP2_VARNAME, iu%Ip2_varid))
            call check( nf90_inq_varid(iu%ncid, IP3_VARNAME, iu%Ip3_varid))
            call check( nf90_inq_varid(iu%ncid, ROTX_VARNAME, iu%rotx_varid))
            call check( nf90_inq_varid(iu%ncid, ROTY_VARNAME, iu%roty_varid))
            call check( nf90_inq_varid(iu%ncid, ROTZ_VARNAME, iu%rotz_varid))
   
            call check( nf90_put_var(iu%ncid, iu%Ip1_varid, self%Ip(1), start=[ioutput + 1, id]) )
            call check( nf90_put_var(iu%ncid, iu%Ip2_varid, self%Ip(2), start=[ioutput + 1, id]) )
            call check( nf90_put_var(iu%ncid, iu%Ip3_varid, self%Ip(3), start=[ioutput + 1, id]) )
            call check( nf90_put_var(iu%ncid, iu%rotx_varid, self%rot(1), start=[ioutput + 1, id]) )
            call check( nf90_put_var(iu%ncid, iu%roty_varid, self%rot(2), start=[ioutput + 1, id]) )
            call check( nf90_put_var(iu%ncid, iu%rotz_varid, self%rot(3), start=[ioutput + 1, id]) )
         end if
         if (param%ltides) then
            call check( nf90_inq_varid(iu%ncid, K2_VARNAME, iu%k2_varid))
            call check( nf90_inq_varid(iu%ncid, Q_VARNAME, iu%Q_varid))

            call check( nf90_put_var(iu%ncid, iu%k2_varid, self%k2, start=[ioutput + 1, id]) )
            call check( nf90_put_var(iu%ncid, iu%Q_varid, self%Q, start=[ioutput + 1, id]) )
         end if
      end select

      ! Close the netCDF file
      call check( nf90_close(iu%ncid) )

      return
   end subroutine netcdf_write_frame_base


   module subroutine netcdf_initialize_output(self, param)
      !! author: Carlisle A. Wishard, Dana Singh, and David A. Mintont
      !!
      !! Initialize a NetCDF file system
      !! There is no direct file output from this subroutine
      !!
      !! Adapted from David E. Kaufmann's Swifter routine  io_write_frame.f90
      !! Adapted from Hal Levison's Swift routine io_write_frame.F
      implicit none
      ! Arguments
      class(netcdf_parameters),   intent(inout) :: self    !! Parameters used to identify a particular NetCDF dataset
      class(swiftest_parameters), intent(in)    :: param           !! Current run configuration parameters 
      ! Internals
      logical :: fileExists

      !! Create the new output file, deleting any previously existing output file of the same name
      call check( nf90_create(param%outfile, NF90_HDF5, self%ncid) )
      !call check( nf90_set_fill(self%ncid, nf90_nofill, oldMode) )

      ! Define the NetCDF dimensions with particle name as the record dimension
      call check( nf90_def_dim(self%ncid, ID_DIMNAME, NF90_UNLIMITED, self%id_dimid) )     ! 'x' dimension
      call check( nf90_def_dim(self%ncid, TIME_DIMNAME, NF90_UNLIMITED, self%time_dimid) ) ! 'y' dimension
      self%dimids = [self%time_dimid, self%id_dimid ]

      select case (param%out_type)
      case(NETCDF_FLOAT_TYPE)
         self%out_type = NF90_FLOAT
      case(NETCDF_DOUBLE_TYPE)
         self%out_type = NF90_DOUBLE
      end select

      !! Define the variables
      if ((param%out_form == XV) .or. (param%out_form == XVEL)) then
         call check( nf90_def_var(self%ncid, A_VARNAME, self%out_type, self%dimids, self%a_varid) )
         call check( nf90_def_var(self%ncid, E_VARNAME, self%out_type, self%dimids, self%e_varid) )
         call check( nf90_def_var(self%ncid, INC_VARNAME, self%out_type, self%dimids, self%inc_varid) )
         call check( nf90_def_var(self%ncid, CAPOM_VARNAME, self%out_type, self%dimids, self%capom_varid) )
         call check( nf90_def_var(self%ncid, OMEGA_VARNAME, self%out_type, self%dimids, self%omega_varid) )
         call check( nf90_def_var(self%ncid, CAPM_VARNAME, self%out_type, self%dimids, self%capm_varid) )
      end if
      if ((param%out_form == EL) .or. (param%out_form == XVEL)) then
         call check( nf90_def_var(self%ncid, XHX_VARNAME, self%out_type, self%dimids, self%xhx_varid) )
         call check( nf90_def_var(self%ncid, XHY_VARNAME, self%out_type, self%dimids, self%xhy_varid) )
         call check( nf90_def_var(self%ncid, XHZ_VARNAME, self%out_type, self%dimids, self%xhz_varid) )
         call check( nf90_def_var(self%ncid, VHX_VARNAME, self%out_type, self%dimids, self%vhx_varid) )
         call check( nf90_def_var(self%ncid, VHY_VARNAME, self%out_type, self%dimids, self%vhy_varid) )
         call check( nf90_def_var(self%ncid, VHZ_VARNAME, self%out_type, self%dimids, self%vhz_varid) )
      end if

      call check( nf90_def_var(self%ncid, GMASS_VARNAME, self%out_type, self%dimids, self%Gmass_varid) )
      if (param%lrhill_present) call check( nf90_def_var(self%ncid, RHILL_VARNAME, self%out_type, self%dimids, self%rhill_varid) )
      if (param%lclose) call check( nf90_def_var(self%ncid, RADIUS_VARNAME, self%out_type, self%dimids, self%radius_varid) )
      if (param%lrotation) then
         call check( nf90_def_var(self%ncid, IP1_VARNAME, self%out_type, self%dimids, self%Ip1_varid) )
         call check( nf90_def_var(self%ncid, IP2_VARNAME, self%out_type, self%dimids, self%Ip2_varid) )
         call check( nf90_def_var(self%ncid, IP3_VARNAME, self%out_type, self%dimids, self%Ip3_varid) )
         call check( nf90_def_var(self%ncid, ROTX_VARNAME, self%out_type, self%dimids, self%rotx_varid) )
         call check( nf90_def_var(self%ncid, ROTY_VARNAME, self%out_type, self%dimids, self%roty_varid) )
         call check( nf90_def_var(self%ncid, ROTZ_VARNAME, self%out_type, self%dimids, self%rotz_varid) )
      end if
      if (param%ltides) then
         call check( nf90_def_var(self%ncid, K2_VARNAME, self%out_type, self%dimids, self%k2_varid) )
         call check( nf90_def_var(self%ncid, Q_VARNAME, self%out_type, self%dimids, self%Q_varid) )
      end if

      call check( nf90_close(self%ncid) )

      return
   end subroutine netcdf_initialize_output


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



end submodule s_netcdf