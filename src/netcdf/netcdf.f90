submodule (swiftest_classes) s_netcdf
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

   module subroutine netcdf_close(self, param)
      !! author: Carlisle A. Wishard, Dana Singh, and David A. Minton
      !!
      !! Closes a NetCDF file
      implicit none
      ! Arguments
      class(netcdf_parameters),   intent(inout) :: self   !! Parameters used to identify a particular NetCDF dataset
      class(swiftest_parameters), intent(in)    :: param  !! Current run configuration parameters

      call check( nf90_close(self%ncid) )

      return
   end subroutine netcdf_close


   module subroutine netcdf_initialize_output(self, param)
      !! author: Carlisle A. Wishard, Dana Singh, and David A. Minton
      !!
      !! Initialize a NetCDF file system and defines all variables.
      implicit none
      ! Arguments
      class(netcdf_parameters),   intent(inout) :: self    !! Parameters used to identify a particular NetCDF dataset
      class(swiftest_parameters), intent(in)    :: param           !! Current run configuration parameters 
      ! Internals
      logical :: fileExists
      integer(I4B) :: old_mode

      !! Create the new output file, deleting any previously existing output file of the same name
      call check( nf90_create(param%outfile, NF90_NETCDF4, self%ncid) )
      call check( nf90_set_fill(self%ncid, nf90_nofill, old_mode) )

      ! Define the NetCDF dimensions with particle name as the record dimension
      call check( nf90_def_dim(self%ncid, ID_DIMNAME, NF90_UNLIMITED, self%id_dimid) )     ! 'x' dimension
      call check( nf90_def_dim(self%ncid, TIME_DIMNAME, NF90_UNLIMITED, self%time_dimid) ) ! 'y' dimension
      call check( nf90_def_dim(self%ncid, STR_DIMNAME, NAMELEN, self%str_dimid) ) ! Dimension for string variables (aka character arrays)

      select case (param%out_type)
      case(NETCDF_FLOAT_TYPE)
         self%out_type = NF90_FLOAT
      case(NETCDF_DOUBLE_TYPE)
         self%out_type = NF90_DOUBLE
      end select

      !! Define the variables
      call check( nf90_def_var(self%ncid, TIME_DIMNAME, self%out_type, self%time_dimid, self%time_varid) )
      call check( nf90_def_var(self%ncid, ID_DIMNAME, NF90_INT, self%id_dimid, self%id_varid) )
      call check( nf90_def_var(self%ncid, NPL_VARNAME, NF90_INT, self%time_dimid, self%npl_varid) )
      call check( nf90_def_var(self%ncid, NTP_VARNAME, NF90_INT, self%time_dimid, self%ntp_varid) )
      call check( nf90_def_var(self%ncid, NAME_VARNAME, NF90_CHAR, [self%str_dimid, self%id_dimid], self%name_varid) )
      call check( nf90_def_var(self%ncid, PTYPE_VARNAME, NF90_CHAR, [self%str_dimid, self%id_dimid], self%ptype_varid) )
      if ((param%out_form == XV) .or. (param%out_form == XVEL)) then
         call check( nf90_def_var(self%ncid, XHX_VARNAME, self%out_type, [self%id_dimid, self%time_dimid], self%xhx_varid) )
         call check( nf90_def_var(self%ncid, XHY_VARNAME, self%out_type, [self%id_dimid, self%time_dimid], self%xhy_varid) )
         call check( nf90_def_var(self%ncid, XHZ_VARNAME, self%out_type, [self%id_dimid, self%time_dimid], self%xhz_varid) )
         call check( nf90_def_var(self%ncid, VHX_VARNAME, self%out_type, [self%id_dimid, self%time_dimid], self%vhx_varid) )
         call check( nf90_def_var(self%ncid, VHY_VARNAME, self%out_type, [self%id_dimid, self%time_dimid], self%vhy_varid) )
         call check( nf90_def_var(self%ncid, VHZ_VARNAME, self%out_type, [self%id_dimid, self%time_dimid], self%vhz_varid) )
      end if
   
      if ((param%out_form == EL) .or. (param%out_form == XVEL)) then
         call check( nf90_def_var(self%ncid, A_VARNAME, self%out_type, [self%id_dimid, self%time_dimid], self%a_varid) )
         call check( nf90_def_var(self%ncid, E_VARNAME, self%out_type, [self%id_dimid, self%time_dimid], self%e_varid) )
         call check( nf90_def_var(self%ncid, INC_VARNAME, self%out_type, [self%id_dimid, self%time_dimid], self%inc_varid) )
         call check( nf90_def_var(self%ncid, CAPOM_VARNAME, self%out_type, [self%id_dimid, self%time_dimid], self%capom_varid) )
         call check( nf90_def_var(self%ncid, OMEGA_VARNAME, self%out_type, [self%id_dimid, self%time_dimid], self%omega_varid) )
         call check( nf90_def_var(self%ncid, CAPM_VARNAME, self%out_type, [self%id_dimid, self%time_dimid], self%capm_varid) )
      end if

      call check( nf90_def_var(self%ncid, GMASS_VARNAME, self%out_type, [self%id_dimid, self%time_dimid], self%Gmass_varid) )
      if (param%lrhill_present) call check( nf90_def_var(self%ncid, RHILL_VARNAME, self%out_type, [self%id_dimid, self%time_dimid], self%rhill_varid) )
      if (param%lclose) call check( nf90_def_var(self%ncid, RADIUS_VARNAME, self%out_type, [self%id_dimid, self%time_dimid], self%radius_varid) )
      if (param%lrotation) then
         call check( nf90_def_var(self%ncid, IP1_VARNAME, self%out_type, [self%id_dimid, self%time_dimid], self%Ip1_varid) )
         call check( nf90_def_var(self%ncid, IP2_VARNAME, self%out_type, [self%id_dimid, self%time_dimid], self%Ip2_varid) )
         call check( nf90_def_var(self%ncid, IP3_VARNAME, self%out_type, [self%id_dimid, self%time_dimid], self%Ip3_varid) )
         call check( nf90_def_var(self%ncid, ROTX_VARNAME, self%out_type, [self%id_dimid, self%time_dimid], self%rotx_varid) )
         call check( nf90_def_var(self%ncid, ROTY_VARNAME, self%out_type, [self%id_dimid, self%time_dimid], self%roty_varid) )
         call check( nf90_def_var(self%ncid, ROTZ_VARNAME, self%out_type, [self%id_dimid, self%time_dimid], self%rotz_varid) )
      end if
      if (param%ltides) then
         call check( nf90_def_var(self%ncid, K2_VARNAME, self%out_type, [self%id_dimid, self%time_dimid], self%k2_varid) )
         call check( nf90_def_var(self%ncid, Q_VARNAME, self%out_type, [self%id_dimid, self%time_dimid], self%Q_varid) )
      end if
      if (param%lenergy) then
         call check( nf90_def_var(self%ncid, KE_ORB_VARNAME, self%out_type, self%time_dimid, self%KE_orb_varid) )
         call check( nf90_def_var(self%ncid, KE_SPIN_VARNAME, self%out_type, self%time_dimid, self%KE_spin_varid) )
         call check( nf90_def_var(self%ncid, PE_VARNAME, self%out_type, self%time_dimid, self%PE_varid) )
         call check( nf90_def_var(self%ncid, L_ORBX_VARNAME, self%out_type, self%time_dimid, self%L_orbx_varid) )
         call check( nf90_def_var(self%ncid, L_ORBY_VARNAME, self%out_type, self%time_dimid, self%L_orby_varid) )
         call check( nf90_def_var(self%ncid, L_ORBZ_VARNAME, self%out_type, self%time_dimid, self%L_orbz_varid) )
         call check( nf90_def_var(self%ncid, L_SPINX_VARNAME, self%out_type, self%time_dimid, self%L_spinx_varid) )
         call check( nf90_def_var(self%ncid, L_SPINY_VARNAME, self%out_type, self%time_dimid, self%L_spiny_varid) )
         call check( nf90_def_var(self%ncid, L_SPINZ_VARNAME, self%out_type, self%time_dimid, self%L_spinz_varid) )
         call check( nf90_def_var(self%ncid, L_ESCAPEX_VARNAME, self%out_type, self%time_dimid, self%L_escapex_varid) )
         call check( nf90_def_var(self%ncid, L_ESCAPEY_VARNAME, self%out_type, self%time_dimid, self%L_escapey_varid) )
         call check( nf90_def_var(self%ncid, L_ESCAPEZ_VARNAME, self%out_type, self%time_dimid, self%L_escapez_varid) )
         call check( nf90_def_var(self%ncid, ECOLLISIONS_VARNAME, self%out_type, self%time_dimid, self%Ecollisions_varid) )
         call check( nf90_def_var(self%ncid, EUNTRACKED_VARNAME, self%out_type, self%time_dimid, self%Euntracked_varid) )
         call check( nf90_def_var(self%ncid, GMESCAPE_VARNAME, self%out_type, self%time_dimid, self%GMescape_varid) )
      end if

      return
   end subroutine netcdf_initialize_output


   module subroutine netcdf_open(self, param)
      !! author: Carlisle A. Wishard, Dana Singh, and David A. Minton
      !!
      !! Opens a NetCDF file and does the variable inquiries to activate variable ids
      implicit none
      ! Arguments
      class(netcdf_parameters),   intent(inout) :: self   !! Parameters used to identify a particular NetCDF dataset
      class(swiftest_parameters), intent(in)    :: param  !! Current run configuration parameters
      ! Internals
      integer(I4B) :: old_mode

      call check( nf90_open(param%outfile, nf90_write, self%ncid) )
      call check( nf90_set_fill(self%ncid, nf90_nofill, old_mode) )

      call check( nf90_inq_varid(self%ncid, TIME_DIMNAME, self%time_varid))
      call check( nf90_inq_varid(self%ncid, ID_DIMNAME, self%id_varid))
      call check( nf90_inq_varid(self%ncid, NPL_VARNAME, self%npl_varid))
      call check( nf90_inq_varid(self%ncid, NTP_VARNAME, self%ntp_varid))
      call check( nf90_inq_varid(self%ncid, NAME_VARNAME, self%name_varid))
      call check( nf90_inq_varid(self%ncid, PTYPE_VARNAME, self%ptype_varid))

      if ((param%out_form == XV) .or. (param%out_form == XVEL)) then
         call check( nf90_inq_varid(self%ncid, XHX_VARNAME, self%xhx_varid))
         call check( nf90_inq_varid(self%ncid, XHY_VARNAME, self%xhy_varid))
         call check( nf90_inq_varid(self%ncid, XHZ_VARNAME, self%xhz_varid))
         call check( nf90_inq_varid(self%ncid, VHX_VARNAME, self%vhx_varid))
         call check( nf90_inq_varid(self%ncid, VHY_VARNAME, self%vhy_varid))
         call check( nf90_inq_varid(self%ncid, VHZ_VARNAME, self%vhz_varid))
      end if
      if ((param%out_form == EL) .or. (param%out_form == XVEL)) then
         call check( nf90_inq_varid(self%ncid, A_VARNAME, self%a_varid))
         call check( nf90_inq_varid(self%ncid, E_VARNAME, self%e_varid))
         call check( nf90_inq_varid(self%ncid, INC_VARNAME, self%inc_varid))
         call check( nf90_inq_varid(self%ncid, CAPOM_VARNAME, self%capom_varid))
         call check( nf90_inq_varid(self%ncid, OMEGA_VARNAME, self%omega_varid))
         call check( nf90_inq_varid(self%ncid, CAPM_VARNAME, self%capm_varid))
      end if
      call check( nf90_inq_varid(self%ncid, GMASS_VARNAME, self%Gmass_varid))
      if (param%lclose) call check( nf90_inq_varid(self%ncid, RADIUS_VARNAME, self%radius_varid))
      if (param%lrhill_present) call check( nf90_inq_varid(self%ncid, RHILL_VARNAME, self%rhill_varid))

      if (param%lrotation) then
         call check( nf90_inq_varid(self%ncid, IP1_VARNAME, self%Ip1_varid))
         call check( nf90_inq_varid(self%ncid, IP2_VARNAME, self%Ip2_varid))
         call check( nf90_inq_varid(self%ncid, IP3_VARNAME, self%Ip3_varid))
         call check( nf90_inq_varid(self%ncid, ROTX_VARNAME, self%rotx_varid))
         call check( nf90_inq_varid(self%ncid, ROTY_VARNAME, self%roty_varid))
         call check( nf90_inq_varid(self%ncid, ROTZ_VARNAME, self%rotz_varid))
      end if

      if (param%ltides) then
         call check( nf90_inq_varid(self%ncid, K2_VARNAME, self%k2_varid))
         call check( nf90_inq_varid(self%ncid, Q_VARNAME, self%Q_varid))
      end if

      if (param%lenergy) then
         call check( nf90_inq_varid(self%ncid, KE_ORB_VARNAME, self%KE_orb_varid) )
         call check( nf90_inq_varid(self%ncid, KE_SPIN_VARNAME, self%KE_spin_varid) )
         call check( nf90_inq_varid(self%ncid, PE_VARNAME, self%PE_varid) )
         call check( nf90_inq_varid(self%ncid, L_ORBX_VARNAME, self%L_orbx_varid) )
         call check( nf90_inq_varid(self%ncid, L_ORBY_VARNAME, self%L_orby_varid) )
         call check( nf90_inq_varid(self%ncid, L_ORBZ_VARNAME, self%L_orbz_varid) )
         call check( nf90_inq_varid(self%ncid, L_SPINX_VARNAME, self%L_spinx_varid) )
         call check( nf90_inq_varid(self%ncid, L_SPINY_VARNAME, self%L_spiny_varid) )
         call check( nf90_inq_varid(self%ncid, L_SPINZ_VARNAME, self%L_spinz_varid) )
         call check( nf90_inq_varid(self%ncid, L_ESCAPEX_VARNAME, self%L_escapex_varid) )
         call check( nf90_inq_varid(self%ncid, L_ESCAPEY_VARNAME, self%L_escapey_varid) )
         call check( nf90_inq_varid(self%ncid, L_ESCAPEZ_VARNAME, self%L_escapez_varid) )
         call check( nf90_inq_varid(self%ncid, ECOLLISIONS_VARNAME, self%Ecollisions_varid) )
         call check( nf90_inq_varid(self%ncid, EUNTRACKED_VARNAME, self%Euntracked_varid) )
         call check( nf90_inq_varid(self%ncid, GMESCAPE_VARNAME, self%GMescape_varid) )
      end if

      return
   end subroutine netcdf_open


   module subroutine netcdf_write_frame_base(self, iu, param)
      !! author: Carlisle A. Wishard, Dana Singh, and David A. Minton
      !!
      !! Write a frame of output of either test particle or massive body data to the binary output file
      !!    Note: If outputting to orbital elements, but sure that the conversion is done prior to calling this method
      implicit none
      ! Arguments
      class(swiftest_base),       intent(in)    :: self   !! Swiftest particle object
      class(netcdf_parameters),   intent(inout) :: iu     !! Parameters used to identify a particular NetCDF dataset
      class(swiftest_parameters), intent(in)    :: param  !! Current run configuration parameters
      ! Internals
      integer(I4B)                              :: i, j, tslot, strlen, idslot
      integer(I4B), dimension(:), allocatable   :: ind
      character(len=:), allocatable             :: charstring

      tslot = int(param%ioutput, kind=I4B) + 1

      select type(self)
      class is (swiftest_body)
         associate(n => self%nbody)
            if (n == 0) return
            allocate(ind(n))
            call util_sort(self%id(1:n), ind)

            do i = 1, n
               j = ind(i)
               idslot = self%id(j) + 1
               call check( nf90_put_var(iu%ncid, iu%id_varid, self%id(j), start=[idslot]) )
               charstring = trim(adjustl(self%info(j)%name))
               strlen = len(charstring)
               call check( nf90_put_var(iu%ncid, iu%name_varid, charstring, start=[1, idslot], count=[strlen, 1]) )

               charstring = trim(adjustl(self%info(j)%particle_type))
               strlen = len(charstring)
               call check( nf90_put_var(iu%ncid, iu%ptype_varid, charstring, start=[1, idslot], count=[strlen, 1]) )

               if ((param%out_form == XV) .or. (param%out_form == XVEL)) then
                  call check( nf90_put_var(iu%ncid, iu%xhx_varid, self%xh(1, j), start=[idslot, tslot]) )
                  call check( nf90_put_var(iu%ncid, iu%xhy_varid, self%xh(2, j), start=[idslot, tslot]) )
                  call check( nf90_put_var(iu%ncid, iu%xhz_varid, self%xh(3, j), start=[idslot, tslot]) )
                  call check( nf90_put_var(iu%ncid, iu%vhx_varid, self%vh(1, j), start=[idslot, tslot]) )
                  call check( nf90_put_var(iu%ncid, iu%vhy_varid, self%vh(2, j), start=[idslot, tslot]) )
                  call check( nf90_put_var(iu%ncid, iu%vhz_varid, self%vh(3, j), start=[idslot, tslot]) )
               end if
               if ((param%out_form == EL) .or. (param%out_form == XVEL)) then
                  call check( nf90_put_var(iu%ncid, iu%a_varid, self%a(j), start=[idslot, tslot]) )
                  call check( nf90_put_var(iu%ncid, iu%e_varid, self%e(j), start=[idslot, tslot]) )
                  call check( nf90_put_var(iu%ncid, iu%inc_varid, self%inc(j), start=[idslot, tslot]) )
                  call check( nf90_put_var(iu%ncid, iu%capom_varid, self%capom(j), start=[idslot, tslot]) )
                  call check( nf90_put_var(iu%ncid, iu%omega_varid, self%omega(j), start=[idslot, tslot]) )
                  call check( nf90_put_var(iu%ncid, iu%capm_varid, self%capm(j), start=[idslot, tslot]) ) 
               end if
               select type(pl => self)  
               class is (swiftest_pl)  ! Additional output if the passed polymorphic object is a massive body
                  call check( nf90_put_var(iu%ncid, iu%Gmass_varid, pl%Gmass(j), start=[idslot, tslot]) )
                  if (param%lrhill_present) then 
                     call check( nf90_put_var(iu%ncid, iu%rhill_varid, pl%rhill(j), start=[idslot, tslot]) )
                  end if
                  if (param%lclose) then
                     call check( nf90_put_var(iu%ncid, iu%radius_varid, pl%radius(j), start=[idslot, tslot]) )
                  end if
                  if (param%lrotation) then
                     call check( nf90_put_var(iu%ncid, iu%Ip1_varid, pl%Ip(1, j), start=[idslot, tslot]) )
                     call check( nf90_put_var(iu%ncid, iu%Ip2_varid, pl%Ip(2, j), start=[idslot, tslot]) )
                     call check( nf90_put_var(iu%ncid, iu%Ip3_varid, pl%Ip(3, j), start=[idslot, tslot]) )
                     call check( nf90_put_var(iu%ncid, iu%rotx_varid, pl%rot(1, j), start=[idslot, tslot]) )
                     call check( nf90_put_var(iu%ncid, iu%roty_varid, pl%rot(2, j), start=[idslot, tslot]) )
                     call check( nf90_put_var(iu%ncid, iu%rotz_varid, pl%rot(3, j), start=[idslot, tslot]) )
                  end if
                  if (param%ltides) then
                     call check( nf90_put_var(iu%ncid, iu%k2_varid, pl%k2(j), start=[idslot, tslot]) )
                     call check( nf90_put_var(iu%ncid, iu%Q_varid, pl%Q(j), start=[idslot, tslot]) )
                  end if
               end select
            end do
         end associate
      class is (swiftest_cb)
         idslot = self%id + 1
         call check( nf90_put_var(iu%ncid, iu%id_varid, self%id, start=[idslot]) )

         charstring = trim(adjustl(self%info%name))
         strlen = len(charstring)
         call check( nf90_put_var(iu%ncid, iu%name_varid, charstring, start=[1, idslot], count=[strlen, 1]) )

         charstring = trim(adjustl(self%info%particle_type))
         strlen = len(charstring)
         call check( nf90_put_var(iu%ncid, iu%ptype_varid, charstring, start=[1, idslot], count=[strlen, 1]) )

         call check( nf90_put_var(iu%ncid, iu%Gmass_varid, self%Gmass, start=[idslot, tslot]) )
         call check( nf90_put_var(iu%ncid, iu%radius_varid, self%radius, start=[idslot, tslot]) )
         if (param%lrotation) then
            call check( nf90_put_var(iu%ncid, iu%Ip1_varid, self%Ip(1), start=[idslot, tslot]) )
            call check( nf90_put_var(iu%ncid, iu%Ip2_varid, self%Ip(2), start=[idslot, tslot]) )
            call check( nf90_put_var(iu%ncid, iu%Ip3_varid, self%Ip(3), start=[idslot, tslot]) )
            call check( nf90_put_var(iu%ncid, iu%rotx_varid, self%rot(1), start=[idslot, tslot]) )
            call check( nf90_put_var(iu%ncid, iu%roty_varid, self%rot(2), start=[idslot, tslot]) )
            call check( nf90_put_var(iu%ncid, iu%rotz_varid, self%rot(3), start=[idslot, tslot]) )
         end if
         if (param%ltides) then
            call check( nf90_put_var(iu%ncid, iu%k2_varid, self%k2, start=[idslot, tslot]) )
            call check( nf90_put_var(iu%ncid, iu%Q_varid, self%Q, start=[idslot, tslot]) )
         end if
      end select

      return
   end subroutine netcdf_write_frame_base

   module subroutine netcdf_write_hdr_system(self, iu, param) 
      !! author: David A. Minton
      !!
      !! Writes header information (variables that change with time, but not particle id). 
      !! This subroutine significantly improves the output over the original binary file, allowing us to track energy, momentum, and other quantities that 
      !! previously were handled as separate output files.
      implicit none
      ! Arguments
      class(swiftest_nbody_system), intent(in)    :: self  !! Swiftest nbody system object
      class(netcdf_parameters),     intent(inout) :: iu    !! Parameters used to for writing a NetCDF dataset to file
      class(swiftest_parameters),   intent(in)    :: param !! Current run configuration parameters
      ! Internals
      integer(I4B) :: tslot, old_mode

      tslot = int(param%ioutput, kind=I4B) + 1

      call check( nf90_open(param%outfile, nf90_write, iu%ncid) )
      call check( nf90_set_fill(iu%ncid, nf90_nofill, old_mode) )

      call check( nf90_put_var(iu%ncid, iu%time_varid, param%t, start=[tslot]) )
      call check( nf90_put_var(iu%ncid, iu%npl_varid, self%pl%nbody, start=[tslot]) )
      call check( nf90_put_var(iu%ncid, iu%ntp_varid, self%tp%nbody, start=[tslot]) )

      if (param%lenergy) then
         call check( nf90_put_var(iu%ncid, iu%KE_orb_varid, self%ke_orbit, start=[tslot]) )
         call check( nf90_put_var(iu%ncid, iu%KE_spin_varid, self%ke_spin, start=[tslot]) )
         call check( nf90_put_var(iu%ncid, iu%PE_varid, self%pe, start=[tslot]) )
         call check( nf90_put_var(iu%ncid, iu%L_orbx_varid, self%Lorbit(1), start=[tslot]) )
         call check( nf90_put_var(iu%ncid, iu%L_orby_varid, self%Lorbit(2), start=[tslot]) )
         call check( nf90_put_var(iu%ncid, iu%L_orbz_varid, self%Lorbit(3), start=[tslot]) )
         call check( nf90_put_var(iu%ncid, iu%L_spinx_varid, self%Lspin(1), start=[tslot]) )
         call check( nf90_put_var(iu%ncid, iu%L_spiny_varid, self%Lspin(2), start=[tslot]) )
         call check( nf90_put_var(iu%ncid, iu%L_spinz_varid, self%Lspin(3), start=[tslot]) )
         call check( nf90_put_var(iu%ncid, iu%L_escapex_varid, param%Lescape(1), start=[tslot]) )
         call check( nf90_put_var(iu%ncid, iu%L_escapey_varid, param%Lescape(2), start=[tslot]) )
         call check( nf90_put_var(iu%ncid, iu%L_escapez_varid, param%Lescape(3), start=[tslot]) )
         call check( nf90_put_var(iu%ncid, iu%Ecollisions_varid, param%Ecollisions, start=[tslot]) )
         call check( nf90_put_var(iu%ncid, iu%Euntracked_varid, param%Euntracked, start=[tslot]) )
         call check( nf90_put_var(iu%ncid, iu%GMescape_varid, param%GMescape, start=[tslot]) )
      end if

      return
   end subroutine netcdf_write_hdr_system


end submodule s_netcdf