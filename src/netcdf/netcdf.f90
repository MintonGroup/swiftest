!! Copyright 2022 - David Minton, Carlisle Wishard, Jennifer Pouplin, Jake Elliott, & Dana Singh
!! This file is part of Swiftest.
!! Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
!! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!! Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
!! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!! You should have received a copy of the GNU General Public License along with Swiftest. 
!! If not, see: https://www.gnu.org/licenses. 

submodule (swiftest_classes) s_netcdf
   use swiftest
   use netcdf
contains

   subroutine check(status, call_identifier)
      !! author: Carlisle A. Wishard, Dana Singh, and David A. Minton
      !!
      !! Checks the status of all NetCDF operations to catch errors
      implicit none
      ! Arguments
      integer, intent (in) :: status !! The status code returned by a NetCDF function
      character(len=*), intent(in), optional :: call_identifier !! String that indicates which calling function caused the error for diagnostic purposes

      if(status /= nf90_noerr) then
         if (present(call_identifier)) write(*,*) "NetCDF error in ",trim(call_identifier)
         write(*,*) trim(nf90_strerror(status))
         call util_exit(FAILURE)
      end if

      return
   end subroutine check


   module subroutine netcdf_close(self)
      !! author: Carlisle A. Wishard, Dana Singh, and David A. Minton
      !!
      !! Closes a NetCDF file
      implicit none
      ! Arguments
      class(netcdf_parameters),   intent(inout) :: self   !! Parameters used to identify a particular NetCDF dataset

      call check( nf90_close(self%ncid), "netcdf_close" )

      return
   end subroutine netcdf_close


   module subroutine netcdf_flush(self, param)
      !! author: David A. Minton
      !!
      !! Flushes the current buffer to disk by closing and re-opening the file.
      !!    
      implicit none
      ! Arguments
      class(netcdf_parameters),   intent(inout) :: self !! Parameters used to identify a particular NetCDF dataset
      class(swiftest_parameters), intent(inout) :: param !! Current run configuration parameters 

      call self%close()
      call self%open(param)

      return
   end subroutine netcdf_flush


   module function netcdf_get_old_t_final_system(self, param) result(old_t_final)
      !! author: David A. Minton
      !!
      !! Validates the dump file to check whether the dump file initial conditions duplicate the last frame of the netcdf output.
      !!
      implicit none
      ! Arguments
      class(swiftest_nbody_system), intent(inout) :: self
      class(swiftest_parameters),   intent(inout) :: param
      ! Result
      real(DP)                                    :: old_t_final
      ! Internals
      integer(I4B)                              :: itmax, idmax
      real(DP), dimension(:), allocatable       :: vals
      real(DP), dimension(1)                    :: val
      real(DP), dimension(NDIM)                 :: rot0, Ip0, Lnow
      real(DP) :: KE_orb_orig, KE_spin_orig, PE_orig

      call param%nciu%open(param)
      call check( nf90_inquire_dimension(param%nciu%ncid, param%nciu%time_dimid, len=itmax), "netcdf_get_old_t_final_system time_dimid" )
      call check( nf90_inquire_dimension(param%nciu%ncid, param%nciu%id_dimid, len=idmax), "netcdf_get_old_t_final_system id_dimid" )
      allocate(vals(idmax))
      call check( nf90_get_var(param%nciu%ncid, param%nciu%time_varid, val, start=[1], count=[1]), "netcdf_get_old_t_final_system time_varid" )

      !old_t_final = val(1)
      old_t_final = param%t0 ! For NetCDF it is safe to overwrite the final t value on a restart

      if (param%lenergy) then
         call check( nf90_get_var(param%nciu%ncid, param%nciu%KE_orb_varid, val, start=[1], count=[1]), "netcdf_get_old_t_final_system KE_orb_varid" )
         KE_orb_orig = val(1)

         call check( nf90_get_var(param%nciu%ncid, param%nciu%KE_spin_varid, val, start=[1], count=[1]), "netcdf_get_old_t_final_system KE_spin_varid" )
         KE_spin_orig = val(1)

         call check( nf90_get_var(param%nciu%ncid, param%nciu%PE_varid, val, start=[1], count=[1]), "netcdf_get_old_t_final_system PE_varid" )
         PE_orig = val(1)

         call check( nf90_get_var(param%nciu%ncid, param%nciu%Ecollisions_varid, self%Ecollisions, start=[1]), "netcdf_get_old_t_final_system Ecollisions_varid" )
         call check( nf90_get_var(param%nciu%ncid, param%nciu%Euntracked_varid,  self%Euntracked,  start=[1]), "netcdf_get_old_t_final_system Euntracked_varid" )

         self%Eorbit_orig = KE_orb_orig + KE_spin_orig + PE_orig + self%Ecollisions + self%Euntracked

         call check( nf90_get_var(param%nciu%ncid, param%nciu%L_orbx_varid, val, start=[1], count=[1]), "netcdf_get_old_t_final_system L_orbx_varid" )
         self%Lorbit_orig(1) = val(1)
         call check( nf90_get_var(param%nciu%ncid, param%nciu%L_orby_varid, val, start=[1], count=[1]), "netcdf_get_old_t_final_system L_orby_varid" )
         self%Lorbit_orig(2) = val(1)
         call check( nf90_get_var(param%nciu%ncid, param%nciu%L_orbz_varid, val, start=[1], count=[1]), "netcdf_get_old_t_final_system L_orbz_varid" )
         self%Lorbit_orig(3) = val(1)

         call check( nf90_get_var(param%nciu%ncid, param%nciu%L_spinx_varid, val, start=[1], count=[1]), "netcdf_get_old_t_final_system L_spinx_varid" )
         self%Lspin_orig(1) = val(1)
         call check( nf90_get_var(param%nciu%ncid, param%nciu%L_spiny_varid, val, start=[1], count=[1]), "netcdf_get_old_t_final_system L_spiny_varid" )
         self%Lspin_orig(2) = val(1)
         call check( nf90_get_var(param%nciu%ncid, param%nciu%L_spinz_varid, val, start=[1], count=[1]), "netcdf_get_old_t_final_system L_spinz_varid" )
         self%Lspin_orig(3) = val(1)

         call check( nf90_get_var(param%nciu%ncid, param%nciu%L_escapex_varid, self%Lescape(1),  start=[1]), "netcdf_get_old_t_final_system L_escapex_varid" )
         call check( nf90_get_var(param%nciu%ncid, param%nciu%L_escapey_varid, self%Lescape(2),  start=[1]), "netcdf_get_old_t_final_system L_escapey_varid" )
         call check( nf90_get_var(param%nciu%ncid, param%nciu%L_escapez_varid, self%Lescape(3),  start=[1]), "netcdf_get_old_t_final_system L_escapez_varid" )

         self%Ltot_orig(:) = self%Lorbit_orig(:) + self%Lspin_orig(:) + self%Lescape(:)

         call check( nf90_get_var(param%nciu%ncid, param%nciu%Gmass_varid, vals, start=[1,1], count=[idmax,1]), "netcdf_get_old_t_final_system Gmass_varid" )
         call check( nf90_get_var(param%nciu%ncid, param%nciu%GMescape_varid,    self%GMescape,    start=[1]), "netcdf_get_old_t_final_system GMescape_varid" )
         self%GMtot_orig = vals(1) + sum(vals(2:idmax), vals(2:idmax) == vals(2:idmax)) + self%GMescape

         select type(cb => self%cb)
         class is (symba_cb)
            cb%GM0 = vals(1)
            cb%dGM = cb%Gmass - cb%GM0

            call check( nf90_get_var(param%nciu%ncid, param%nciu%radius_varid, val, start=[1,1], count=[1,1]), "netcdf_get_old_t_final_system radius_varid" )
            cb%R0 = val(1) 

            if (param%lrotation) then

               call check( nf90_get_var(param%nciu%ncid, param%nciu%rotx_varid, val, start=[1,1], count=[1,1]), "netcdf_get_old_t_final_system rotx_varid" )
               rot0(1) = val(1)
               call check( nf90_get_var(param%nciu%ncid, param%nciu%roty_varid, val, start=[1,1], count=[1,1]), "netcdf_get_old_t_final_system roty_varid" )
               rot0(2) = val(1)
               call check( nf90_get_var(param%nciu%ncid, param%nciu%rotz_varid, val, start=[1,1], count=[1,1]), "netcdf_get_old_t_final_system rotz_varid" )
               rot0(3) = val(1)

               call check( nf90_get_var(param%nciu%ncid, param%nciu%Ip1_varid, val, start=[1,1], count=[1,1]), "netcdf_get_old_t_final_system Ip1_varid" )
               Ip0(1) = val(1)
               call check( nf90_get_var(param%nciu%ncid, param%nciu%Ip2_varid, val, start=[1,1], count=[1,1]), "netcdf_get_old_t_final_system Ip2_varid" )
               Ip0(2) = val(1)
               call check( nf90_get_var(param%nciu%ncid, param%nciu%Ip3_varid, val, start=[1,1], count=[1,1]), "netcdf_get_old_t_final_system Ip3_varid" )
               Ip0(3) = val(1)

               cb%L0(:) = Ip0(3) * cb%GM0 * cb%R0**2 * rot0(:)

               Lnow(:) = cb%Ip(3) * cb%Gmass * cb%radius**2 * cb%rot(:)
               cb%dL(:) = Lnow(:) - cb%L0(:)
            end if
         end select

      end if

      deallocate(vals)
      
      return
   end function netcdf_get_old_t_final_system


   module subroutine netcdf_initialize_output(self, param)
      !! author: Carlisle A. Wishard, Dana Singh, and David A. Minton
      !!
      !! Initialize a NetCDF file system and defines all variables.
      use, intrinsic :: ieee_arithmetic
      implicit none
      ! Arguments
      class(netcdf_parameters),   intent(inout) :: self    !! Parameters used to identify a particular NetCDF dataset
      class(swiftest_parameters), intent(in)    :: param           !! Current run configuration parameters 
      ! Internals
      integer(I4B) :: nvar, varid, vartype
      real(DP) :: dfill
      real(SP) :: sfill
      logical :: fileExists
      character(len=STRMAX) :: errmsg
      integer(I4B) :: ndims

      dfill = ieee_value(dfill, IEEE_QUIET_NAN)
      sfill = ieee_value(sfill, IEEE_QUIET_NAN)

      ! Check if the file exists, and if it does, delete it
      inquire(file=param%outfile, exist=fileExists)
      if (fileExists) then
         open(unit=LUN, file=param%outfile, status="old", err=667, iomsg=errmsg)
         close(unit=LUN, status="delete")
      end if

      call check( nf90_create(param%outfile, NF90_NETCDF4, self%ncid), "netcdf_initialize_output nf90_create" )

      ! Define the NetCDF dimensions with particle name as the record dimension
      call check( nf90_def_dim(self%ncid, ID_DIMNAME, NF90_UNLIMITED, self%id_dimid), "netcdf_initialize_output nf90_def_dim id_dimid" )     ! 'x' dimension
      call check( nf90_def_dim(self%ncid, STR_DIMNAME, NAMELEN, self%str_dimid), "netcdf_initialize_output nf90_def_dim str_dimid"  ) ! Dimension for string variables (aka character arrays)
      call check( nf90_def_dim(self%ncid, TIME_DIMNAME, NF90_UNLIMITED, self%time_dimid), "netcdf_initialize_output nf90_def_dim time_dimid" ) ! 'y' dimension

      select case (param%out_type)
      case(NETCDF_FLOAT_TYPE)
         self%out_type = NF90_FLOAT
      case(NETCDF_DOUBLE_TYPE)
         self%out_type = NF90_DOUBLE
      end select

      !! Define the variables
      call check( nf90_def_var(self%ncid, TIME_DIMNAME, self%out_type, self%time_dimid, self%time_varid), "netcdf_initialize_output nf90_def_var time_varid"  )
      call check( nf90_def_var(self%ncid, ID_DIMNAME, NF90_INT, self%id_dimid, self%id_varid), "netcdf_initialize_output nf90_def_var id_varid"  )
      call check( nf90_def_var(self%ncid, NPL_VARNAME, NF90_INT, self%time_dimid, self%npl_varid), "netcdf_initialize_output nf90_def_var npl_varid"  )
      call check( nf90_def_var(self%ncid, NTP_VARNAME, NF90_INT, self%time_dimid, self%ntp_varid), "netcdf_initialize_output nf90_def_var ntp_varid"  )
      if (param%integrator == SYMBA) call check( nf90_def_var(self%ncid, NPLM_VARNAME, NF90_INT, self%time_dimid, self%nplm_varid), "netcdf_initialize_output nf90_def_var nplm_varid"  )
      call check( nf90_def_var(self%ncid, NAME_VARNAME, NF90_CHAR, [self%str_dimid, self%id_dimid], self%name_varid), "netcdf_initialize_output nf90_def_var name_varid"  )
      call check( nf90_def_var(self%ncid, PTYPE_VARNAME, NF90_CHAR, [self%str_dimid, self%id_dimid], self%ptype_varid), "netcdf_initialize_output nf90_def_var ptype_varid"  )
      call check( nf90_def_var(self%ncid, STATUS_VARNAME, NF90_CHAR, [self%str_dimid, self%id_dimid], self%status_varid), "netcdf_initialize_output nf90_def_var status_varid"  )

      if ((param%out_form == XV) .or. (param%out_form == XVEL)) then
         call check( nf90_def_var(self%ncid, XHX_VARNAME, self%out_type, [self%id_dimid, self%time_dimid], self%xhx_varid), "netcdf_initialize_output nf90_def_var xhx_varid"  )
         call check( nf90_def_var(self%ncid, XHY_VARNAME, self%out_type, [self%id_dimid, self%time_dimid], self%xhy_varid), "netcdf_initialize_output nf90_def_var xhy_varid"  )
         call check( nf90_def_var(self%ncid, XHZ_VARNAME, self%out_type, [self%id_dimid, self%time_dimid], self%xhz_varid), "netcdf_initialize_output nf90_def_var xhz_varid"  )
         call check( nf90_def_var(self%ncid, VHX_VARNAME, self%out_type, [self%id_dimid, self%time_dimid], self%vhx_varid), "netcdf_initialize_output nf90_def_var vhx_varid"  )
         call check( nf90_def_var(self%ncid, VHY_VARNAME, self%out_type, [self%id_dimid, self%time_dimid], self%vhy_varid), "netcdf_initialize_output nf90_def_var vhy_varid"  )
         call check( nf90_def_var(self%ncid, VHZ_VARNAME, self%out_type, [self%id_dimid, self%time_dimid], self%vhz_varid), "netcdf_initialize_output nf90_def_var vhz_varid"  )

         !! When GR is enabled, we need to save the pseudovelocity vectors in addition to the true heliocentric velocity vectors, otherwise
         !! we cannnot expect bit-identical runs from restarted runs with GR enabled due to floating point errors during the conversion.
         if (param%lgr) then
            call check( nf90_def_var(self%ncid, GR_PSEUDO_VHX_VARNAME, self%out_type, [self%id_dimid, self%time_dimid], self%gr_pseudo_vhx_varid), "netcdf_initialize_output nf90_def_var gr_psuedo_vhx_varid"  )
            call check( nf90_def_var(self%ncid, GR_PSEUDO_VHY_VARNAME, self%out_type, [self%id_dimid, self%time_dimid], self%gr_pseudo_vhy_varid), "netcdf_initialize_output nf90_def_var gr_psuedo_vhy_varid"  )
            call check( nf90_def_var(self%ncid, GR_PSEUDO_VHZ_VARNAME, self%out_type, [self%id_dimid, self%time_dimid], self%gr_pseudo_vhz_varid), "netcdf_initialize_output nf90_def_var gr_psuedo_vhz_varid"  )
            self%lpseudo_vel_exists = .true.
         end if

      end if
   
      if ((param%out_form == EL) .or. (param%out_form == XVEL)) then
         call check( nf90_def_var(self%ncid, A_VARNAME, self%out_type, [self%id_dimid, self%time_dimid], self%a_varid), "netcdf_initialize_output nf90_def_var a_varid"  )
         call check( nf90_def_var(self%ncid, E_VARNAME, self%out_type, [self%id_dimid, self%time_dimid], self%e_varid), "netcdf_initialize_output nf90_def_var e_varid"  )
         call check( nf90_def_var(self%ncid, INC_VARNAME, self%out_type, [self%id_dimid, self%time_dimid], self%inc_varid), "netcdf_initialize_output nf90_def_var inc_varid"  )
         call check( nf90_def_var(self%ncid, CAPOM_VARNAME, self%out_type, [self%id_dimid, self%time_dimid], self%capom_varid), "netcdf_initialize_output nf90_def_var capom_varid"  )
         call check( nf90_def_var(self%ncid, OMEGA_VARNAME, self%out_type, [self%id_dimid, self%time_dimid], self%omega_varid), "netcdf_initialize_output nf90_def_var omega_varid"  )
         call check( nf90_def_var(self%ncid, CAPM_VARNAME, self%out_type, [self%id_dimid, self%time_dimid], self%capm_varid), "netcdf_initialize_output nf90_def_var capm_varid"  )
      end if

      call check( nf90_def_var(self%ncid, GMASS_VARNAME, self%out_type, [self%id_dimid, self%time_dimid], self%Gmass_varid), "netcdf_initialize_output nf90_def_var Gmass_varid"  )

      if (param%lrhill_present) then
         call check( nf90_def_var(self%ncid, RHILL_VARNAME, self%out_type, [self%id_dimid, self%time_dimid], self%rhill_varid), "netcdf_initialize_output nf90_def_var rhill_varid"  )
      end if

      if (param%lclose) then
         call check( nf90_def_var(self%ncid, RADIUS_VARNAME, self%out_type, [self%id_dimid, self%time_dimid], self%radius_varid), "netcdf_initialize_output nf90_def_var radius_varid"  )

         call check( nf90_def_var(self%ncid, ORIGIN_TIME_VARNAME, self%out_type, self%id_dimid, self%origin_time_varid), "netcdf_initialize_output nf90_def_var origin_time_varid"  )
         call check( nf90_def_var(self%ncid, ORIGIN_TYPE_VARNAME, NF90_CHAR, [self%str_dimid, self%id_dimid], &
                                  self%origin_type_varid), "netcdf_initialize_output nf90_create"  )
         call check( nf90_def_var(self%ncid, ORIGIN_XHX_VARNAME, self%out_type, self%id_dimid, self%origin_xhx_varid), "netcdf_initialize_output nf90_def_var origin_xhx_varid"  )
         call check( nf90_def_var(self%ncid, ORIGIN_XHY_VARNAME, self%out_type, self%id_dimid, self%origin_xhy_varid), "netcdf_initialize_output nf90_def_var origin_xhy_varid"  )
         call check( nf90_def_var(self%ncid, ORIGIN_XHZ_VARNAME, self%out_type, self%id_dimid, self%origin_xhz_varid), "netcdf_initialize_output nf90_def_var origin_xhz_varid"  )
         call check( nf90_def_var(self%ncid, ORIGIN_VHX_VARNAME, self%out_type, self%id_dimid, self%origin_vhx_varid), "netcdf_initialize_output nf90_def_var origin_vhx_varid"  )
         call check( nf90_def_var(self%ncid, ORIGIN_VHY_VARNAME, self%out_type, self%id_dimid, self%origin_vhy_varid), "netcdf_initialize_output nf90_def_var origin_vhy_varid"  )
         call check( nf90_def_var(self%ncid, ORIGIN_VHZ_VARNAME, self%out_type, self%id_dimid,  self%origin_vhz_varid), "netcdf_initialize_output nf90_def_var origin_vhz_varid"  )

         call check( nf90_def_var(self%ncid, COLLISION_ID_VARNAME, NF90_INT, self%id_dimid, self%collision_id_varid), "netcdf_initialize_output nf90_def_var collision_id_varid"  )
         call check( nf90_def_var(self%ncid, DISCARD_TIME_VARNAME, self%out_type, self%id_dimid, self%discard_time_varid), "netcdf_initialize_output nf90_def_var discard_time_varid"  )
         call check( nf90_def_var(self%ncid, DISCARD_XHX_VARNAME, self%out_type, self%id_dimid, self%discard_xhx_varid), "netcdf_initialize_output nf90_def_var discard_xhx_varid"  )
         call check( nf90_def_var(self%ncid, DISCARD_XHY_VARNAME, self%out_type, self%id_dimid, self%discard_xhy_varid), "netcdf_initialize_output nf90_def_var discard_xhy_varid"  )
         call check( nf90_def_var(self%ncid, DISCARD_XHZ_VARNAME, self%out_type, self%id_dimid, self%discard_xhz_varid), "netcdf_initialize_output nf90_def_var discard_xhz_varid"  )
         call check( nf90_def_var(self%ncid, DISCARD_VHX_VARNAME, self%out_type, self%id_dimid, self%discard_vhx_varid), "netcdf_initialize_output nf90_def_var discard_vhx_varid"  )
         call check( nf90_def_var(self%ncid, DISCARD_VHY_VARNAME, self%out_type, self%id_dimid, self%discard_vhy_varid), "netcdf_initialize_output nf90_def_var discard_vhy_varid"  )
         call check( nf90_def_var(self%ncid, DISCARD_VHZ_VARNAME, self%out_type, self%id_dimid,  self%discard_vhz_varid), "netcdf_initialize_output nf90_def_var discard_vhz_varid"  )
         call check( nf90_def_var(self%ncid, DISCARD_BODY_ID_VARNAME, NF90_INT, self%id_dimid, self%discard_body_id_varid), "netcdf_initialize_output nf90_def_var discard_body_id_varid"  )
      end if

      if (param%lrotation) then
         call check( nf90_def_var(self%ncid, IP1_VARNAME, self%out_type, [self%id_dimid, self%time_dimid], self%Ip1_varid), "netcdf_initialize_output nf90_def_var Ip1_varid"  )
         call check( nf90_def_var(self%ncid, IP2_VARNAME, self%out_type, [self%id_dimid, self%time_dimid], self%Ip2_varid), "netcdf_initialize_output nf90_def_var Ip2_varid"  )
         call check( nf90_def_var(self%ncid, IP3_VARNAME, self%out_type, [self%id_dimid, self%time_dimid], self%Ip3_varid), "netcdf_initialize_output nf90_def_var Ip3_varid"  )
         call check( nf90_def_var(self%ncid, ROTX_VARNAME, self%out_type, [self%id_dimid, self%time_dimid], self%rotx_varid), "netcdf_initialize_output nf90_def_var rotx_varid"  )
         call check( nf90_def_var(self%ncid, ROTY_VARNAME, self%out_type, [self%id_dimid, self%time_dimid], self%roty_varid), "netcdf_initialize_output nf90_def_var roty_varid"  )
         call check( nf90_def_var(self%ncid, ROTZ_VARNAME, self%out_type, [self%id_dimid, self%time_dimid], self%rotz_varid), "netcdf_initialize_output nf90_def_var rotz_varid"  )
      end if

      ! if (param%ltides) then
      !    call check( nf90_def_var(self%ncid, K2_VARNAME, self%out_type, [self%id_dimid, self%time_dimid], self%k2_varid), "netcdf_initialize_output nf90_def_var k2_varid"  )
      !    call check( nf90_def_var(self%ncid, Q_VARNAME, self%out_type, [self%id_dimid, self%time_dimid], self%Q_varid), "netcdf_initialize_output nf90_def_var Q_varid"  )
      ! end if

      if (param%lenergy) then
         call check( nf90_def_var(self%ncid, KE_ORB_VARNAME, self%out_type, self%time_dimid, self%KE_orb_varid), "netcdf_initialize_output nf90_def_var KE_orb_varid"  )
         call check( nf90_def_var(self%ncid, KE_SPIN_VARNAME, self%out_type, self%time_dimid, self%KE_spin_varid), "netcdf_initialize_output nf90_def_var KE_spin_varid"  )
         call check( nf90_def_var(self%ncid, PE_VARNAME, self%out_type, self%time_dimid, self%PE_varid), "netcdf_initialize_output nf90_def_var PE_varid"  )
         call check( nf90_def_var(self%ncid, L_ORBX_VARNAME, self%out_type, self%time_dimid, self%L_orbx_varid), "netcdf_initialize_output nf90_def_var L_orbx_varid"  )
         call check( nf90_def_var(self%ncid, L_ORBY_VARNAME, self%out_type, self%time_dimid, self%L_orby_varid), "netcdf_initialize_output nf90_def_var L_orby_varid"  )
         call check( nf90_def_var(self%ncid, L_ORBZ_VARNAME, self%out_type, self%time_dimid, self%L_orbz_varid), "netcdf_initialize_output nf90_def_var L_orbz_varid"  )
         call check( nf90_def_var(self%ncid, L_SPINX_VARNAME, self%out_type, self%time_dimid, self%L_spinx_varid), "netcdf_initialize_output nf90_def_var L_spinx_varid"  )
         call check( nf90_def_var(self%ncid, L_SPINY_VARNAME, self%out_type, self%time_dimid, self%L_spiny_varid), "netcdf_initialize_output nf90_def_var L_spiny_varid"  )
         call check( nf90_def_var(self%ncid, L_SPINZ_VARNAME, self%out_type, self%time_dimid, self%L_spinz_varid), "netcdf_initialize_output nf90_def_var L_spinz_varid"  )
         call check( nf90_def_var(self%ncid, L_ESCAPEX_VARNAME, self%out_type, self%time_dimid, self%L_escapex_varid), "netcdf_initialize_output nf90_def_var L_escapex_varid"  )
         call check( nf90_def_var(self%ncid, L_ESCAPEY_VARNAME, self%out_type, self%time_dimid, self%L_escapey_varid), "netcdf_initialize_output nf90_def_var L_escapey_varid"  )
         call check( nf90_def_var(self%ncid, L_ESCAPEZ_VARNAME, self%out_type, self%time_dimid, self%L_escapez_varid), "netcdf_initialize_output nf90_def_var L_escapez_varid"  )
         call check( nf90_def_var(self%ncid, ECOLLISIONS_VARNAME, self%out_type, self%time_dimid, self%Ecollisions_varid), "netcdf_initialize_output nf90_def_var Ecollisions_varid"  )
         call check( nf90_def_var(self%ncid, EUNTRACKED_VARNAME, self%out_type, self%time_dimid, self%Euntracked_varid), "netcdf_initialize_output nf90_def_var Euntracked_varid"  )
         call check( nf90_def_var(self%ncid, GMESCAPE_VARNAME, self%out_type, self%time_dimid, self%GMescape_varid), "netcdf_initialize_output nf90_def_var GMescape_varid"  )
      end if

      call check( nf90_def_var(self%ncid, J2RP2_VARNAME, self%out_type, self%time_dimid, self%j2rp2_varid), "netcdf_initialize_output nf90_def_var j2rp2_varid"  )
      call check( nf90_def_var(self%ncid, J4RP4_VARNAME, self%out_type, self%time_dimid, self%j4rp4_varid), "netcdf_initialize_output nf90_def_var j4rp4_varid"  )


      ! Set fill mode to NaN for all variables
      call check( nf90_inquire(self%ncid, nVariables=nvar), "netcdf_initialize_output nf90_inquire nVariables"  )
      do varid = 1, nvar
         call check( nf90_inquire_variable(self%ncid, varid, xtype=vartype, ndims=ndims), "netcdf_initialize_output nf90_inquire_variable"  )
         select case(vartype)
         case(NF90_INT)
            call check( nf90_def_var_fill(self%ncid, varid, 0, NF90_FILL_INT), "netcdf_initialize_output nf90_def_var_fill NF90_INT"  )
         case(NF90_FLOAT)
            call check( nf90_def_var_fill(self%ncid, varid, 0, sfill), "netcdf_initialize_output nf90_def_var_fill NF90_FLOAT"  )
         case(NF90_DOUBLE)
            call check( nf90_def_var_fill(self%ncid, varid, 0, dfill), "netcdf_initialize_output nf90_def_var_fill NF90_DOUBLE"  )
         case(NF90_CHAR)
            call check( nf90_def_var_fill(self%ncid, varid, 0, 0), "netcdf_initialize_output nf90_def_var_fill NF90_CHAR"  )
         end select
      end do

      ! Take the file out of define mode
      call check( nf90_enddef(self%ncid), "netcdf_initialize_output nf90_enddef"  )

      return

      667 continue
      write(*,*) "Error creating NetCDF output file. " // trim(adjustl(errmsg))
      call util_exit(FAILURE)
   end subroutine netcdf_initialize_output


   module subroutine netcdf_open(self, param, readonly)
      !! author: Carlisle A. Wishard, Dana Singh, and David A. Minton
      !!
      !! Opens a NetCDF file and does the variable inquiries to activate variable ids
      implicit none
      ! Arguments
      class(netcdf_parameters),   intent(inout) :: self     !! Parameters used to identify a particular NetCDF dataset
      class(swiftest_parameters), intent(in)    :: param    !! Current run configuration parameters
      logical, optional,          intent(in)    :: readonly !! Logical flag indicating that this should be open read only
      ! Internals
      integer(I4B) :: mode, status
      character(len=NF90_MAX_NAME) :: str_dim_name
      character(len=STRMAX) :: errmsg

      mode = NF90_WRITE
      if (present(readonly)) then
         if (readonly) mode = NF90_NOWRITE
      end if

      write(errmsg,*) "netcdf_open nf90_open ",trim(adjustl(param%outfile))
      call check( nf90_open(param%outfile, mode, self%ncid), errmsg)

      call check( nf90_inq_dimid(self%ncid, TIME_DIMNAME, self%time_dimid), "netcdf_open nf90_inq_dimid time_dimid"  )
      call check( nf90_inq_dimid(self%ncid, ID_DIMNAME, self%id_dimid), "netcdf_open nf90_inq_dimid id_dimid"  )
      if (max(self%time_dimid,self%id_dimid) == 2) then
         self%str_dimid = 3
      else if (min(self%time_dimid,self%id_dimid) == 0) then
         self%str_dimid = 1
      else
         self%str_dimid = 2
      end if 
      call check( nf90_inquire_dimension(self%ncid, self%str_dimid, name=str_dim_name), "netcdf_open nf90_inquire_dimension str_dim_name"  )
      call check( nf90_inq_dimid(self%ncid, str_dim_name, self%str_dimid), "netcdf_open nf90_inq_dimid str_dimid"  )

      ! Required Variables

      call check( nf90_inq_varid(self%ncid, TIME_DIMNAME, self%time_varid), "netcdf_open nf90_inq_varid time_varid" )
      call check( nf90_inq_varid(self%ncid, ID_DIMNAME, self%id_varid), "netcdf_open nf90_inq_varid id_varid" )
      call check( nf90_inq_varid(self%ncid, NAME_VARNAME, self%name_varid), "netcdf_open nf90_inq_varid name_varid" )
      call check( nf90_inq_varid(self%ncid, PTYPE_VARNAME, self%ptype_varid), "netcdf_open nf90_inq_varid ptype_varid" )
      call check( nf90_inq_varid(self%ncid, GMASS_VARNAME, self%Gmass_varid), "netcdf_open nf90_inq_varid Gmass_varid" )

      if ((param%out_form == XV) .or. (param%out_form == XVEL)) then
         call check( nf90_inq_varid(self%ncid, XHX_VARNAME, self%xhx_varid), "netcdf_open nf90_inq_varid xhx_varid" )
         call check( nf90_inq_varid(self%ncid, XHY_VARNAME, self%xhy_varid), "netcdf_open nf90_inq_varid xhy_varid" )
         call check( nf90_inq_varid(self%ncid, XHZ_VARNAME, self%xhz_varid), "netcdf_open nf90_inq_varid xhz_varid" )
         call check( nf90_inq_varid(self%ncid, VHX_VARNAME, self%vhx_varid), "netcdf_open nf90_inq_varid vhx_varid" )
         call check( nf90_inq_varid(self%ncid, VHY_VARNAME, self%vhy_varid), "netcdf_open nf90_inq_varid vhy_varid" )
         call check( nf90_inq_varid(self%ncid, VHZ_VARNAME, self%vhz_varid), "netcdf_open nf90_inq_varid vhz_varid" )

         if (param%lgr) then
            !! check if pseudovelocity vectors exist in this file. If they are, set the correct flag so we know whe should not do the conversion.
            status = nf90_inq_varid(self%ncid, GR_PSEUDO_VHX_VARNAME, self%gr_pseudo_vhx_varid)
            self%lpseudo_vel_exists = (status == nf90_noerr)
            if (self%lpseudo_vel_exists) then
               status = nf90_inq_varid(self%ncid, GR_PSEUDO_VHY_VARNAME, self%gr_pseudo_vhy_varid)
               self%lpseudo_vel_exists = (status == nf90_noerr)
               if (self%lpseudo_vel_exists) then
                  status = nf90_inq_varid(self%ncid, GR_PSEUDO_VHZ_VARNAME, self%gr_pseudo_vhz_varid)
                  self%lpseudo_vel_exists = (status == nf90_noerr)
               end if
            end if
            if (.not.self%lpseudo_vel_exists) then
               write(*,*) "Warning! Pseudovelocity not found in input file for GR enabled run. If this is a restarted run, bit-identical trajectories are not guarunteed!"
            end if

         end if
      end if

      if ((param%out_form == EL) .or. (param%out_form == XVEL)) then
         call check( nf90_inq_varid(self%ncid, A_VARNAME, self%a_varid), "netcdf_open nf90_inq_varid a_varid" )
         call check( nf90_inq_varid(self%ncid, E_VARNAME, self%e_varid), "netcdf_open nf90_inq_varid e_varid" )
         call check( nf90_inq_varid(self%ncid, INC_VARNAME, self%inc_varid), "netcdf_open nf90_inq_varid inc_varid" )
         call check( nf90_inq_varid(self%ncid, CAPOM_VARNAME, self%capom_varid), "netcdf_open nf90_inq_varid capom_varid" )
         call check( nf90_inq_varid(self%ncid, OMEGA_VARNAME, self%omega_varid), "netcdf_open nf90_inq_varid omega_varid" )
         call check( nf90_inq_varid(self%ncid, CAPM_VARNAME, self%capm_varid), "netcdf_open nf90_inq_varid capm_varid" )
      end if

      if (param%lclose) then
         call check( nf90_inq_varid(self%ncid, RADIUS_VARNAME, self%radius_varid), "netcdf_open nf90_inq_varid radius_varid" )
      end if 

      if (param%lrotation) then
         call check( nf90_inq_varid(self%ncid, IP1_VARNAME, self%Ip1_varid), "netcdf_open nf90_inq_varid Ip1_varid" )
         call check( nf90_inq_varid(self%ncid, IP2_VARNAME, self%Ip2_varid), "netcdf_open nf90_inq_varid Ip2_varid" )
         call check( nf90_inq_varid(self%ncid, IP3_VARNAME, self%Ip3_varid), "netcdf_open nf90_inq_varid Ip3_varid" )
         call check( nf90_inq_varid(self%ncid, ROTX_VARNAME, self%rotx_varid), "netcdf_open nf90_inq_varid rotx_varid" )
         call check( nf90_inq_varid(self%ncid, ROTY_VARNAME, self%roty_varid), "netcdf_open nf90_inq_varid roty_varid" )
         call check( nf90_inq_varid(self%ncid, ROTZ_VARNAME, self%rotz_varid), "netcdf_open nf90_inq_varid rotz_varid" )
      end if

      ! if (param%ltides) then
      !    call check( nf90_inq_varid(self%ncid, K2_VARNAME, self%k2_varid), "netcdf_open nf90_inq_varid k2_varid" )
      !    call check( nf90_inq_varid(self%ncid, Q_VARNAME, self%Q_varid), "netcdf_open nf90_inq_varid Q_varid" )
      ! end if

      ! Optional Variables
      if (param%lrhill_present) then
         status = nf90_inq_varid(self%ncid, RHILL_VARNAME, self%rhill_varid)
         if (status /= nf90_noerr) write(*,*) "Warning! RHILL variable not set in input file. Calculating."
      end if

      ! Optional variables The User Doesn't Need to Know About
      status = nf90_inq_varid(self%ncid, NPL_VARNAME, self%npl_varid)
      status = nf90_inq_varid(self%ncid, NTP_VARNAME, self%ntp_varid)
      status = nf90_inq_varid(self%ncid, STATUS_VARNAME, self%status_varid)
      status = nf90_inq_varid(self%ncid, J2RP2_VARNAME, self%j2rp2_varid)
      status = nf90_inq_varid(self%ncid, J4RP4_VARNAME, self%j4rp4_varid)

      if (param%integrator == SYMBA) then
         status = nf90_inq_varid(self%ncid, NPLM_VARNAME, self%nplm_varid)
      end if

      if (param%lclose) then
         status = nf90_inq_varid(self%ncid, ORIGIN_TYPE_VARNAME, self%origin_type_varid)
         status = nf90_inq_varid(self%ncid, ORIGIN_TIME_VARNAME, self%origin_time_varid)
         status = nf90_inq_varid(self%ncid, ORIGIN_XHX_VARNAME, self%origin_xhx_varid)
         status = nf90_inq_varid(self%ncid, ORIGIN_XHY_VARNAME, self%origin_xhy_varid)
         status = nf90_inq_varid(self%ncid, ORIGIN_XHZ_VARNAME, self%origin_xhz_varid)
         status = nf90_inq_varid(self%ncid, ORIGIN_VHX_VARNAME, self%origin_vhx_varid)
         status = nf90_inq_varid(self%ncid, ORIGIN_VHY_VARNAME, self%origin_vhy_varid)
         status = nf90_inq_varid(self%ncid, ORIGIN_VHZ_VARNAME, self%origin_vhz_varid)
         status = nf90_inq_varid(self%ncid, COLLISION_ID_VARNAME, self%collision_id_varid)
         status = nf90_inq_varid(self%ncid, DISCARD_TIME_VARNAME, self%discard_time_varid)
         status = nf90_inq_varid(self%ncid, DISCARD_XHX_VARNAME, self%discard_xhx_varid)
         status = nf90_inq_varid(self%ncid, DISCARD_XHY_VARNAME, self%discard_xhy_varid)
         status = nf90_inq_varid(self%ncid, DISCARD_XHZ_VARNAME, self%discard_xhz_varid)
         status = nf90_inq_varid(self%ncid, DISCARD_VHX_VARNAME, self%discard_vhx_varid)
         status = nf90_inq_varid(self%ncid, DISCARD_VHY_VARNAME, self%discard_vhy_varid)
         status = nf90_inq_varid(self%ncid, DISCARD_VHZ_VARNAME, self%discard_vhz_varid)
         status = nf90_inq_varid(self%ncid, DISCARD_BODY_ID_VARNAME, self%discard_body_id_varid)
      end if

      if (param%lenergy) then
         status = nf90_inq_varid(self%ncid, KE_ORB_VARNAME, self%KE_orb_varid)
         status = nf90_inq_varid(self%ncid, KE_SPIN_VARNAME, self%KE_spin_varid)
         status = nf90_inq_varid(self%ncid, PE_VARNAME, self%PE_varid)
         status = nf90_inq_varid(self%ncid, L_ORBX_VARNAME, self%L_orbx_varid)
         status = nf90_inq_varid(self%ncid, L_ORBY_VARNAME, self%L_orby_varid)
         status = nf90_inq_varid(self%ncid, L_ORBZ_VARNAME, self%L_orbz_varid)
         status = nf90_inq_varid(self%ncid, L_SPINX_VARNAME, self%L_spinx_varid)
         status = nf90_inq_varid(self%ncid, L_SPINY_VARNAME, self%L_spiny_varid)
         status = nf90_inq_varid(self%ncid, L_SPINZ_VARNAME, self%L_spinz_varid)
         status = nf90_inq_varid(self%ncid, L_ESCAPEX_VARNAME, self%L_escapex_varid)
         status = nf90_inq_varid(self%ncid, L_ESCAPEY_VARNAME, self%L_escapey_varid)
         status = nf90_inq_varid(self%ncid, L_ESCAPEZ_VARNAME, self%L_escapez_varid)
         status = nf90_inq_varid(self%ncid, ECOLLISIONS_VARNAME, self%Ecollisions_varid)
         status = nf90_inq_varid(self%ncid, EUNTRACKED_VARNAME, self%Euntracked_varid)
         status = nf90_inq_varid(self%ncid, GMESCAPE_VARNAME, self%GMescape_varid)
      end if

      return
   end subroutine netcdf_open


   module function netcdf_read_frame_system(self, iu, param) result(ierr)
      !! author: The Purdue Swiftest Team - David A. Minton, Carlisle A. Wishard, Jennifer L.L. Pouplin, and Jacob R. Elliott
      !!
      !! Read a frame (header plus records for each massive body and active test particle) from an output binary file
      implicit none
      ! Arguments
      class(swiftest_nbody_system), intent(inout) :: self  !! Swiftest system object
      class(netcdf_parameters),     intent(inout) :: iu    !! Parameters used to identify a particular NetCDF dataset
      class(swiftest_parameters),   intent(inout) :: param !! Current run configuration parameters 
      ! Return
      integer(I4B)                                :: ierr  !! Error code: returns 0 if the read is successful
      ! Internals
      integer(I4B)                              :: tslot, idmax, npl_check, ntp_check, nplm_check, t_max, str_max, status
      real(DP), dimension(:), allocatable       :: rtemp
      integer(I4B), dimension(:), allocatable   :: itemp
      logical, dimension(:), allocatable        :: validmask, tpmask, plmask

      call iu%open(param, readonly=.true.)
      call self%read_hdr(iu, param)

      associate(cb => self%cb, pl => self%pl, tp => self%tp, npl => self%pl%nbody, ntp => self%tp%nbody)

         call pl%setup(npl, param)
         call tp%setup(ntp, param)

         tslot = int(param%ioutput, kind=I4B) + 1

         call check( nf90_inquire_dimension(iu%ncid, iu%id_dimid, len=idmax), "netcdf_read_frame_system nf90_inquire_dimension id_dimid"  )
         allocate(rtemp(idmax))
         allocate(itemp(idmax))
         allocate(validmask(idmax))
         allocate(tpmask(idmax))
         allocate(plmask(idmax))
         call check( nf90_inquire_dimension(iu%ncid, iu%time_dimid, len=t_max), "netcdf_read_frame_system nf90_inquire_dimension time_dimid"  )
         call check( nf90_inquire_dimension(iu%ncid, iu%str_dimid, len=str_max), "netcdf_read_frame_system nf90_inquire_dimension str_dimid"  )

         ! First filter out only the id slots that contain valid bodies
         if (param%in_form == XV) then
            call check( nf90_get_var(iu%ncid, iu%xhx_varid, rtemp(:), start=[1, tslot]), "netcdf_read_frame_system filter pass nf90_getvar xhx_varid"  )
         else
            call check( nf90_get_var(iu%ncid, iu%a_varid, rtemp(:), start=[1, tslot]), "netcdf_read_frame_system filter pass nf90_getvar a_varid"  )
         end if

         validmask(:) = rtemp(:) == rtemp(:)

         ! Next, filter only bodies that don't have mass (test particles)
         call check( nf90_get_var(iu%ncid, iu%Gmass_varid, rtemp(:), start=[1, tslot]), "netcdf_read_frame_system nf90_getvar Gmass_varid"  )
         plmask(:) = rtemp(:) == rtemp(:)  .and. validmask(:)
         tpmask(:) = .not. plmask(:) .and. validmask(:)
         plmask(1) = .false. ! This is the central body

         ! Check to make sure the number of bodies is correct
         npl_check = count(plmask(:))
         ntp_check = count(tpmask(:))

         if (npl_check /= npl) then
            write(*,*) "Error reading in NetCDF file: The recorded value of npl does not match the number of active massive bodies"
            call util_exit(failure)
         end if

         if (ntp_check /= ntp) then
            write(*,*) "Error reading in NetCDF file: The recorded value of ntp does not match the number of active test particles"
            call util_exit(failure)
         end if

         select type (pl)
         class is (symba_pl)
            select type (param)
            class is (symba_parameters)
               nplm_check = count(pack(rtemp,plmask) > param%GMTINY )
               if (nplm_check /= pl%nplm) then
                  write(*,*) "Error reading in NetCDF file: The recorded value of nplm does not match the number of active fully interacting massive bodies"
                  call util_exit(failure)
               end if
            end select
         end select

         ! Now read in each variable and split the outputs by body type
         if ((param%in_form == XV) .or. (param%in_form == XVEL)) then
            call check( nf90_get_var(iu%ncid, iu%xhx_varid, rtemp, start=[1, tslot]), "netcdf_read_frame_system nf90_getvar xhx_varid"  )
            if (npl > 0) pl%xh(1,:) = pack(rtemp, plmask)
            if (ntp > 0) tp%xh(1,:) = pack(rtemp, tpmask)

            call check( nf90_get_var(iu%ncid, iu%xhy_varid, rtemp, start=[1, tslot]), "netcdf_read_frame_system nf90_getvar xhy_varid"  )
            if (npl > 0) pl%xh(2,:) = pack(rtemp, plmask)
            if (ntp > 0) tp%xh(2,:) = pack(rtemp, tpmask)

            call check( nf90_get_var(iu%ncid, iu%xhz_varid, rtemp, start=[1, tslot]), "netcdf_read_frame_system nf90_getvar xhz_varid"  )
            if (npl > 0) pl%xh(3,:) = pack(rtemp, plmask)
            if (ntp > 0) tp%xh(3,:) = pack(rtemp, tpmask)

            if (param%lgr .and. iu%lpseudo_vel_exists) then
               call check( nf90_get_var(iu%ncid, iu%gr_pseudo_vhx_varid, rtemp, start=[1, tslot]), "netcdf_read_frame_system nf90_getvar gr_pseudo_vhx_varid"  )
               if (npl > 0) pl%vh(1,:) = pack(rtemp, plmask)
               if (ntp > 0) tp%vh(1,:) = pack(rtemp, tpmask)

               call check( nf90_get_var(iu%ncid, iu%gr_pseudo_vhy_varid, rtemp, start=[1, tslot]), "netcdf_read_frame_system nf90_getvar gr_pseudo_vhy_varid"  )
               if (npl > 0) pl%vh(2,:) = pack(rtemp, plmask)
               if (ntp > 0) tp%vh(2,:) = pack(rtemp, tpmask)

               call check( nf90_get_var(iu%ncid, iu%gr_pseudo_vhz_varid, rtemp, start=[1, tslot]), "netcdf_read_frame_system nf90_getvar gr_pseudo_vhz_varid"  )
               if (npl > 0) pl%vh(3,:) = pack(rtemp, plmask)
               if (ntp > 0) tp%vh(3,:) = pack(rtemp, tpmask)
            else
               call check( nf90_get_var(iu%ncid, iu%vhx_varid, rtemp, start=[1, tslot]), "netcdf_read_frame_system nf90_getvar vhx_varid"  )
               if (npl > 0) pl%vh(1,:) = pack(rtemp, plmask)
               if (ntp > 0) tp%vh(1,:) = pack(rtemp, tpmask)

               call check( nf90_get_var(iu%ncid, iu%vhy_varid, rtemp, start=[1, tslot]), "netcdf_read_frame_system nf90_getvar vhy_varid"  )
               if (npl > 0) pl%vh(2,:) = pack(rtemp, plmask)
               if (ntp > 0) tp%vh(2,:) = pack(rtemp, tpmask)

               call check( nf90_get_var(iu%ncid, iu%vhz_varid, rtemp, start=[1, tslot]), "netcdf_read_frame_system nf90_getvar vhz_varid"  )
               if (npl > 0) pl%vh(3,:) = pack(rtemp, plmask)
               if (ntp > 0) tp%vh(3,:) = pack(rtemp, tpmask)
            end if
         end if

         if ((param%in_form == EL)  .or. (param%in_form == XVEL)) then
            call check( nf90_get_var(iu%ncid, iu%a_varid, rtemp, start=[1, tslot]), "netcdf_read_frame_system nf90_getvar a_varid"  )
            if (.not.allocated(pl%a)) allocate(pl%a(npl))
            if (.not.allocated(tp%a)) allocate(tp%a(ntp))
            if (npl > 0) pl%a(:) = pack(rtemp, plmask)
            if (ntp > 0) tp%a(:) = pack(rtemp, tpmask)

            call check( nf90_get_var(iu%ncid, iu%e_varid, rtemp, start=[1, tslot]), "netcdf_read_frame_system nf90_getvar e_varid"  )
            if (.not.allocated(pl%e)) allocate(pl%e(npl))
            if (.not.allocated(tp%e)) allocate(tp%e(ntp))
            if (npl > 0) pl%e(:) = pack(rtemp, plmask)
            if (ntp > 0) tp%e(:) = pack(rtemp, tpmask)

            call check( nf90_get_var(iu%ncid, iu%inc_varid, rtemp, start=[1, tslot]), "netcdf_read_frame_system nf90_getvar inc_varid"  )
            if (.not.allocated(pl%inc)) allocate(pl%inc(npl))
            if (.not.allocated(tp%inc)) allocate(tp%inc(ntp))
            if (npl > 0) pl%inc(:) = pack(rtemp, plmask)
            if (ntp > 0) tp%inc(:) = pack(rtemp, tpmask)

            call check( nf90_get_var(iu%ncid, iu%capom_varid, rtemp, start=[1, tslot]), "netcdf_read_frame_system nf90_getvar capom_varid"  )
            if (.not.allocated(pl%capom)) allocate(pl%capom(npl))
            if (.not.allocated(tp%capom)) allocate(tp%capom(ntp))
            if (npl > 0) pl%capom(:) = pack(rtemp, plmask)
            if (ntp > 0) tp%capom(:) = pack(rtemp, tpmask)

            call check( nf90_get_var(iu%ncid, iu%omega_varid, rtemp, start=[1, tslot]), "netcdf_read_frame_system nf90_getvar omega_varid"  )
            if (.not.allocated(pl%omega)) allocate(pl%omega(npl))
            if (.not.allocated(tp%omega)) allocate(tp%omega(ntp))
            if (npl > 0) pl%omega(:) = pack(rtemp, plmask)
            if (ntp > 0) tp%omega(:) = pack(rtemp, tpmask)

            call check( nf90_get_var(iu%ncid, iu%capm_varid, rtemp, start=[1, tslot]), "netcdf_read_frame_system nf90_getvar capm_varid"  )
            if (.not.allocated(pl%capm)) allocate(pl%capm(npl))
            if (.not.allocated(tp%capm)) allocate(tp%capm(ntp))
            if (npl > 0) pl%capm(:) = pack(rtemp, plmask)
            if (ntp > 0) tp%capm(:) = pack(rtemp, tpmask)

         end if
      
         call check( nf90_get_var(iu%ncid, iu%Gmass_varid, rtemp, start=[1, tslot]), "netcdf_read_frame_system nf90_getvar Gmass_varid"  )
         cb%Gmass = rtemp(1)
         cb%mass = cb%Gmass / param%GU

         ! Set initial central body mass for Helio bookkeeping
         select type(cb)
            class is (symba_cb)
               cb%GM0 = cb%Gmass
         end select
            

         if (npl > 0) then
            pl%Gmass(:) = pack(rtemp, plmask)
            pl%mass(:) = pl%Gmass(:) / param%GU

            if (param%lrhill_present) then 
               call check( nf90_get_var(iu%ncid, iu%rhill_varid, rtemp, start=[1, tslot]), "netcdf_read_frame_system nf90_getvar rhill_varid"  )
               pl%rhill(:) = pack(rtemp, plmask)
            end if
         end if

         if (param%lclose) then
            call check( nf90_get_var(iu%ncid, iu%radius_varid, rtemp, start=[1, tslot]), "netcdf_read_frame_system nf90_getvar radius_varid"  )
            cb%radius = rtemp(1)

            ! Set initial central body radius for SyMBA bookkeeping
            select type(cb)
               class is (symba_cb)
                  cb%R0 = cb%radius
            end select
            if (npl > 0) pl%radius(:) = pack(rtemp, plmask)
         else
            cb%radius = param%rmin
            if (npl > 0) pl%radius(:) = 0.0_DP
         end if

         if (param%lrotation) then
            call check( nf90_get_var(iu%ncid, iu%Ip1_varid,  rtemp, start=[1, tslot]), "netcdf_read_frame_system nf90_getvar Ip1_varid"  )
            cb%Ip(1) = rtemp(1)
            if (npl > 0) pl%Ip(1,:) = pack(rtemp, plmask)

            call check( nf90_get_var(iu%ncid, iu%Ip2_varid,  rtemp, start=[1, tslot]), "netcdf_read_frame_system nf90_getvar Ip2_varid"  )
            cb%Ip(2) = rtemp(1)
            if (npl > 0) pl%Ip(2,:) = pack(rtemp, plmask)

            call check( nf90_get_var(iu%ncid, iu%Ip3_varid,  rtemp, start=[1, tslot]), "netcdf_read_frame_system nf90_getvar Ip3_varid"  )
            cb%Ip(3) = rtemp(1)
            if (npl > 0) pl%Ip(3,:) = pack(rtemp, plmask)

            call check( nf90_get_var(iu%ncid, iu%rotx_varid, rtemp, start=[1, tslot]), "netcdf_read_frame_system nf90_getvar rotx_varid"  )
            cb%rot(1) = rtemp(1)
            if (npl > 0) pl%rot(1,:) = pack(rtemp, plmask)

            call check( nf90_get_var(iu%ncid, iu%roty_varid, rtemp, start=[1, tslot]), "netcdf_read_frame_system nf90_getvar roty_varid"  )
            cb%rot(2) = rtemp(1)
            if (npl > 0) pl%rot(2,:) = pack(rtemp, plmask)

            call check( nf90_get_var(iu%ncid, iu%rotz_varid, rtemp, start=[1, tslot]), "netcdf_read_frame_system nf90_getvar rotz_varid"  )
            cb%rot(3) = rtemp(1)
            if (npl > 0) pl%rot(3,:) = pack(rtemp, plmask)

            ! Set initial central body angular momentum for Helio bookkeeping
            select type(cb)
               class is (symba_cb)
                  cb%L0(:) = cb%Ip(3) * cb%GM0 * cb%R0**2 * cb%rot(:)         
            end select
         end if

         ! if (param%ltides) then
         !    call check( nf90_get_var(iu%ncid, iu%k2_varid, rtemp, start=[1, tslot]), "netcdf_read_frame_system nf90_getvar k2_varid"  )
         !    cb%k2 = rtemp(1)
         !    if (npl > 0) pl%k2(:) = pack(rtemp, plmask)

         !    call check( nf90_get_var(iu%ncid, iu%Q_varid,  rtemp,  start=[1, tslot]), "netcdf_read_frame_system nf90_getvar Q_varid"  )
         !    cb%Q = rtemp(1)
         !    if (npl > 0) pl%Q(:) = pack(rtemp, plmask)
         ! end if

         status = nf90_inq_varid(iu%ncid, J2RP2_VARNAME, iu%j2rp2_varid)
         if (status == nf90_noerr) then
            call check( nf90_get_var(iu%ncid, iu%j2rp2_varid, cb%j2rp2, start=[tslot]), "netcdf_read_frame_system nf90_getvar j2rp2_varid"  )
         else 
            cb%j2rp2 = 0.0_DP
         end if

         status = nf90_inq_varid(iu%ncid, J4RP4_VARNAME, iu%j4rp4_varid)   
         if (status == nf90_noerr) then      
            call check( nf90_get_var(iu%ncid, iu%j4rp4_varid, cb%j4rp4, start=[tslot]), "netcdf_read_frame_system nf90_getvar j4rp4_varid"  )
         else 
            cb%j4rp4 = 0.0_DP
         end if

         call self%read_particle_info(iu, param, plmask, tpmask) 

         if (param%in_form == "EL") then
            call pl%el2xv(cb)
            call tp%el2xv(cb)
         end if
         ! if this is a GR-enabled run, check to see if we got the pseudovelocities in. Otherwise, we'll need to generate them.
         if (param%lgr .and. .not.(iu%lpseudo_vel_exists)) then
            call pl%set_mu(cb)
            call tp%set_mu(cb)
            call pl%v2pv(param)
            call tp%v2pv(param)
         end if
         
      end associate

      call iu%close()

      ierr = 0
      return

      667 continue
      write(*,*) "Error reading system frame in netcdf_read_frame_system"

   end function netcdf_read_frame_system


   module subroutine netcdf_read_hdr_system(self, iu, param) 
      !! author: David A. Minton
      !!
      !! Reads header information (variables that change with time, but not particle id). 
      !! This subroutine significantly improves the output over the original binary file, allowing us to track energy, momentum, and other quantities that 
      !! previously were handled as separate output files.
      implicit none
      ! Arguments
      class(swiftest_nbody_system), intent(inout) :: self  !! Swiftest nbody system object
      class(netcdf_parameters),     intent(inout) :: iu    !! Parameters used to for writing a NetCDF dataset to file
      class(swiftest_parameters),   intent(inout) :: param !! Current run configuration parameters
      ! Internals
      integer(I4B) :: tslot, status, idmax
      real(DP), dimension(:), allocatable       :: gmtemp
      logical, dimension(:), allocatable        :: plmask, tpmask, plmmask


      tslot = int(param%ioutput, kind=I4B) + 1
      call check( nf90_inquire_dimension(iu%ncid, iu%id_dimid, len=idmax), "netcdf_read_frame_system nf90_inquire_dimension id_dimid"  )
      call check( nf90_get_var(iu%ncid, iu%time_varid, param%t,       start=[tslot]), "netcdf_read_hdr_system nf90_getvar time_varid"  )

      allocate(gmtemp(idmax))
      allocate(tpmask(idmax))
      allocate(plmask(idmax))
      allocate(plmmask(idmax))

      call check( nf90_get_var(iu%ncid, iu%Gmass_varid, gmtemp, start=[1,1]), "netcdf_read_frame_system nf90_getvar Gmass_varid"  )

      plmask(:) = gmtemp(:) == gmtemp(:)
      tpmask(:) = .not. plmask(:)
      plmask(1) = .false. ! This is the central body
      select type (param)
      class is (symba_parameters)
         plmmask(:) = plmask(:)
         where(plmask(:))
            plmmask(:) = gmtemp(:) > param%GMTINY
         endwhere
      end select

      status = nf90_inq_varid(iu%ncid, NPL_VARNAME, iu%npl_varid)
      if (status == nf90_noerr) then
         call check( nf90_get_var(iu%ncid, iu%npl_varid,  self%pl%nbody, start=[tslot]), "netcdf_read_hdr_system nf90_getvar npl_varid"  )
      else
         self%pl%nbody = count(plmask(:))
      end if

      status = nf90_inq_varid(iu%ncid, NTP_VARNAME, iu%ntp_varid)
      if (status == nf90_noerr) then
         call check( nf90_get_var(iu%ncid, iu%ntp_varid,  self%tp%nbody, start=[tslot]), "netcdf_read_hdr_system nf90_getvar ntp_varid"  )
      else
         self%tp%nbody = count(tpmask(:))
      end if

      if (param%integrator == SYMBA) then
         status = nf90_inq_varid(iu%ncid, NPLM_VARNAME, iu%nplm_varid)
         select type(pl => self%pl)
         class is (symba_pl)
            if (status == nf90_noerr) then
               call check( nf90_get_var(iu%ncid, iu%nplm_varid,  pl%nplm, start=[tslot]), "netcdf_read_hdr_system nf90_getvar nplm_varid"  )
            else
               pl%nplm = count(plmmask(:))
            end if
         end select
      end if

      if (param%lenergy) then
         status = nf90_inq_varid(iu%ncid, KE_ORB_VARNAME, iu%KE_orb_varid)
         if (status == nf90_noerr) call check( nf90_get_var(iu%ncid, iu%KE_orb_varid,      self%ke_orbit,    start=[tslot]), "netcdf_read_hdr_system nf90_getvar KE_orb_varid"  )
         status = nf90_inq_varid(iu%ncid, KE_SPIN_VARNAME, iu%KE_spin_varid)
         if (status == nf90_noerr) call check( nf90_get_var(iu%ncid, iu%KE_spin_varid,     self%ke_spin,     start=[tslot]), "netcdf_read_hdr_system nf90_getvar KE_spin_varid"  )
         status = nf90_inq_varid(iu%ncid, PE_VARNAME, iu%PE_varid)
         if (status == nf90_noerr) call check( nf90_get_var(iu%ncid, iu%PE_varid,          self%pe,          start=[tslot]), "netcdf_read_hdr_system nf90_getvar PE_varid"  )
         status = nf90_inq_varid(iu%ncid, L_ORBX_VARNAME, iu%L_orbx_varid)
         if (status == nf90_noerr) call check( nf90_get_var(iu%ncid, iu%L_orbx_varid,      self%Lorbit(1),   start=[tslot]), "netcdf_read_hdr_system nf90_getvar L_orbx_varid"  )
         status = nf90_inq_varid(iu%ncid, L_ORBY_VARNAME, iu%L_orby_varid)
         if (status == nf90_noerr) call check( nf90_get_var(iu%ncid, iu%L_orby_varid,      self%Lorbit(2),   start=[tslot]), "netcdf_read_hdr_system nf90_getvar L_orby_varid"  )
         status = nf90_inq_varid(iu%ncid, L_ORBZ_VARNAME, iu%L_orbz_varid)
         if (status == nf90_noerr) call check( nf90_get_var(iu%ncid, iu%L_orbz_varid,      self%Lorbit(3),   start=[tslot]), "netcdf_read_hdr_system nf90_getvar L_orbz_varid"  )
         status = nf90_inq_varid(iu%ncid, L_SPINX_VARNAME, iu%L_spinx_varid)
         if (status == nf90_noerr) call check( nf90_get_var(iu%ncid, iu%L_spinx_varid,     self%Lspin(1),    start=[tslot]), "netcdf_read_hdr_system nf90_getvar L_spinx_varid"  )
         status = nf90_inq_varid(iu%ncid, L_SPINY_VARNAME, iu%L_spiny_varid)
         if (status == nf90_noerr) call check( nf90_get_var(iu%ncid, iu%L_spiny_varid,     self%Lspin(2),    start=[tslot]), "netcdf_read_hdr_system nf90_getvar L_spiny_varid"  )
         status = nf90_inq_varid(iu%ncid, L_SPINZ_VARNAME, iu%L_spinz_varid)
         if (status == nf90_noerr) call check( nf90_get_var(iu%ncid, iu%L_spinz_varid,     self%Lspin(3),    start=[tslot]), "netcdf_read_hdr_system nf90_getvar L_spinz_varid"  )
         status = nf90_inq_varid(iu%ncid, L_ESCAPEX_VARNAME, iu%L_escapex_varid)
         if (status == nf90_noerr) call check( nf90_get_var(iu%ncid, iu%L_escapex_varid,   self%Lescape(1),  start=[tslot]), "netcdf_read_hdr_system nf90_getvar L_escapex_varid"  )
         status = nf90_inq_varid(iu%ncid, L_ESCAPEY_VARNAME, iu%L_escapey_varid)
         if (status == nf90_noerr) call check( nf90_get_var(iu%ncid, iu%L_escapey_varid,   self%Lescape(2),  start=[tslot]), "netcdf_read_hdr_system nf90_getvar L_escapey_varid"  )
         status = nf90_inq_varid(iu%ncid, L_ESCAPEZ_VARNAME, iu%L_escapez_varid)
         if (status == nf90_noerr) call check( nf90_get_var(iu%ncid, iu%L_escapez_varid,   self%Lescape(3),  start=[tslot]), "netcdf_read_hdr_system nf90_getvar L_escapez_varid"  )
         status = nf90_inq_varid(iu%ncid, ECOLLISIONS_VARNAME, iu%Ecollisions_varid)
         if (status == nf90_noerr) call check( nf90_get_var(iu%ncid, iu%Ecollisions_varid, self%Ecollisions, start=[tslot]), "netcdf_read_hdr_system nf90_getvar Ecollisions_varid"  )
         status = nf90_inq_varid(iu%ncid, EUNTRACKED_VARNAME, iu%Euntracked_varid)
         if (status == nf90_noerr) call check( nf90_get_var(iu%ncid, iu%Euntracked_varid,  self%Euntracked,  start=[tslot]), "netcdf_read_hdr_system nf90_getvar Euntracked_varid"  )
         status = nf90_inq_varid(iu%ncid, GMESCAPE_VARNAME, iu%GMescape_varid)
         if (status == nf90_noerr)  call check( nf90_get_var(iu%ncid, iu%GMescape_varid,    self%GMescape,    start=[tslot]), "netcdf_read_hdr_system nf90_getvar GMescape_varid"  )
      end if

      return
   end subroutine netcdf_read_hdr_system


   module subroutine netcdf_read_particle_info_system(self, iu, param, plmask, tpmask)
      !! author: Carlisle A. Wishard, Dana Singh, and David A. Minton
      !!
      !! Reads particle information metadata from file
      implicit none
      ! Arguments
      class(swiftest_nbody_system), intent(inout) :: self   !! Swiftest nbody system object
      class(netcdf_parameters),     intent(inout) :: iu     !! Parameters used to identify a particular NetCDF dataset
      class(swiftest_parameters),   intent(inout) :: param  !! Current run configuration parameters
      logical, dimension(:),        intent(in)    :: plmask !! Logical array indicating which index values belong to massive bodies
      logical, dimension(:),        intent(in)    :: tpmask !! Logical array indicating which index values belong to test particles
      ! Internals
      integer(I4B)                                :: i, idmax, status
      real(DP), dimension(:), allocatable         :: rtemp
      real(DP), dimension(:,:), allocatable       :: rtemp_arr
      integer(I4B), dimension(:), allocatable     :: itemp
      character(len=NAMELEN), dimension(:), allocatable :: ctemp
      integer(I4B), dimension(:), allocatable     :: plind, tpind

      ! This string of spaces of length NAMELEN is used to clear out any old data left behind inside the string variables
      idmax = size(plmask)
      allocate(rtemp(idmax))
      allocate(rtemp_arr(NDIM,idmax))
      allocate(itemp(idmax))
      allocate(ctemp(idmax))

      associate(cb => self%cb, pl => self%pl, tp => self%tp, npl => self%pl%nbody, ntp => self%tp%nbody)

         if (npl > 0) then
            pl%status(:) = ACTIVE
            pl%lmask(:) = .true.
            do i = 1, npl
               call pl%info(i)%set_value(status="ACTIVE")
            end do
            allocate(plind(npl))
            plind(:) = pack([(i, i = 1, idmax)], plmask(:))
         end if
         if (ntp > 0) then
            tp%status(:) = ACTIVE
            tp%lmask(:) = .true.
            do i = 1, ntp
               call tp%info(i)%set_value(status="ACTIVE")
            end do
            allocate(tpind(ntp))
            tpind(:) = pack([(i, i = 1, idmax)], tpmask(:))
         end if

         call check( nf90_get_var(iu%ncid, iu%id_varid, itemp), "netcdf_read_particle_info_system nf90_getvar id_varid"  )
         cb%id = itemp(1)
         pl%id(:) = pack(itemp, plmask)
         tp%id(:) = pack(itemp, tpmask)

         call check( nf90_get_var(iu%ncid, iu%name_varid, ctemp, count=[NAMELEN, idmax]), "netcdf_read_particle_info_system nf90_getvar name_varid"  )
         call cb%info%set_value(name=ctemp(1))
         do i = 1, npl
            call pl%info(i)%set_value(name=ctemp(plind(i)))
         end do
         do i = 1, ntp
            call tp%info(i)%set_value(name=ctemp(tpind(i)))
         end do

         call check( nf90_get_var(iu%ncid, iu%ptype_varid, ctemp, count=[NAMELEN, idmax]), "netcdf_read_particle_info_system nf90_getvar ptype_varid"  )
         call cb%info%set_value(particle_type=ctemp(1))
         do i = 1, npl
            call pl%info(i)%set_value(particle_type=ctemp(plind(i)))
         end do
         do i = 1, ntp
            call tp%info(i)%set_value(particle_type=ctemp(tpind(i)))
         end do

         status = nf90_inq_varid(iu%ncid, STATUS_VARNAME, iu%status_varid) 
         if (status == nf90_noerr) then
            call check( nf90_get_var(iu%ncid, iu%status_varid, ctemp, count=[NAMELEN, idmax]), "netcdf_read_particle_info_system nf90_getvar status_varid")
            call cb%info%set_value(status=ctemp(1))
         else
            call cb%info%set_value(status="ACTIVE")
         end if
         do i = 1, npl
            call pl%info(i)%set_value(status=ctemp(plind(i)))
         end do
         do i = 1, ntp
            call tp%info(i)%set_value(status=ctemp(tpind(i)))
         end do

         if (param%lclose) then

            status = nf90_inq_varid(iu%ncid, ORIGIN_TYPE_VARNAME, iu%origin_type_varid)
            if (status == nf90_noerr) then
               call check( nf90_get_var(iu%ncid, iu%origin_type_varid, ctemp, count=[NAMELEN, idmax]), "netcdf_read_particle_info_system nf90_getvar origin_type_varid"  )
            else
               ctemp = "Initial Conditions"
            end if

            call cb%info%set_value(origin_type=ctemp(1))
            do i = 1, npl
               call pl%info(i)%set_value(origin_type=ctemp(plind(i)))
            end do
            do i = 1, ntp
               call tp%info(i)%set_value(origin_type=ctemp(tpind(i)))
            end do

            status = nf90_inq_varid(iu%ncid, ORIGIN_TIME_VARNAME, iu%origin_time_varid)
            if (status == nf90_noerr) then
               call check( nf90_get_var(iu%ncid, iu%origin_time_varid, rtemp), "netcdf_read_particle_info_system nf90_getvar origin_time_varid"  )
            else
               rtemp = param%t0
            end if

            call cb%info%set_value(origin_time=rtemp(1))
            do i = 1, npl
               call pl%info(i)%set_value(origin_time=rtemp(plind(i)))
            end do
            do i = 1, ntp
               call tp%info(i)%set_value(origin_time=rtemp(tpind(i)))
            end do

            status = nf90_inq_varid(iu%ncid, ORIGIN_XHX_VARNAME, iu%origin_xhx_varid)
            if (status == nf90_noerr) then
               call check( nf90_get_var(iu%ncid, iu%origin_xhx_varid, rtemp_arr(1,:)), "netcdf_read_particle_info_system nf90_getvar origin_xhx_varid"  )
            else if ((param%out_form == XV) .or. (param%out_form == XVEL)) then
               call check( nf90_get_var(iu%ncid, iu%xhx_varid, rtemp_arr(1,:)), "netcdf_read_particle_info_system nf90_getvar xhx_varid"  )
            else 
               rtemp_arr(1,:) = 0._DP
            end if 

            status = nf90_inq_varid(iu%ncid, ORIGIN_XHY_VARNAME, iu%origin_xhy_varid)
            if (status == nf90_noerr) then
               call check( nf90_get_var(iu%ncid, iu%origin_xhy_varid, rtemp_arr(2,:)), "netcdf_read_particle_info_system nf90_getvar origin_xhy_varid"  )
            else if ((param%out_form == XV) .or. (param%out_form == XVEL)) then
               call check( nf90_get_var(iu%ncid, iu%xhy_varid, rtemp_arr(2,:)), "netcdf_read_particle_info_system nf90_getvar xhx_varid"  )
            else 
               rtemp_arr(2,:) = 0._DP
            end if 

            status = nf90_inq_varid(iu%ncid, ORIGIN_XHZ_VARNAME, iu%origin_xhz_varid)
            if (status == nf90_noerr) then
               call check( nf90_get_var(iu%ncid, iu%origin_xhz_varid, rtemp_arr(3,:)), "netcdf_read_particle_info_system nf90_getvar origin_xhz_varid"  )
            else if ((param%out_form == XV) .or. (param%out_form == XVEL)) then
               call check( nf90_get_var(iu%ncid, iu%xhz_varid, rtemp_arr(3,:)), "netcdf_read_particle_info_system nf90_getvar xhz_varid"  )
            else
               rtemp_arr(3,:) = 0._DP
            end if 

            do i = 1, npl
               call pl%info(i)%set_value(origin_xh=rtemp_arr(:,plind(i)))
            end do
            do i = 1, ntp
               call tp%info(i)%set_value(origin_xh=rtemp_arr(:,tpind(i)))
            end do

            status = nf90_inq_varid(iu%ncid, ORIGIN_VHX_VARNAME, iu%origin_vhx_varid)
            if (status == nf90_noerr) then
               call check( nf90_get_var(iu%ncid, iu%origin_vhx_varid, rtemp_arr(1,:)), "netcdf_read_particle_info_system nf90_getvar origin_vhx_varid"  )
            else if ((param%out_form == XV) .or. (param%out_form == XVEL)) then
               call check( nf90_get_var(iu%ncid, iu%vhx_varid, rtemp_arr(1,:)), "netcdf_read_particle_info_system nf90_getvar vhx_varid"  )
            else
               rtemp_arr(1,:) = 0._DP
            end if 
            
            status = nf90_inq_varid(iu%ncid, ORIGIN_VHY_VARNAME, iu%origin_vhy_varid)
            if (status == nf90_noerr) then
               call check( nf90_get_var(iu%ncid, iu%origin_vhy_varid, rtemp_arr(2,:)), "netcdf_read_particle_info_system nf90_getvar origin_vhy_varid"  )
            else if ((param%out_form == XV) .or. (param%out_form == XVEL)) then
               call check( nf90_get_var(iu%ncid, iu%vhy_varid, rtemp_arr(2,:)), "netcdf_read_particle_info_system nf90_getvar vhy_varid"  )
            else
               rtemp_arr(2,:) = 0._DP
            end if 

            status = nf90_inq_varid(iu%ncid, ORIGIN_VHZ_VARNAME, iu%origin_vhz_varid)
            if (status == nf90_noerr) then
               call check( nf90_get_var(iu%ncid, iu%origin_vhz_varid, rtemp_arr(3,:)), "netcdf_read_particle_info_system nf90_getvar origin_vhz_varid"  )
            else if ((param%out_form == XV) .or. (param%out_form == XVEL)) then
               call check( nf90_get_var(iu%ncid, iu%vhz_varid, rtemp_arr(3,:)), "netcdf_read_particle_info_system nf90_getvar vhz_varid" )
            else
               rtemp_arr(3,:) = 0._DP
            end if 

            do i = 1, npl
               call pl%info(i)%set_value(origin_vh=rtemp_arr(:,plind(i)))
            end do
            do i = 1, ntp
               call tp%info(i)%set_value(origin_vh=rtemp_arr(:,tpind(i)))
            end do

            status = nf90_inq_varid(iu%ncid, COLLISION_ID_VARNAME, iu%collision_id_varid)
            if (status == nf90_noerr) then
               call check( nf90_get_var(iu%ncid, iu%collision_id_varid, itemp), "netcdf_read_particle_info_system nf90_getvar collision_id_varid"  )
            else
               itemp = 0.0_DP
            end if 

            do i = 1, npl
               call pl%info(i)%set_value(collision_id=itemp(plind(i)))
            end do
            do i = 1, ntp
               call tp%info(i)%set_value(collision_id=itemp(tpind(i)))
            end do

            status = nf90_inq_varid(iu%ncid, DISCARD_TIME_VARNAME, iu%discard_time_varid)
            if (status == nf90_noerr) then
               call check( nf90_get_var(iu%ncid, iu%discard_time_varid, rtemp), "netcdf_read_particle_info_system nf90_getvar discard_time_varid"  )
            else
               rtemp = 0.0_DP
            end if 

            call cb%info%set_value(discard_time=rtemp(1))
            do i = 1, npl
               call pl%info(i)%set_value(discard_time=rtemp(plind(i)))
            end do
            do i = 1, ntp
               call tp%info(i)%set_value(discard_time=rtemp(tpind(i)))
            end do

            status = nf90_inq_varid(iu%ncid, DISCARD_XHX_VARNAME, iu%discard_xhx_varid)
            if (status == nf90_noerr) then
               call check( nf90_get_var(iu%ncid, iu%discard_xhx_varid, rtemp_arr(1,:)), "netcdf_read_particle_info_system nf90_getvar discard_xhx_varid"  )
            else
               rtemp_arr(1,:) = 0.0_DP
            end if 

            status = nf90_inq_varid(iu%ncid, DISCARD_XHY_VARNAME, iu%discard_xhy_varid)
            if (status == nf90_noerr) then
               call check( nf90_get_var(iu%ncid, iu%discard_xhy_varid, rtemp_arr(2,:)), "netcdf_read_particle_info_system nf90_getvar discard_xhy_varid"  )
            else
               rtemp_arr(2,:) = 0.0_DP
            end if 

            status = nf90_inq_varid(iu%ncid, DISCARD_XHZ_VARNAME, iu%discard_xhz_varid)
            if (status == nf90_noerr) then
               call check( nf90_get_var(iu%ncid, iu%discard_xhz_varid, rtemp_arr(3,:)), "netcdf_read_particle_info_system nf90_getvar discard_xhz_varid"  )
            else
               rtemp_arr(3,:) = 0.0_DP
            end if 

            do i = 1, npl
               call pl%info(i)%set_value(discard_xh=rtemp_arr(:,plind(i)))
            end do
            do i = 1, ntp
               call tp%info(i)%set_value(discard_xh=rtemp_arr(:,tpind(i)))
            end do

            status = nf90_inq_varid(iu%ncid, DISCARD_VHX_VARNAME, iu%discard_vhx_varid)
            if (status == nf90_noerr) then
               call check( nf90_get_var(iu%ncid, iu%discard_vhx_varid, rtemp_arr(1,:)), "netcdf_read_particle_info_system nf90_getvar discard_vhx_varid"  )
            else
               rtemp_arr(1,:) = 0.0_DP
            end if 

            status = nf90_inq_varid(iu%ncid, DISCARD_VHY_VARNAME, iu%discard_vhy_varid)
            if (status == nf90_noerr) then
               call check( nf90_get_var(iu%ncid, iu%discard_vhy_varid, rtemp_arr(2,:)), "netcdf_read_particle_info_system nf90_getvar discard_vhy_varid"  )
            else
               rtemp_arr(2,:) = 0.0_DP
            end if 

            status = nf90_inq_varid(iu%ncid, DISCARD_VHZ_VARNAME, iu%discard_vhz_varid)
            if (status == nf90_noerr) then
               call check( nf90_get_var(iu%ncid, iu%discard_vhz_varid, rtemp_arr(3,:)), "netcdf_read_particle_info_system nf90_getvar discard_vhz_varid"  )
            else
               rtemp_arr(3,:) = 0.0_DP
            end if 

            do i = 1, npl
               call pl%info(i)%set_value(discard_vh=rtemp_arr(:,plind(i)))
            end do
            do i = 1, ntp
               call tp%info(i)%set_value(discard_vh=rtemp_arr(:,tpind(i)))
            end do
         end if

      end associate

      return
   end subroutine netcdf_read_particle_info_system


   module subroutine netcdf_sync(self)
      !! author: David A. Minton
      !!
      !! Syncrhonize the disk and memory buffer of the NetCDF file (e.g. commit the frame files stored in memory to disk) 
      !!    
      implicit none
      ! Arguments
      class(netcdf_parameters),   intent(inout) :: self !! Parameters used to identify a particular NetCDF dataset

      call check( nf90_sync(self%ncid), "netcdf_sync nf90_sync"  )

      return
   end subroutine netcdf_sync


   module subroutine netcdf_write_frame_base(self, iu, param)
      !! author: Carlisle A. Wishard, Dana Singh, and David A. Minton
      !!
      !! Write a frame of output of either test particle or massive body data to the binary output file
      !!    Note: If outputting to orbital elements, but sure that the conversion is done prior to calling this method
      implicit none
      ! Arguments
      class(swiftest_base),       intent(in)    :: self   !! Swiftest particle object
      class(netcdf_parameters),   intent(inout) :: iu     !! Parameters used to identify a particular NetCDF dataset
      class(swiftest_parameters), intent(inout) :: param  !! Current run configuration parameters
      ! Internals
      integer(I4B)                              :: i, j, tslot, idslot, old_mode
      integer(I4B), dimension(:), allocatable   :: ind
      real(DP), dimension(NDIM)                 :: vh !! Temporary variable to store heliocentric velocity values when converting from pseudovelocity in GR-enabled runs
      real(DP)                                  :: a, e, inc, omega, capom, capm

      call self%write_particle_info(iu, param)

      tslot = int(param%ioutput, kind=I4B) + 1

      call check( nf90_set_fill(iu%ncid, nf90_nofill, old_mode), "netcdf_write_frame_base nf90_set_fill"  )
      select type(self)
         class is (swiftest_body)
         associate(n => self%nbody)
            if (n == 0) return

            call util_sort(self%id(1:n), ind)

            do i = 1, n
               j = ind(i)
               idslot = self%id(j) + 1

               !! Convert from pseudovelocity to heliocentric without replacing the current value of pseudovelocity 
               if (param%lgr) call gr_pseudovel2vel(param, self%mu(j), self%xh(:, j), self%vh(:, j), vh(:))

               if ((param%out_form == XV) .or. (param%out_form == XVEL)) then
                  call check( nf90_put_var(iu%ncid, iu%xhx_varid, self%xh(1, j), start=[idslot, tslot]), "netcdf_write_frame_base nf90_put_var xhx_varid"  )
                  call check( nf90_put_var(iu%ncid, iu%xhy_varid, self%xh(2, j), start=[idslot, tslot]), "netcdf_write_frame_base nf90_put_var xhy_varid"  )
                  call check( nf90_put_var(iu%ncid, iu%xhz_varid, self%xh(3, j), start=[idslot, tslot]), "netcdf_write_frame_base nf90_put_var xhz_varid"  )
                  if (param%lgr) then !! Convert from pseudovelocity to heliocentric without replacing the current value of pseudovelocity
                     call check( nf90_put_var(iu%ncid, iu%vhx_varid, vh(1), start=[idslot, tslot]), "netcdf_write_frame_base nf90_put_var vhx_varid (gr case)"  )
                     call check( nf90_put_var(iu%ncid, iu%vhy_varid, vh(2), start=[idslot, tslot]), "netcdf_write_frame_base nf90_put_var vhy_varid (gr case)"  )
                     call check( nf90_put_var(iu%ncid, iu%vhz_varid, vh(3), start=[idslot, tslot]), "netcdf_write_frame_base nf90_put_var vhz_varid (gr case)"  )
                     call check( nf90_put_var(iu%ncid, iu%gr_pseudo_vhx_varid, self%vh(1, j), start=[idslot, tslot]), "netcdf_write_frame_base nf90_put_var gr_pseudo_vhx_varid"  )
                     call check( nf90_put_var(iu%ncid, iu%gr_pseudo_vhy_varid, self%vh(2, j), start=[idslot, tslot]), "netcdf_write_frame_base nf90_put_var gr_pseudo_vhy_varid"  )
                     call check( nf90_put_var(iu%ncid, iu%gr_pseudo_vhz_varid, self%vh(3, j), start=[idslot, tslot]), "netcdf_write_frame_base nf90_put_var gr_pseudo_vhz_varid"  )

                  else
                     call check( nf90_put_var(iu%ncid, iu%vhx_varid, self%vh(1, j), start=[idslot, tslot]), "netcdf_write_frame_base nf90_put_var vhx_varid"  )
                     call check( nf90_put_var(iu%ncid, iu%vhy_varid, self%vh(2, j), start=[idslot, tslot]), "netcdf_write_frame_base nf90_put_var vhy_varid"  )
                     call check( nf90_put_var(iu%ncid, iu%vhz_varid, self%vh(3, j), start=[idslot, tslot]), "netcdf_write_frame_base nf90_put_var vhz_varid"  )
                  end if
               end if

               if ((param%out_form == EL) .or. (param%out_form == XVEL)) then
                  if (param%lgr) then !! For GR-enabled runs, use the true value of velocity computed above
                     call orbel_xv2el(self%mu(j), self%xh(1,j), self%xh(2,j), self%xh(3,j), &
                                       vh(1), vh(2), vh(3), &
                                       a, e, inc, capom, omega, capm)
                  else !! For non-GR runs just convert from the velocity we have
                     call orbel_xv2el(self%mu(j), self%xh(1,j), self%xh(2,j), self%xh(3,j), &
                                       self%vh(1,j), self%vh(2,j), self%vh(3,j), &
                                       a, e, inc, capom, omega, capm)
                  end if
                  call check( nf90_put_var(iu%ncid, iu%a_varid, a, start=[idslot, tslot]), "netcdf_write_frame_base nf90_put_var a_varid"  )
                  call check( nf90_put_var(iu%ncid, iu%e_varid, e, start=[idslot, tslot]), "netcdf_write_frame_base nf90_put_var e_varid"  )
                  call check( nf90_put_var(iu%ncid, iu%inc_varid, inc * RAD2DEG, start=[idslot, tslot]), "netcdf_write_frame_base nf90_put_var inc_varid"  )
                  call check( nf90_put_var(iu%ncid, iu%capom_varid, capom * RAD2DEG, start=[idslot, tslot]), "netcdf_write_frame_base nf90_put_var capom_varid"  )
                  call check( nf90_put_var(iu%ncid, iu%omega_varid, omega * RAD2DEG, start=[idslot, tslot]), "netcdf_write_frame_base nf90_put_var omega_varid"  )
                  call check( nf90_put_var(iu%ncid, iu%capm_varid, capm * RAD2DEG, start=[idslot, tslot]), "netcdf_write_frame_base nf90_put_var capm_varid"  ) 
               end if

               select type(self)  
               class is (swiftest_pl)  ! Additional output if the passed polymorphic object is a massive body
                  call check( nf90_put_var(iu%ncid, iu%Gmass_varid, self%Gmass(j), start=[idslot, tslot]), "netcdf_write_frame_base nf90_put_var Gmass_varid"  )
                  if (param%lrhill_present) then
                     call check( nf90_put_var(iu%ncid, iu%rhill_varid, self%rhill(j), start=[idslot, tslot]), "netcdf_write_frame_base nf90_put_var rhill_varid"  )
                  end if
                  if (param%lclose) call check( nf90_put_var(iu%ncid, iu%radius_varid, self%radius(j), start=[idslot, tslot]), "netcdf_write_frame_base nf90_put_var radius_varid"  )
                  if (param%lrotation) then
                     call check( nf90_put_var(iu%ncid, iu%Ip1_varid, self%Ip(1, j), start=[idslot, tslot]), "netcdf_write_frame_base nf90_put_var Ip1_varid"  )
                     call check( nf90_put_var(iu%ncid, iu%Ip2_varid, self%Ip(2, j), start=[idslot, tslot]), "netcdf_write_frame_base nf90_put_var Ip2_varid"  )
                     call check( nf90_put_var(iu%ncid, iu%Ip3_varid, self%Ip(3, j), start=[idslot, tslot]), "netcdf_write_frame_base nf90_put_var Ip3_varid"  )
                     call check( nf90_put_var(iu%ncid, iu%rotx_varid, self%rot(1, j), start=[idslot, tslot]), "netcdf_write_frame_base nf90_put_var rotx_varid"  )
                     call check( nf90_put_var(iu%ncid, iu%roty_varid, self%rot(2, j), start=[idslot, tslot]), "netcdf_write_frame_base nf90_put_var roty_varid"  )
                     call check( nf90_put_var(iu%ncid, iu%rotz_varid, self%rot(3, j), start=[idslot, tslot]), "netcdf_write_frame_base nf90_put_var rotz_varid"  )
                  end if
                  ! if (param%ltides) then
                  !    call check( nf90_put_var(iu%ncid, iu%k2_varid, self%k2(j), start=[idslot, tslot]), "netcdf_write_frame_base nf90_put_var k2_varid"  )
                  !    call check( nf90_put_var(iu%ncid, iu%Q_varid, self%Q(j), start=[idslot, tslot]), "netcdf_write_frame_base nf90_put_var Q_varid"  )
                  ! end if

               end select
            end do
         end associate
      class is (swiftest_cb)
         idslot = self%id + 1
         call check( nf90_put_var(iu%ncid, iu%id_varid, self%id, start=[idslot]), "netcdf_write_frame_base nf90_put_var cb id_varid"  )

         call check( nf90_put_var(iu%ncid, iu%Gmass_varid, self%Gmass, start=[idslot, tslot]), "netcdf_write_frame_base nf90_put_var cb Gmass_varid"  )
         if (param%lclose) call check( nf90_put_var(iu%ncid, iu%radius_varid, self%radius, start=[idslot, tslot]), "netcdf_write_frame_base nf90_put_var cb radius_varid"  )
         call check( nf90_put_var(iu%ncid, iu%j2rp2_varid, self%j2rp2, start=[tslot]), "netcdf_write_frame_base nf90_put_var cb j2rp2_varid" )
         call check( nf90_put_var(iu%ncid, iu%j4rp4_varid, self%j4rp4, start=[tslot]), "netcdf_write_frame_base nf90_put_var cb j4rp4_varid" )
         if (param%lrotation) then
            call check( nf90_put_var(iu%ncid, iu%Ip1_varid, self%Ip(1), start=[idslot, tslot]), "netcdf_write_frame_base nf90_put_var cb Ip1_varid"  )
            call check( nf90_put_var(iu%ncid, iu%Ip2_varid, self%Ip(2), start=[idslot, tslot]), "netcdf_write_frame_base nf90_put_var cb Ip2_varid"  )
            call check( nf90_put_var(iu%ncid, iu%Ip3_varid, self%Ip(3), start=[idslot, tslot]), "netcdf_write_frame_base nf90_put_var cb Ip3_varid"  )
            call check( nf90_put_var(iu%ncid, iu%rotx_varid, self%rot(1), start=[idslot, tslot]), "netcdf_write_frame_base nf90_put_var cb rotx_varid"  )
            call check( nf90_put_var(iu%ncid, iu%roty_varid, self%rot(2), start=[idslot, tslot]), "netcdf_write_frame_base nf90_put_var cb roty_varid"  )
            call check( nf90_put_var(iu%ncid, iu%rotz_varid, self%rot(3), start=[idslot, tslot]), "netcdf_write_frame_base nf90_put_var cb rotz_varid"  )
         end if
         ! if (param%ltides) then
         !    call check( nf90_put_var(iu%ncid, iu%k2_varid, self%k2, start=[idslot, tslot]), "netcdf_write_frame_base nf90_put_var cb k2_varid"  )
         !    call check( nf90_put_var(iu%ncid, iu%Q_varid, self%Q, start=[idslot, tslot]), "netcdf_write_frame_base nf90_put_var cb Q_varid"  )
         ! end if

      end select
      call check( nf90_set_fill(iu%ncid, old_mode, old_mode), "netcdf_write_frame_base nf90_set_fill old_mode"  )

      return
   end subroutine netcdf_write_frame_base


   module subroutine netcdf_write_frame_system(self, iu, param)
      !! author: The Purdue Swiftest Team - David A. Minton, Carlisle A. Wishard, Jennifer L.L. Pouplin, and Jacob R. Elliott
      !!
      !! Write a frame (header plus records for each massive body and active test particle) to a output binary file
      implicit none
      ! Arguments
      class(swiftest_nbody_system), intent(inout) :: self  !! Swiftest system object
      class(netcdf_parameters),     intent(inout) :: iu    !! Parameters used to identify a particular NetCDF dataset
      class(swiftest_parameters),   intent(inout) :: param !! Current run configuration parameters 

      call self%write_hdr(iu, param)
      call self%cb%write_frame(iu, param)
      call self%pl%write_frame(iu, param)
      call self%tp%write_frame(iu, param)

      return
   end subroutine netcdf_write_frame_system


   module subroutine netcdf_write_particle_info_base(self, iu, param)
      !! author: Carlisle A. Wishard, Dana Singh, and David A. Minton
      !!
      !! Write all current particle to file
      implicit none
      ! Arguments
      class(swiftest_base),       intent(in)    :: self   !! Swiftest particle object
      class(netcdf_parameters),   intent(inout) :: iu     !! Parameters used to identify a particular NetCDF dataset
      class(swiftest_parameters), intent(inout) :: param  !! Current run configuration parameters
      ! Internals
      integer(I4B)                              :: i, j, idslot, old_mode
      integer(I4B), dimension(:), allocatable   :: ind
      character(len=NAMELEN)                    :: charstring

      ! This string of spaces of length NAMELEN is used to clear out any old data left behind inside the string variables
      call check( nf90_set_fill(iu%ncid, nf90_nofill, old_mode), "netcdf_write_particle_info_base nf90_set_fill nf90_nofill"  )

      select type(self)
         class is (swiftest_body)
         associate(n => self%nbody)
            if (n == 0) return
            call util_sort(self%id(1:n), ind)

            do i = 1, n
               j = ind(i)
               idslot = self%id(j) + 1
               call check( nf90_put_var(iu%ncid, iu%id_varid, self%id(j), start=[idslot]), "netcdf_write_particle_info_base nf90_put_var id_varid"  )

               charstring = trim(adjustl(self%info(j)%name))
               call check( nf90_put_var(iu%ncid, iu%name_varid, charstring, start=[1, idslot], count=[NAMELEN, 1]), "netcdf_write_particle_info_base nf90_put_var name_varid"  )

               charstring = trim(adjustl(self%info(j)%particle_type))
               call check( nf90_put_var(iu%ncid, iu%ptype_varid, charstring, start=[1, idslot], count=[NAMELEN, 1]), "netcdf_write_particle_info_base nf90_put_var particle_type_varid"  )

               charstring = trim(adjustl(self%info(j)%status))
               call check( nf90_put_var(iu%ncid, iu%status_varid, charstring, start=[1, idslot], count=[NAMELEN, 1]), "netcdf_write_particle_info_base nf90_put_var status_varid"  )

               if (param%lclose) then
                  charstring = trim(adjustl(self%info(j)%origin_type))
                  call check( nf90_put_var(iu%ncid, iu%origin_type_varid, charstring, start=[1, idslot], count=[NAMELEN, 1]), "netcdf_write_particle_info_base nf90_put_var origin_type_varid"  )
                  call check( nf90_put_var(iu%ncid, iu%origin_time_varid,  self%info(j)%origin_time,  start=[idslot]), "netcdf_write_particle_info_base nf90_put_var origin_time_varid"  )
                  call check( nf90_put_var(iu%ncid, iu%origin_xhx_varid,   self%info(j)%origin_xh(1), start=[idslot]), "netcdf_write_particle_info_base nf90_put_var origin_xhx_varid"  )
                  call check( nf90_put_var(iu%ncid, iu%origin_xhy_varid,   self%info(j)%origin_xh(2), start=[idslot]), "netcdf_write_particle_info_base nf90_put_var origin_xhy_varid"  )
                  call check( nf90_put_var(iu%ncid, iu%origin_xhz_varid,   self%info(j)%origin_xh(3), start=[idslot]), "netcdf_write_particle_info_base nf90_put_var origin_xhz_varid"  )
                  call check( nf90_put_var(iu%ncid, iu%origin_vhx_varid,   self%info(j)%origin_vh(1), start=[idslot]), "netcdf_write_particle_info_base nf90_put_var origin_vhx_varid"  )
                  call check( nf90_put_var(iu%ncid, iu%origin_vhy_varid,   self%info(j)%origin_vh(2), start=[idslot]), "netcdf_write_particle_info_base nf90_put_var origin_vhy_varid"  )
                  call check( nf90_put_var(iu%ncid, iu%origin_vhz_varid,   self%info(j)%origin_vh(3), start=[idslot]), "netcdf_write_particle_info_base nf90_put_var origin_vhz_varid"  )
   
                  call check( nf90_put_var(iu%ncid, iu%collision_id_varid, self%info(j)%collision_id, start=[idslot]), "netcdf_write_particle_info_base nf90_put_var collision_id_varid"  )
                  call check( nf90_put_var(iu%ncid, iu%discard_time_varid, self%info(j)%discard_time, start=[idslot]), "netcdf_write_particle_info_base nf90_put_var discard_time_varid"  )
                  call check( nf90_put_var(iu%ncid, iu%discard_xhx_varid, self%info(j)%discard_xh(1), start=[idslot]), "netcdf_write_particle_info_base nf90_put_var discard_xhx_varid"  )
                  call check( nf90_put_var(iu%ncid, iu%discard_xhy_varid, self%info(j)%discard_xh(2), start=[idslot]), "netcdf_write_particle_info_base nf90_put_var discard_xhy_varid"  )
                  call check( nf90_put_var(iu%ncid, iu%discard_xhz_varid, self%info(j)%discard_xh(3), start=[idslot]), "netcdf_write_particle_info_base nf90_put_var discard_xhz_varid"  )
                  call check( nf90_put_var(iu%ncid, iu%discard_vhx_varid, self%info(j)%discard_vh(1), start=[idslot]), "netcdf_write_particle_info_base nf90_put_var discard_vhx_varid"  )
                  call check( nf90_put_var(iu%ncid, iu%discard_vhy_varid, self%info(j)%discard_vh(2), start=[idslot]), "netcdf_write_particle_info_base nf90_put_var discard_vhy_varid"  )
                  call check( nf90_put_var(iu%ncid, iu%discard_vhz_varid, self%info(j)%discard_vh(3), start=[idslot]), "netcdf_write_particle_info_base nf90_put_var discard_vhz_varid"  )
               end if

            end do
         end associate

      class is (swiftest_cb)
         idslot = self%id + 1
         call check( nf90_put_var(iu%ncid, iu%id_varid, self%id, start=[idslot]), "netcdf_write_particle_info_base nf90_put_var cb id_varid"  )

         charstring = trim(adjustl(self%info%name))
         call check( nf90_put_var(iu%ncid, iu%name_varid, charstring, start=[1, idslot], count=[NAMELEN, 1]), "netcdf_write_particle_info_base nf90_put_var cb name_varid"  )

         charstring = trim(adjustl(self%info%particle_type))
         call check( nf90_put_var(iu%ncid, iu%ptype_varid, charstring, start=[1, idslot], count=[NAMELEN, 1]), "netcdf_write_particle_info_base nf90_put_var cb ptype_varid"  )

         charstring = trim(adjustl(self%info%status))
         call check( nf90_put_var(iu%ncid, iu%status_varid, charstring, start=[1, idslot], count=[NAMELEN, 1]), "netcdf_write_particle_info_base nf90_put_var cb status_varid"  )

         if (param%lclose) then
            charstring = trim(adjustl(self%info%origin_type))
            call check( nf90_put_var(iu%ncid, iu%origin_type_varid, charstring, start=[1, idslot], count=[NAMELEN, 1]), "netcdf_write_particle_info_base nf90_put_var cb origin_type_varid"  )

            call check( nf90_put_var(iu%ncid, iu%origin_time_varid, self%info%origin_time, start=[idslot]), "netcdf_write_particle_info_base nf90_put_var cb origin_time_varid"  )
            call check( nf90_put_var(iu%ncid, iu%origin_xhx_varid, self%info%origin_xh(1), start=[idslot]), "netcdf_write_particle_info_base nf90_put_var cb origin_xhx_varid"  )
            call check( nf90_put_var(iu%ncid, iu%origin_xhy_varid, self%info%origin_xh(2), start=[idslot]), "netcdf_write_particle_info_base nf90_put_var cb origin_xhy_varid"  )
            call check( nf90_put_var(iu%ncid, iu%origin_xhz_varid, self%info%origin_xh(3), start=[idslot]), "netcdf_write_particle_info_base nf90_put_var cb origin_xhz_varid"  )
            call check( nf90_put_var(iu%ncid, iu%origin_vhx_varid, self%info%origin_vh(1), start=[idslot]), "netcdf_write_particle_info_base nf90_put_var cb origin_vhx_varid"  )
            call check( nf90_put_var(iu%ncid, iu%origin_vhy_varid, self%info%origin_vh(2), start=[idslot]), "netcdf_write_particle_info_base nf90_put_var cb origin_vhy_varid"  )
            call check( nf90_put_var(iu%ncid, iu%origin_vhz_varid, self%info%origin_vh(3), start=[idslot]), "netcdf_write_particle_info_base nf90_put_var cb origin_vhz_varid"  )
   
            call check( nf90_put_var(iu%ncid, iu%collision_id_varid, self%info%collision_id, start=[idslot]), "netcdf_write_particle_info_base nf90_put_var cb collision_id_varid"  )
            call check( nf90_put_var(iu%ncid, iu%discard_time_varid, self%info%discard_time, start=[idslot]), "netcdf_write_particle_info_base nf90_put_var cb discard_time_varid"  )
            call check( nf90_put_var(iu%ncid, iu%discard_xhx_varid, self%info%discard_xh(1), start=[idslot]), "netcdf_write_particle_info_base nf90_put_var cb discard_xhx_varid"  )
            call check( nf90_put_var(iu%ncid, iu%discard_xhy_varid, self%info%discard_xh(2), start=[idslot]), "netcdf_write_particle_info_base nf90_put_var cb discard_xhy_varid"  )
            call check( nf90_put_var(iu%ncid, iu%discard_xhz_varid, self%info%discard_xh(3), start=[idslot]), "netcdf_write_particle_info_base nf90_put_var cb discard_xhz_varid"  )
            call check( nf90_put_var(iu%ncid, iu%discard_vhx_varid, self%info%discard_vh(1), start=[idslot]), "netcdf_write_particle_info_base nf90_put_var cb discard_vhx_varid"  )
            call check( nf90_put_var(iu%ncid, iu%discard_vhy_varid, self%info%discard_vh(2), start=[idslot]), "netcdf_write_particle_info_base nf90_put_var cb discard_vhy_varid"  )
            call check( nf90_put_var(iu%ncid, iu%discard_vhz_varid, self%info%discard_vh(3), start=[idslot]), "netcdf_write_particle_info_base nf90_put_var cb discard_vhz_varid"  )
         end if

      end select

      call check( nf90_set_fill(iu%ncid, old_mode, old_mode) )
      return
   end subroutine netcdf_write_particle_info_base


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
      class(swiftest_parameters),   intent(inout) :: param !! Current run configuration parameters
      ! Internals
      integer(I4B) :: tslot

      tslot = int(param%ioutput, kind=I4B) + 1

      call check( nf90_put_var(iu%ncid, iu%time_varid, param%t, start=[tslot]), "netcdf_write_hdr_system nf90_put_var time_varid"  )
      call check( nf90_put_var(iu%ncid, iu%npl_varid, self%pl%nbody, start=[tslot]), "netcdf_write_hdr_system nf90_put_var npl_varid"  )
      call check( nf90_put_var(iu%ncid, iu%ntp_varid, self%tp%nbody, start=[tslot]), "netcdf_write_hdr_system nf90_put_var ntp_varid"  )
      select type(pl => self%pl)
      class is (symba_pl)
         call check( nf90_put_var(iu%ncid, iu%nplm_varid, pl%nplm, start=[tslot]), "netcdf_write_hdr_system nf90_put_var nplm_varid"  )
      end select

      if (param%lenergy) then
         call check( nf90_put_var(iu%ncid, iu%KE_orb_varid, self%ke_orbit, start=[tslot]), "netcdf_write_hdr_system nf90_put_var KE_orb_varid"  )
         call check( nf90_put_var(iu%ncid, iu%KE_spin_varid, self%ke_spin, start=[tslot]), "netcdf_write_hdr_system nf90_put_var KE_spin_varid"  )
         call check( nf90_put_var(iu%ncid, iu%PE_varid, self%pe, start=[tslot]), "netcdf_write_hdr_system nf90_put_var PE_varid"  )
         call check( nf90_put_var(iu%ncid, iu%L_orbx_varid, self%Lorbit(1), start=[tslot]), "netcdf_write_hdr_system nf90_put_var L_orbx_varid"  )
         call check( nf90_put_var(iu%ncid, iu%L_orby_varid, self%Lorbit(2), start=[tslot]), "netcdf_write_hdr_system nf90_put_var L_orby_varid"  )
         call check( nf90_put_var(iu%ncid, iu%L_orbz_varid, self%Lorbit(3), start=[tslot]), "netcdf_write_hdr_system nf90_put_var L_orbz_varid"  )
         call check( nf90_put_var(iu%ncid, iu%L_spinx_varid, self%Lspin(1), start=[tslot]), "netcdf_write_hdr_system nf90_put_var L_spinx_varid"  )
         call check( nf90_put_var(iu%ncid, iu%L_spiny_varid, self%Lspin(2), start=[tslot]), "netcdf_write_hdr_system nf90_put_var L_spiny_varid"  )
         call check( nf90_put_var(iu%ncid, iu%L_spinz_varid, self%Lspin(3), start=[tslot]), "netcdf_write_hdr_system nf90_put_var L_spinz_varid"  )
         call check( nf90_put_var(iu%ncid, iu%L_escapex_varid, self%Lescape(1), start=[tslot]), "netcdf_write_hdr_system nf90_put_var L_escapex_varid"  )
         call check( nf90_put_var(iu%ncid, iu%L_escapey_varid, self%Lescape(2), start=[tslot]), "netcdf_write_hdr_system nf90_put_var L_escapey_varid"  )
         call check( nf90_put_var(iu%ncid, iu%L_escapez_varid, self%Lescape(3), start=[tslot]), "netcdf_write_hdr_system nf90_put_var L_escapez_varid"  )
         call check( nf90_put_var(iu%ncid, iu%Ecollisions_varid, self%Ecollisions, start=[tslot]), "netcdf_write_hdr_system nf90_put_var Ecollisions_varid"  )
         call check( nf90_put_var(iu%ncid, iu%Euntracked_varid, self%Euntracked, start=[tslot]), "netcdf_write_hdr_system nf90_put_var Euntracked_varid"  )
         call check( nf90_put_var(iu%ncid, iu%GMescape_varid, self%GMescape, start=[tslot]), "netcdf_write_hdr_system nf90_put_var GMescape_varid"  )
      end if

      return
   end subroutine netcdf_write_hdr_system

end submodule s_netcdf
