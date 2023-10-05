!! Copyright 2022 - David Minton, Carlisle Wishard, Jennifer Pouplin, Jake Elliott, & Dana Singh
!! This file is part of Swiftest.
!! Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
!! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!! Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
!! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!! You should have received a copy of the GNU General Public License along with Swiftest. 
!! If not, see: https://www.gnu.org/licenses. 

submodule (swiftest) s_swiftest_io
   use symba
   use netcdf
contains

   module subroutine swiftest_io_compact_output(self, param, timer)
      !! author: David Minton
      !!
      !! Generates the terminal output displayed when display_style is set to COMPACT. This is used by the Python driver to 
      !! make nice-looking progress reports.
      implicit none

      interface fmt
         !! author: David Minton
         !!
         !! Formats a pair of variables and corresponding values for the compact display output. Generic interface for different 
         !! variable types to format.
         procedure :: fmt_I4B, fmt_I8B, fmt_DP
      end interface

      ! Arguments
      class(swiftest_nbody_system), intent(in) :: self  !! Swiftest nbody system object   
      class(swiftest_parameters),   intent(in) :: param !! Input colleciton of user-defined parameters
      class(*),                     intent(in) :: timer !! Object used for computing elapsed wall time  (must be unlimited 
                                                        !! polymorphic because the walltimer module requires base)
      ! Internals
      character(len=:), allocatable :: formatted_output

      select type(timer)
      class is (walltimer)
         formatted_output = fmt("ILOOP",param%iloop) // fmt("T",self%t) // fmt("NPL",self%pl%nbody) // fmt("NTP",self%tp%nbody) 
         select type(pl => self%pl)
         class is (symba_pl)
            formatted_output = formatted_output // fmt("NPLM",pl%nplm)
         end select
         if (param%lenergy) then
            formatted_output = formatted_output // fmt("LTOTERR",self%L_total_error) // fmt("ETOTERR",self%te_error) &
                             // fmt("MTOTERR",self%Mtot_error) // fmt("KEOERR",self%ke_orbit_error) // fmt("PEERR",self%pe_error) &
                             // fmt("EORBERR",self%E_orbit_error) // fmt("EUNTRERR",self%E_untracked_error) &
                             // fmt("LESCERR",self%L_escape_error) // fmt("MESCERR",self%Mescape_error)
            if (param%lclose) formatted_output = formatted_output // fmt("ECOLLERR",self%Ecoll_error)
            if (param%lrotation) formatted_output = formatted_output // fmt("KESPINERR",self%ke_spin_error) &
                             // fmt("LSPINERR",self%L_spin_error) 
         end if

         if (.not. timer%main_is_started) then ! This is the start of a new run
            formatted_output =  formatted_output // fmt("WT",0.0_DP) // fmt("IWT",0.0_DP) // fmt("WTPS",0.0_DP) 
         else
            formatted_output = formatted_output // fmt("WT",timer%wall_main) // fmt("IWT",timer%wall_step) &
                                               // fmt("WTPS",timer%wall_per_substep)
         end if
         write(*,*) formatted_output
      end select
      return

      contains

         function fmt_I4B(varname,val) result(pair_string)
            implicit none
            ! Arguments
            character(*), intent(in) :: varname !! The variable name of the pair
            integer(I4B), intent(in) :: val !! A 4-byte integer value
            ! Result
            character(len=:), allocatable :: pair_string
            ! Internals
            character(len=24) :: str_value
      
            write(str_value,*) val
            pair_string = trim(adjustl(varname)) // " " // trim(adjustl(str_value)) // ";"

            return 
         end function fmt_I4B

         function fmt_I8B(varname, val) result(pair_string)
            implicit none
            ! Arguments
            character(*), intent(in) :: varname !! The variable name of the pair
            integer(I8B), intent(in) :: val     !! An 8-byte integer value
            ! Result
            character(len=:), allocatable :: pair_string
            ! Internals
            character(len=24) :: str_value
      
            write(str_value,*) val
            pair_string = trim(adjustl(varname)) // " " // trim(adjustl(str_value)) // ";"

            return 
         end function fmt_I8B

         function fmt_DP(varname, val) result(pair_string)
            implicit none
            ! Arguments
            character(*), intent(in) :: varname !! The variable name of the pair
            real(DP),     intent(in) :: val     !! A double precision floating point value
            ! Result
            character(len=:), allocatable :: pair_string
            ! Internals
            character(len=24) :: str_value
      
            write(str_value,'(ES24.16)') val
            pair_string = trim(adjustl(varname)) // " " // trim(adjustl(str_value)) // ";"

            return 
         end function fmt_DP

   end subroutine swiftest_io_compact_output


   module subroutine swiftest_io_conservation_report(self, param, lterminal)
      !! author: The Purdue Swiftest Team -  David A. Minton, Carlisle A. Wishard, Jennifer L.L. Pouplin, and Jacob R. Elliott
      !!
      !! Reports the current state of energy, mass, and angular momentum conservation in a run
      implicit none
      ! Arguments
      class(swiftest_nbody_system), intent(inout) :: self      !! Swiftest nbody system object
      class(swiftest_parameters),   intent(inout) :: param     !! Input colleciton of user-defined parameters
      logical,                      intent(in)    :: lterminal !! Indicates whether to output information to the terminal screen
      ! Internals
      real(DP), dimension(NDIM)       :: L_total_now,  L_orbit_now,  L_spin_now
      real(DP)                        :: ke_orbit_now,  ke_spin_now,  pe_now,  E_orbit_now, be_now, be_cb_now, be_cb_orig, te_now
      real(DP)                        :: GMtot_now
      character(len=*), parameter     :: EGYTERMFMT = '(" DL/L0 = ", ES12.5, "; DE_orbit/|E0| = ", ES12.5,' &
                                                     //'"; DE_total/|E0| = ", ES12.5, "; DM/M0 = ", ES12.5)'

      associate(nbody_system => self, pl => self%pl, cb => self%cb, npl => self%pl%nbody, display_unit => param%display_unit)

         select type(self)
         class is (helio_nbody_system) ! Don't convert vh to vb for Helio-based integrators, because it's already calculated
         class is (whm_nbody_system)
            call pl%vh2vb(cb)
         end select
         call pl%rh2rb(cb)

         call nbody_system%get_energy_and_momentum(param) 
         ke_orbit_now = nbody_system%ke_orbit
         ke_spin_now = nbody_system%ke_spin
         pe_now = nbody_system%pe
         be_now = nbody_system%be
         be_cb_now = nbody_system%be_cb
         te_now = nbody_system%te
         L_orbit_now(:) = nbody_system%L_orbit(:)
         L_spin_now(:) = nbody_system%L_spin(:)
         E_orbit_now = ke_orbit_now + pe_now
         L_total_now(:) = nbody_system%L_total(:) + nbody_system%L_escape(:)
         GMtot_now = nbody_system%GMtot + nbody_system%GMescape 

         if (param%lfirstenergy) then
            nbody_system%ke_orbit_orig = ke_orbit_now
            nbody_system%ke_spin_orig = ke_spin_now
            nbody_system%pe_orig = pe_now
            nbody_system%be_orig = be_now
            nbody_system%te_orig = te_now
            nbody_system%E_orbit_orig = E_orbit_now
            nbody_system%GMtot_orig = GMtot_now
            nbody_system%L_orbit_orig(:) = L_orbit_now(:)
            nbody_system%L_spin_orig(:) = L_spin_now(:)
            nbody_system%L_total_orig(:) = L_total_now(:)
            param%lfirstenergy = .false.
         end if

         if (.not.param%lfirstenergy) then 

            nbody_system%ke_orbit_error = (ke_orbit_now - nbody_system%ke_orbit_orig) / abs(nbody_system%E_orbit_orig)
            nbody_system%ke_spin_error = (ke_spin_now - nbody_system%ke_spin_orig) / abs(nbody_system%E_orbit_orig)
            nbody_system%pe_error = (pe_now - nbody_system%pe_orig) / abs(nbody_system%E_orbit_orig)

            be_cb_orig = -(3 * cb%GM0**2 / param%GU) / (5 * cb%R0)
            nbody_system%be_error = (be_now - nbody_system%be_orig) / abs(nbody_system%te_orig) + (be_cb_now - be_cb_orig) & 
                                    / abs(nbody_system%te_orig)

            nbody_system%E_orbit_error = (E_orbit_now - nbody_system%E_orbit_orig) / abs(nbody_system%E_orbit_orig)
            nbody_system%Ecoll_error = nbody_system%E_collisions / abs(nbody_system%te_orig)
            nbody_system%E_untracked_error = nbody_system%E_untracked / abs(nbody_system%te_orig)
            nbody_system%te_error = (nbody_system%te - nbody_system%te_orig - nbody_system%E_collisions - nbody_system%E_untracked)&
                                   / abs(nbody_system%te_orig) + (be_cb_now - be_cb_orig) / abs(nbody_system%te_orig)

            nbody_system%L_orbit_error = norm2(L_orbit_now(:) - nbody_system%L_orbit_orig(:)) / norm2(nbody_system%L_total_orig(:))
            nbody_system%L_spin_error = norm2(L_spin_now(:) - nbody_system%L_spin_orig(:)) / norm2(nbody_system%L_total_orig(:))
            nbody_system%L_escape_error = norm2(nbody_system%L_escape(:)) / norm2(nbody_system%L_total_orig(:))
            nbody_system%L_total_error = norm2(L_total_now(:) - nbody_system%L_total_orig(:)) / norm2(nbody_system%L_total_orig(:))

            nbody_system%Mescape_error = nbody_system%GMescape / nbody_system%GMtot_orig
#ifdef COARRAY
   if (this_image() == 1 .or. param%log_output) then
#endif 
            if (lterminal) then
               write(display_unit, EGYTERMFMT) nbody_system%L_total_error, nbody_system%E_orbit_error, nbody_system%te_error, &
                                               nbody_system%Mtot_error
               if (param%log_output) flush(display_unit)
            end if

#ifdef COARRAY
   end if ! (this_image() == 1) then
#endif 

            if (abs(nbody_system%Mtot_error) > 100 * epsilon(nbody_system%Mtot_error)) then
               write(*,*) "Severe error! Mass not conserved! Halting!"
               ! Save the frame of data to the bin file in the slot just after the present one for diagnostics
               call base_util_exit(FAILURE,param%display_unit)
            end if
         end if
      end associate

      return
   end subroutine swiftest_io_conservation_report


   module subroutine swiftest_io_display_run_information(self, param, integration_timer, phase)
      !! author:  David A. Minton
      !!
      !! Displays helpful information while a run is executing. The format of the output depends on user parameters
      implicit none
      ! Arguments
      class(swiftest_nbody_system), intent(inout) :: self !! Swiftest nbody system object
      class(swiftest_parameters),   intent(inout) :: param !! Current run configuration parameters 
      type(walltimer),              intent(inout) :: integration_timer !! Object used for computing elapsed wall time
      character(len=*), optional,   intent(in)    :: phase !! One of "first" or "last" 
      ! Internals
      integer(I4B)                                   :: phase_val   
      real(DP)                                       :: tfrac             !! Fraction of total simulation time completed
      type(progress_bar), save                       :: pbar              !! Object used to print out a progress bar
      character(len=64)                              :: pbarmessage
#ifdef COARRAY
      character(*), parameter :: co_statusfmt = '("Image: ",I4, "; Time = ", ES12.5, "; fraction done = ", F6.3, ' // & 
                                             '"; Number of active pl, tp = ", I6, ", ", I6)'
      character(*), parameter :: co_symbastatfmt = '("Image: ",I4, "; Image: Time = ", ES12.5, "; fraction done = ", F6.3, ' // &
                                                '"; Number of active pl, plm, tp = ", I6, ", ", I6, ", ", I6)'
#endif
      character(*), parameter :: statusfmt = '("Time = ", ES12.5, "; fraction done = ", F6.3, ' // & 
                                             '"; Number of active pl, tp = ", I6, ", ", I6)'
      character(*), parameter :: symbastatfmt = '("Time = ", ES12.5, "; fraction done = ", F6.3, ' // &
                                                '"; Number of active pl, plm, tp = ", I6, ", ", I6, ", ", I6)'
      character(*), parameter :: pbarfmt = '("Time = ", ES12.5," of ",ES12.5)'

! The following will syncronize the images so that they report in order, and only write to file one at at ime


      phase_val = 1
      if (present(phase)) then
         if (phase == "first") then
            phase_val = 0
         else if (phase == "last") then
            phase_val = -1
         end if
      end if

      tfrac = (self%t - param%t0) / (param%tstop - param%t0)

#ifdef COARRAY
         if (this_image() == 1 .or. param%log_output) then
#endif
         if (phase_val == 0) then
            if (param%lrestart) then
               write(param%display_unit, *) " *************** Swiftest " // VERSION // " restart " // &
                                               trim(adjustl(param%integrator)) // " *************** "
            else
               write(param%display_unit, *) " *************** Swiftest " // VERSION // " start " // &
                                              trim(adjustl(param%integrator)) // " *************** "
            end if
            if (param%display_style == "PROGRESS") then
               call pbar%reset(param%nloops)
            else if (param%display_style == "COMPACT") then
               write(*,*) "SWIFTEST START " // trim(adjustl(param%integrator))
            end if
         end if
#ifdef COARRAY
         end if !(this_image() == 1)
#endif

      if (param%display_style == "PROGRESS") then
#ifdef COARRAY
         if (this_image() == 1) then
#endif
            write(pbarmessage,fmt=pbarfmt) self%t, param%tstop
            call pbar%update(param%iloop,message=pbarmessage)
#ifdef COARRAY
         end if !(this_image() == 1)
#endif
      else if (param%display_style == "COMPACT") then
         call self%compact_output(param,integration_timer)
      end if

      if (param%lmtiny_pl) then
#ifdef COARRAY
         if (param%lcoarray) then
            write(param%display_unit, co_symbastatfmt) this_image(),self%t, tfrac, self%pl%nbody, self%pl%nplm, self%tp%nbody
         else
#endif
            write(param%display_unit, symbastatfmt) self%t, tfrac, self%pl%nbody, self%pl%nplm, self%tp%nbody
#ifdef COARRAY
         end if
#endif
      else
#ifdef COARRAY
         if (param%lcoarray) then
            write(param%display_unit, co_statusfmt) this_image(),self%t, tfrac, self%pl%nbody, self%tp%nbody
         else
#endif
            write(param%display_unit, statusfmt) self%t, tfrac, self%pl%nbody, self%tp%nbody
#ifdef COARRAY
         end if
#endif
      end if

#ifdef COARRAY
      if (this_image() == num_images() .or. param%log_output) then
#endif
         if (phase_val == -1) then
            write(param%display_unit, *)" *************** Swiftest stop " // trim(adjustl(param%integrator)) // " *************** "
            if (param%display_style == "COMPACT") write(*,*) "SWIFTEST STOP" // trim(adjustl(param%integrator))
            if (param%log_output) close(param%display_unit)
         end if

#ifdef COARRAY
      end if ! this_image() == num_images()
      if (param%log_output) flush(param%display_unit)
#endif

      return
   end subroutine swiftest_io_display_run_information


   module subroutine swiftest_io_dump_param(self, param_file_name)
      !! author: David A. Minton
      !!
      !! Dump integration parameters to file
      !!
      !! Adapted from David E. Kaufmann's Swifter routine io_dump_param.f90
      !! Adapted from Martin Duncan's Swift routine io_dump_param.f
      implicit none
      ! Arguments
      class(swiftest_parameters),intent(in) :: self    !! Output collection of parameters
      character(len=*),          intent(in) :: param_file_name !! Parameter input file name (i.e. param.in)
      ! Internals
      character(STRMAX)        :: errmsg !! Error message in UDIO procedure
      integer(I4B)             :: ierr

      open(unit = LUN, file = param_file_name, status='replace', form = 'FORMATTED', err = 667, iomsg = errmsg)
      !! todo: Currently this procedure does not work in user-defined derived-type input mode 
      !!    due to compiler incompatabilities
      !write(LUN,'(DT)') param
      call self%writer(LUN, iotype = "none", v_list = [0], iostat = ierr, iomsg = errmsg)
      if (ierr == 0) then
         close(LUN, err = 667, iomsg = errmsg)
         return
      end if

      667 continue
      write(*,*) "Error opening parameter dump file " // trim(adjustl(errmsg))
      call base_util_exit(FAILURE,self%display_unit)
   end subroutine swiftest_io_dump_param


   module subroutine swiftest_io_dump_system(self, param, system_history)
      !! author: David A. Minton
      !!
      !! Dumps the state of the nbody_system to files in case the simulation is interrupted.
      !! As a safety mechanism, there are two dump files that are written in alternating order
      !! so that if a dump file gets corrupted during writing, the user can restart from the older one.
      implicit none
      ! Arguments
      class(swiftest_nbody_system), intent(inout) :: self  !! Swiftest nbody_system object
      class(swiftest_parameters),   intent(inout) :: param !! Current run configuration parameters 
      class(swiftest_storage),      intent(inout) :: system_history    !! Stores the system history between output dumps
      ! Internals
      class(swiftest_parameters), allocatable :: param_restart !! Local parameters variable used to parameters change input file 
                                                               !! names to dump file-specific values without changing the 
                                                               !! user-defined values
      character(len=:), allocatable :: param_file_name
      character(len=STRMAX) :: time_text
   
      ! Dump the encounter history if necessary
      if (param%lenc_save_trajectory &
         .or. param%lenc_save_closest &
         .and. allocated(self%encounter_history)) &
            call self%encounter_history%dump(param)
      if (allocated(self%collision_history)) call self%collision_history%dump(param)

      ! Dump the nbody_system history to file
      call system_history%dump(param)

#ifdef COARRAY
      if (this_image() == 1) then
#endif 
         allocate(param_restart, source=param)
         param_restart%in_form  = "XV"
         param_restart%out_stat = 'APPEND'
         param_restart%in_type = "NETCDF_DOUBLE"
         param_restart%nc_in = param%outfile
         param_restart%lrestart = .true.
         param_restart%tstart = self%t
         param_file_name    = trim(adjustl(PARAM_RESTART_FILE))
         call param_restart%dump(param_file_name)
         write(time_text,'(I0.20)') param%iloop
         param_file_name = "param." // trim(adjustl(time_text)) // ".in"
         call param_restart%dump(param_file_name)
#ifdef COARRAY
      end if ! (this_image() == 1) 
#endif 
      return
   end subroutine swiftest_io_dump_system


   module subroutine swiftest_io_dump_storage(self, param)
      !! author: David A. Minton
      !!
      !! Dumps the time history of the simulation to file. Each time it writes a frame to file, it deallocates the nbody_system
      !! object from inside. It will only dump frames with systems that are allocated, so this can be called at the end of
      !! a simulation for cases when the number of saved frames is not equal to the dump cadence (for instance, if the dump
      !! cadence is not divisible by the total number of loops).
      implicit none
      ! Arguments
      class(swiftest_storage),   intent(inout)        :: self   !! Swiftest simulation history storage object
      class(swiftest_parameters),   intent(inout)     :: param  !! Current run configuration parameters 
      ! Internals
      integer(I4B) :: i
#ifdef COARRAY
      type(walltimer) :: iotimer
#endif

      if (self%iframe == 0) return
      call self%make_index_map()
      associate(nc => self%nc)
#ifdef COARRAY
         sync all
         write(param%display_unit,*) "File output started"
         call iotimer%start()
         critical
#endif
         call nc%open(param)
         do i = 1, self%iframe
            ! Writing files is more efficient if we write out the common frames from each image before going to the next frame
            if (allocated(self%frame(i)%item)) then
               select type(nbody_system => self%frame(i)%item)
               class is (swiftest_nbody_system)
                  call nbody_system%write_frame(nc, param)
               end select
               deallocate(self%frame(i)%item)
            end if
         end do
         call nc%close()
#ifdef COARRAY  
         end critical
         call iotimer%stop()
         sync all
         call iotimer%report(message="File output :", unit=param%display_unit)
         flush(param%display_unit)
#endif
      end associate

      call self%reset()
      return
   end subroutine swiftest_io_dump_storage


   module subroutine swiftest_io_get_args(integrator, param_file_name, display_style, from_cli)
      !! author: David A. Minton
      !!
      !! Reads in the name of the parameter file from command line arguments. 
      implicit none
      ! Arguments
      character(len=:), intent(inout), allocatable :: integrator      !! Symbolic code of the requested integrator  
      character(len=:), intent(inout), allocatable :: param_file_name !! Name of the input parameters file
      character(len=:), intent(inout), allocatable :: display_style   !! Style of the output display 
                                                                      !! {"STANDARD", "COMPACT", "PROGRESS"}). Default is "STANDARD"
      logical,          intent(in)                 :: from_cli        !! If true, get command-line arguments. Otherwise, use the 
                                                                      !! values of the input variables
      ! Internals
      character(len=STRMAX), dimension(:), allocatable :: arg
      character(len=STRMAX), dimension(3)              :: tmp_arg
      integer(I4B), dimension(:), allocatable :: ierr
      integer :: i,narg

      if (from_cli) then
         narg = command_argument_count() 
         if (narg > 0) then
            allocate(arg(narg),ierr(narg))
            do i = 1,narg
               call get_command_argument(i, arg(i), status = ierr(i))
            end do
            if (any(ierr /= 0)) call base_util_exit(USAGE)
         else
            call base_util_exit(USAGE)
         end if
      else
         narg = 0
         if (len(integrator) > 0) then
            narg = narg + 1
            tmp_arg(narg) = integrator
         end if
         if (len(param_file_name) > 0) then
            narg = narg + 1
            tmp_arg(narg) = param_file_name
         end if
         if (len(display_style) > 0) then
            narg = narg + 1
            tmp_arg(narg) = display_style
         end if
         allocate(arg(narg))
         arg(1:narg) = tmp_arg(1:narg)
         deallocate(integrator, param_file_name, display_style)
      end if
   
      if (narg == 1) then
         if (arg(1) == '-v' .or. arg(1) == '--version') then
            call swiftest_util_version() 
         else if (arg(1) == '-h' .or. arg(1) == '--help') then
            call base_util_exit(HELP)
         else
            call base_util_exit(USAGE)
         end if
      else if (narg >= 2) then
         call swiftest_io_toupper(arg(1))
         select case(arg(1))
         case('BS')
            integrator = INT_BS
         case('HELIO')
            integrator = INT_HELIO
         case('RA15')
            integrator = INT_RA15
         case('TU4')
            integrator = INT_TU4
         case('WHM')
            integrator = INT_WHM
         case('RMVS')
            integrator = INT_RMVS
         case('SYMBA')
            integrator = INT_SYMBA
         case('RINGMOONS')
            integrator = INT_RINGMOONS
         case default
            integrator = UNKNOWN_INTEGRATOR
            write(*,*) trim(adjustl(arg(1))) // ' is not a valid integrator.'
            call base_util_exit(USAGE)
         end select
         param_file_name = trim(adjustl(arg(2)))
      end if

      if (narg == 2) then
         display_style = "STANDARD"
      else if (narg == 3) then
         call swiftest_io_toupper(arg(3))
         display_style = trim(adjustl(arg(3)))
      else
         call base_util_exit(USAGE)
      end if

      return
   end subroutine swiftest_io_get_args


   module function swiftest_io_get_token(buffer, ifirst, ilast, ierr) result(token)
      !! author: David A. Minton
      !!
      !! Retrieves a character token from an input string. Here a token is defined as any set of contiguous non-blank characters not
      !! beginning with or containing "!". If "!" is present, any remaining part of the buffer including the "!" is ignored
      !!
      !! Adapted from David E. Kaufmann's Swifter routine swiftest_io_get_token.f90
      implicit none
      ! Arguments
      character(len=*), intent(in)    :: buffer         !! Input string buffer
      integer(I4B),     intent(inout) :: ifirst         !! Index of the buffer at which to start the search for a token
      integer(I4B),     intent(out)   :: ilast          !! Index of the buffer at the end of the returned token
      integer(I4B),     intent(out)   :: ierr           !! Error code
      ! Result
      character(len=:), allocatable   :: token          !! Returned token string
      ! Internals
      integer(I4B) :: i,ilength
   
      ilength = len(buffer)
   
      if (ifirst > ilength) then
          ilast = ifirst
          ierr = -1 !! Bad input
          token = ''
          return
      end if
      do i = ifirst, ilength
          if (buffer(i:i) /= ' ') exit
      end do
      if ((i > ilength) .or. (buffer(i:i) == '!')) then
          ifirst = i
          ilast = i
          ierr = -2 !! No valid token
          token = ''
          return
      end if
      ifirst = i
      do i = ifirst, ilength
          if ((buffer(i:i) == ' ') .or. (buffer(i:i) == '!')) exit
      end do
      ilast = i - 1
      ierr = 0
   
      token = buffer(ifirst:ilast)

      return
   end function swiftest_io_get_token


   module subroutine swiftest_io_log_one_message(file, message)
      !! author: David A. Minton
      !!
      !! Writes a single message to a log file
      implicit none
      ! Arguments
      character(len=*), intent(in) :: file   !! Name of file to log
      character(len=*), intent(in) :: message
      ! Internals
      character(STRMAX) :: errmsg

      open(unit=LUN, file=trim(adjustl(file)), status = 'OLD', position = 'APPEND', form = 'FORMATTED', err = 667, iomsg = errmsg)
      write(LUN, *) trim(adjustl(message)) 
      close(LUN)

      return
      667 continue
      write(*,*) "Error writing message to log file: " // trim(adjustl(errmsg))
   end subroutine swiftest_io_log_one_message


   module subroutine swiftest_io_log_start(param, file, header)
      !! author: David A. Minton
      !!
      !! Checks to see if a log file needs to be created if this is a new run, or appended if this is a restarted run
      implicit none
      ! Arguments
      class(swiftest_parameters), intent(in) :: param  !! Current Swiftest run configuration parameters
      character(len=*),           intent(in) :: file   !! Name of file to log
      character(len=*),           intent(in) :: header !! Header to print at top of log file
      ! Internals
      character(STRMAX) :: errmsg
      logical           :: fileExists

      inquire(file=trim(adjustl(file)), exist=fileExists)
      if (.not.param%lrestart .or. .not.fileExists) then
         open(unit=LUN, file=file, status="REPLACE", err = 667, iomsg = errmsg)
         write(LUN, *, err = 667, iomsg = errmsg) trim(adjustl(header))
      end if
      close(LUN)

      return

      667 continue
      write(*,*) "Error writing log file: " // trim(adjustl(errmsg))
   end subroutine swiftest_io_log_start


   module subroutine swiftest_io_netcdf_flush(self, param)
      !! author: David A. Minton
      !!
      !! Flushes the current buffer to disk by closing and re-opening the file.
      !!    
      implicit none
      ! Arguments
      class(swiftest_netcdf_parameters), intent(inout) :: self !! Parameters used to identify a particular NetCDF dataset
      class(swiftest_parameters),        intent(inout) :: param !! Current run configuration parameters 

      call self%close()
      call self%open(param,readonly=.false.)

      return
   end subroutine swiftest_io_netcdf_flush


   module subroutine swiftest_io_netcdf_get_t0_values_system(self, nc, param) 
      !! author: David A. Minton
      !!
      !! Gets the t0 values of various parameters such as energy and momentum
      !!
      implicit none
      ! Arguments
      class(swiftest_nbody_system),      intent(inout) :: self
      class(swiftest_netcdf_parameters), intent(inout) :: nc     !! Parameters used to identify a particular NetCDF dataset
      class(swiftest_parameters),        intent(inout) :: param
      ! Internals
      integer(I4B)                              :: itmax, idmax, tslot
      real(DP), dimension(:), allocatable       :: vals
      logical, dimension(:), allocatable        :: plmask, tpmask
      real(DP), dimension(1)                    :: rtemp
      real(DP), dimension(NDIM)                 :: rot0, Ip0, L
      real(DP) :: mass0

      associate (cb => self%cb)
         call nc%open(param, readonly=.true.)
         call nc%find_tslot(param%t0, tslot)
         call netcdf_io_check( nf90_inquire_dimension(nc%id, nc%time_dimid, len=itmax), "netcdf_io_get_t0_values_system time_dimid")
         call netcdf_io_check( nf90_inquire_dimension(nc%id, nc%name_dimid, len=idmax), "netcdf_io_get_t0_values_system name_dimid")
         allocate(vals(idmax))
         call netcdf_io_check( nf90_get_var(nc%id, nc%time_varid, rtemp, start=[tslot], count=[1]), &
                              "netcdf_io_get_t0_values_system time_varid" )

         if (param%lenergy) then
            call netcdf_io_check( nf90_get_var(nc%id, nc%KE_orb_varid, rtemp, start=[tslot], count=[1]), & 
                                  "netcdf_io_get_t0_values_system KE_orb_varid" )
            self%ke_orbit_orig = rtemp(1)

            call netcdf_io_check( nf90_get_var(nc%id, nc%KE_spin_varid, rtemp, start=[tslot], count=[1]), &
                                 "netcdf_io_get_t0_values_system KE_spin_varid" )
            self%ke_spin_orig = rtemp(1)

            call netcdf_io_check( nf90_get_var(nc%id, nc%PE_varid, rtemp, start=[tslot], count=[1]), &
                                  "netcdf_io_get_t0_values_system PE_varid" )
            self%pe_orig = rtemp(1)

            call netcdf_io_check( nf90_get_var(nc%id, nc%BE_varid, rtemp, start=[tslot], count=[1]), &
                                  "netcdf_io_get_t0_values_system BE_varid" )
            self%be_orig = rtemp(1)
            
            call netcdf_io_check( nf90_get_var(nc%id, nc%TE_varid, rtemp, start=[tslot], count=[1]), &
                                  "netcdf_io_get_t0_values_system TE_varid" )
            self%te_orig = rtemp(1)

            self%E_orbit_orig = self%ke_orbit_orig + self%pe_orig

            call netcdf_io_check( nf90_get_var(nc%id, nc%L_orbit_varid, self%L_orbit_orig(:), start=[1,tslot], count=[NDIM,1]), &
                                  "netcdf_io_get_t0_values_system L_orbit_varid" )
            call netcdf_io_check( nf90_get_var(nc%id, nc%L_spin_varid, self%L_spin_orig(:), start=[1,tslot], count=[NDIM,1]), &
                                  "netcdf_io_get_t0_values_system L_spin_varid" )

            self%L_total_orig(:) = self%L_orbit_orig(:) + self%L_spin_orig(:) 

            call netcdf_io_check( nf90_get_var(nc%id, nc%Gmass_varid, vals, start=[1,tslot], count=[idmax,1]), &
                                  "netcdf_io_get_t0_values_system Gmass_varid" )
            call nc%get_valid_masks(plmask,tpmask)
            self%GMtot_orig = vals(1) + sum(vals(2:idmax), plmask(:))

            cb%GM0 = vals(1)
            cb%dGM = cb%Gmass - cb%GM0
            mass0 = cb%GM0 / param%GU

            call netcdf_io_check( nf90_get_var(nc%id, nc%radius_varid, rtemp, start=[1,tslot], count=[1,1]), &
                                  "netcdf_io_get_t0_values_system radius_varid" )
            cb%R0 = rtemp(1) 

            if (param%lrotation) then
               call netcdf_io_check( nf90_get_var(nc%id, nc%rot_varid, rot0, start=[1,1,tslot], count=[NDIM,1,1]), &
                                     "netcdf_io_get_t0_values_system rot_varid" )
               rot0(:) = rot0(:) * DEG2RAD
               call netcdf_io_check( nf90_get_var(nc%id, nc%Ip_varid, Ip0, start=[1,1,tslot], count=[NDIM,1,1]), &
                                     "netcdf_io_get_t0_values_system Ip_varid" )
               cb%L0(:) = Ip0(3) * mass0 * cb%R0**2 * rot0(:)
               L(:) = cb%Ip(3) * cb%mass * cb%radius**2 * cb%rot(:)
               cb%dL(:) = L(:) - cb%L0
            end if

            ! Retrieve the current bookkeeping variables
            call nc%find_tslot(self%t, tslot)
            call netcdf_io_check( nf90_get_var(nc%id, nc%L_escape_varid, self%L_escape(:),  start=[1,tslot], count=[NDIM,1]), &
                                  "netcdf_io_get_t0_values_system L_escape_varid" )
            call netcdf_io_check( nf90_get_var(nc%id, nc%GMescape_varid,    self%GMescape,    start=[tslot]), &
                                  "netcdf_io_get_t0_values_system GMescape_varid" )
            call netcdf_io_check( nf90_get_var(nc%id, nc%E_collisions_varid, self%E_collisions, start=[tslot]), &
                                  "netcdf_io_get_t0_values_system E_collisions_varid" )
            call netcdf_io_check( nf90_get_var(nc%id, nc%E_untracked_varid,  self%E_untracked,  start=[tslot]), &
                                   "netcdf_io_get_t0_values_system E_untracked_varid" )

         end if

         deallocate(vals)
         call nc%close()
      end associate
      
      return
   end subroutine swiftest_io_netcdf_get_t0_values_system


   module subroutine swiftest_io_netcdf_initialize_output(self, param)
      !! author: Carlisle A. Wishard, Dana Singh, and David A. Minton
      !!
      !! Initialize a NetCDF file nbody_system and defines all variables.
      use, intrinsic :: ieee_arithmetic
      implicit none
      ! Arguments
      class(swiftest_netcdf_parameters), intent(inout) :: self  !! Parameters used to for writing a NetCDF dataset to file
      class(swiftest_parameters),        intent(in)    :: param !! Current run configuration parameters 
      ! Internals
      integer(I4B) :: nvar, varid, vartype
      real(DP) :: dfill
      real(SP) :: sfill
      integer(I4B), parameter :: NO_FILL = 0
      logical :: fileExists
      character(len=STRMAX) :: errmsg
      integer(I4B) :: ndims

      associate(nc => self)

         dfill = ieee_value(dfill, IEEE_QUIET_NAN)
         sfill = ieee_value(sfill, IEEE_QUIET_NAN)

         select case (param%out_type)
         case("NETCDF_FLOAT")
            nc%out_type = NF90_FLOAT
         case("NETCDF_DOUBLE")
            nc%out_type = NF90_DOUBLE
         case default
            write(*,*) trim(adjustl(param%out_type)), " is an invalid OUT_TYPE"
         end select

         ! Check if the file exists, and if it does, delete it
         inquire(file=nc%file_name, exist=fileExists)
         if (fileExists) then
            open(unit=LUN, file=nc%file_name, status="old", err=667, iomsg=errmsg)
            close(unit=LUN, status="delete")
         end if

         ! Create the file
         call netcdf_io_check( nf90_create(nc%file_name, NF90_NETCDF4, nc%id), "netcdf_io_initialize_output nf90_create" )
         nc%lfile_is_open = .true.

         ! Dimensions
         call netcdf_io_check( nf90_def_dim(nc%id, nc%time_dimname, NF90_UNLIMITED, nc%time_dimid), &
                               "netcdf_io_initialize_output nf90_def_dim time_dimid" ) ! Simulation time dimension
         call netcdf_io_check( nf90_def_dim(nc%id, nc%space_dimname, NDIM, nc%space_dimid), &
                               "netcdf_io_initialize_output nf90_def_dim space_dimid" ) ! 3D space dimension
         call netcdf_io_check( nf90_def_dim(nc%id, nc%name_dimname, NF90_UNLIMITED, nc%name_dimid), &
                               "netcdf_io_initialize_output nf90_def_dim name_dimid" ) ! dimension to store particle id numbers
         call netcdf_io_check( nf90_def_dim(nc%id, nc%str_dimname, NAMELEN, nc%str_dimid), &
                               "netcdf_io_initialize_output nf90_def_dim str_dimid"  )  ! Dimension for string variables 

         ! Dimension coordinates
         call netcdf_io_check( nf90_def_var(nc%id, nc%time_dimname, nc%out_type, nc%time_dimid, nc%time_varid), &
                              "netcdf_io_initialize_output nf90_def_var time_varid"  )
         call netcdf_io_check( nf90_def_var(nc%id, nc%space_dimname, NF90_CHAR, nc%space_dimid, nc%space_varid), &
                               "netcdf_io_initialize_output nf90_def_var space_varid"  )
         call netcdf_io_check( nf90_def_var(nc%id, nc%name_dimname, NF90_CHAR, [nc%str_dimid, nc%name_dimid], nc%name_varid), &
                               "netcdf_io_initialize_output nf90_def_var name_varid"  )

         ! Variables
         call netcdf_io_check( nf90_def_var(nc%id, nc%id_varname, NF90_INT, nc%name_dimid, nc%id_varid), &
                               "netcdf_io_initialize_output nf90_def_var id_varid"  )
         call netcdf_io_check( nf90_def_var(nc%id, nc%status_varname, NF90_INT, [nc%name_dimid, nc%time_dimid], nc%status_varid), &
                               "netcdf_io_initialize_output nf90_def_var status_varid"  )
         call netcdf_io_check( nf90_def_var(nc%id, nc%npl_varname, NF90_INT, nc%time_dimid, nc%npl_varid), &
                               "netcdf_io_initialize_output nf90_def_var npl_varid"  )
         call netcdf_io_check( nf90_def_var(nc%id, nc%ntp_varname, NF90_INT, nc%time_dimid, nc%ntp_varid), &
                               "netcdf_io_initialize_output nf90_def_var ntp_varid"  )
         if (param%lmtiny_pl) call netcdf_io_check( nf90_def_var(nc%id, nc%nplm_varname, NF90_INT, nc%time_dimid, nc%nplm_varid), &
                               "netcdf_io_initialize_output nf90_def_var nplm_varid"  )
         call netcdf_io_check( nf90_def_var(nc%id, nc%ptype_varname, NF90_CHAR, [nc%str_dimid, nc%name_dimid], nc%ptype_varid), &
                               "netcdf_io_initialize_output nf90_def_var ptype_varid"  )

         if ((param%out_form == "XV") .or. (param%out_form == "XVEL")) then
            call netcdf_io_check( nf90_def_var(nc%id, nc%rh_varname,  nc%out_type, [nc%space_dimid, nc%name_dimid, nc%time_dimid], &
                                  nc%rh_varid), "netcdf_io_initialize_output nf90_def_var rh_varid"  )
            call netcdf_io_check( nf90_def_var(nc%id, nc%vh_varname,  nc%out_type, [nc%space_dimid, nc%name_dimid, nc%time_dimid], &
                                 nc%vh_varid), "netcdf_io_initialize_output nf90_def_var vh_varid"  )

            !! When GR is enabled, we need to save the pseudovelocity vectors in addition to the true heliocentric velocity vectors,
            !! otherwise !! we cannnot expect bit-identical runs from restarted runs with GR enabled due to floating point errors 
            !! during the conversion.
            if (param%lgr) then
               call netcdf_io_check( nf90_def_var(nc%id, nc%gr_pseudo_vh_varname, nc%out_type, &
                                     [nc%space_dimid, nc%name_dimid, nc%time_dimid], nc%gr_pseudo_vh_varid), &
                                     "netcdf_io_initialize_output nf90_def_var gr_psuedo_vh_varid"  )
               nc%lpseudo_vel_exists = .true.
            end if

         end if
      
         if ((param%out_form == "EL") .or. (param%out_form == "XVEL")) then
            call netcdf_io_check( nf90_def_var(nc%id, nc%a_varname, nc%out_type, [nc%name_dimid, nc%time_dimid], nc%a_varid), &
                                  "netcdf_io_initialize_output nf90_def_var a_varid"  )
            call netcdf_io_check( nf90_def_var(nc%id, nc%e_varname, nc%out_type, [nc%name_dimid, nc%time_dimid], nc%e_varid), &
                                  "netcdf_io_initialize_output nf90_def_var e_varid"  )
            call netcdf_io_check( nf90_def_var(nc%id, nc%inc_varname, nc%out_type, [nc%name_dimid, nc%time_dimid], nc%inc_varid), &
                                  "netcdf_io_initialize_output nf90_def_var inc_varid"  )
            call netcdf_io_check( nf90_def_var(nc%id, nc%capom_varname, nc%out_type, [nc%name_dimid, nc%time_dimid], &
                                               nc%capom_varid), &
                                  "netcdf_io_initialize_output nf90_def_var capom_varid"  )
            call netcdf_io_check( nf90_def_var(nc%id, nc%omega_varname, nc%out_type, [nc%name_dimid, nc%time_dimid], &
                                               nc%omega_varid), &
                                  "netcdf_io_initialize_output nf90_def_var omega_varid"  )
            call netcdf_io_check( nf90_def_var(nc%id, nc%capm_varname, nc%out_type, [nc%name_dimid, nc%time_dimid], &
                                               nc%capm_varid), &
                                  "netcdf_io_initialize_output nf90_def_var capm_varid"  )
            call netcdf_io_check( nf90_def_var(nc%id, nc%varpi_varname, nc%out_type, [nc%name_dimid, nc%time_dimid], &
                                               nc%varpi_varid), &
                                  "netcdf_io_initialize_output nf90_def_var varpi_varid"  )
            call netcdf_io_check( nf90_def_var(nc%id, nc%lam_varname, nc%out_type, [nc%name_dimid, nc%time_dimid], nc%lam_varid), &
                                  "netcdf_io_initialize_output nf90_def_var lam_varid"  )
            call netcdf_io_check( nf90_def_var(nc%id, nc%f_varname, nc%out_type, [nc%name_dimid, nc%time_dimid], nc%f_varid), &
                                  "netcdf_io_initialize_output nf90_def_var f_varid"  )
            call netcdf_io_check( nf90_def_var(nc%id, nc%cape_varname, nc%out_type, [nc%name_dimid, nc%time_dimid], nc%cape_varid),&
                                  "netcdf_io_initialize_output nf90_def_var cape_varid"  )
         end if

         call netcdf_io_check( nf90_def_var(nc%id, nc%Gmass_varname, nc%out_type, [nc%name_dimid, nc%time_dimid], nc%Gmass_varid), &
                                  "netcdf_io_initialize_output nf90_def_var Gmass_varid"  )
         call netcdf_io_check( nf90_def_var(nc%id, nc%mass_varname, nc%out_type, [nc%name_dimid, nc%time_dimid], nc%mass_varid), &
                                  "netcdf_io_initialize_output nf90_def_var mass_varid"  )
         call netcdf_io_check( nf90_def_var(nc%id, nc%rhill_varname, nc%out_type, [nc%name_dimid, nc%time_dimid], nc%rhill_varid), &
                                  "netcdf_io_initialize_output nf90_def_var rhill_varid"  )

         if (param%lclose) then
            call netcdf_io_check( nf90_def_var(nc%id, nc%radius_varname, nc%out_type, [nc%name_dimid, nc%time_dimid], &
                                               nc%radius_varid), &
                                  "netcdf_io_initialize_output nf90_def_var radius_varid"  )

            call netcdf_io_check( nf90_def_var(nc%id, nc%origin_time_varname, nc%out_type, nc%name_dimid, nc%origin_time_varid), &
                                  "netcdf_io_initialize_output nf90_def_var origin_time_varid"  )
            call netcdf_io_check( nf90_def_var(nc%id, nc%origin_type_varname, NF90_CHAR, [nc%str_dimid, nc%name_dimid], &
                                    nc%origin_type_varid), &
                                  "netcdf_io_initialize_output nf90_create"  )
            call netcdf_io_check( nf90_def_var(nc%id, nc%origin_rh_varname, nc%out_type, [nc%space_dimid, nc%name_dimid], &
                                               nc%origin_rh_varid), &
                                  "netcdf_io_initialize_output nf90_def_var origin_rh_varid"  )
            call netcdf_io_check( nf90_def_var(nc%id, nc%origin_vh_varname, nc%out_type, [nc%space_dimid, nc%name_dimid], &
                                               nc%origin_vh_varid), &
                                  "netcdf_io_initialize_output nf90_def_var origin_vh_varid"  )

            call netcdf_io_check( nf90_def_var(nc%id, nc%collision_id_varname, NF90_INT, nc%name_dimid, nc%collision_id_varid), &
                                  "netcdf_io_initialize_output nf90_def_var collision_id_varid"  )
            call netcdf_io_check( nf90_def_var(nc%id, nc%discard_time_varname, nc%out_type, nc%name_dimid, nc%discard_time_varid), &
                                  "netcdf_io_initialize_output nf90_def_var discard_time_varid"  )
            call netcdf_io_check( nf90_def_var(nc%id, nc%discard_rh_varname, nc%out_type, [nc%space_dimid, nc%name_dimid], &
                                               nc%discard_rh_varid), &
                                  "netcdf_io_initialize_output nf90_def_var discard_rh_varid"  )
            call netcdf_io_check( nf90_def_var(nc%id, nc%discard_vh_varname, nc%out_type, [nc%space_dimid, nc%name_dimid], &
                                               nc%discard_vh_varid), &
                                  "netcdf_io_initialize_output nf90_def_var discard_vh_varid"  )
            call netcdf_io_check( nf90_def_var(nc%id, nc%discard_body_id_varname, NF90_INT, nc%name_dimid, &
                                               nc%discard_body_id_varid), &
                                  "netcdf_io_initialize_output nf90_def_var discard_body_id_varid"  )
         end if

         if (param%lrotation) then
            call netcdf_io_check( nf90_def_var(nc%id, nc%Ip_varname, nc%out_type, [nc%space_dimid, nc%name_dimid, nc%time_dimid], &
                                               nc%Ip_varid), &
                                  "netcdf_io_initialize_output nf90_def_var Ip_varid"  )
            call netcdf_io_check( nf90_def_var(nc%id, nc%rot_varname, nc%out_type, [nc%space_dimid, nc%name_dimid, nc%time_dimid], &
                                               nc%rot_varid), &
                                  "netcdf_io_initialize_output nf90_def_var rot_varid"  )
         end if

         ! if (param%ltides) then
         !    call netcdf_io_check( nf90_def_var(nc%id, nc%k2_varname, nc%out_type, [nc%name_dimid, nc%time_dimid], nc%k2_varid), &
         !                        "netcdf_io_initialize_output nf90_def_var k2_varid"  )
         !    call netcdf_io_check( nf90_def_var(nc%id, nc%q_varname, nc%out_type, [nc%name_dimid, nc%time_dimid], nc%Q_varid), &
         !                        "netcdf_io_initialize_output nf90_def_var Q_varid"  )
         ! end if

         if (param%lenergy) then
            call netcdf_io_check( nf90_def_var(nc%id, nc%ke_orb_varname, nc%out_type, nc%time_dimid, nc%KE_orb_varid), &
                                  "netcdf_io_initialize_output nf90_def_var KE_orb_varid"  )
            call netcdf_io_check( nf90_def_var(nc%id, nc%ke_spin_varname, nc%out_type, nc%time_dimid, nc%KE_spin_varid), &
                                  "netcdf_io_initialize_output nf90_def_var KE_spin_varid"  )
            call netcdf_io_check( nf90_def_var(nc%id, nc%pe_varname, nc%out_type, nc%time_dimid, nc%PE_varid), &
                                  "netcdf_io_initialize_output nf90_def_var PE_varid"  )
            call netcdf_io_check( nf90_def_var(nc%id, nc%be_varname, nc%out_type, nc%time_dimid, nc%BE_varid), &
                                  "netcdf_io_initialize_output nf90_def_var BE_varid"  )
            call netcdf_io_check( nf90_def_var(nc%id, nc%te_varname, nc%out_type, nc%time_dimid, nc%TE_varid), &
                                  "netcdf_io_initialize_output nf90_def_var TE_varid"  )
            call netcdf_io_check( nf90_def_var(nc%id, nc%L_orbit_varname, nc%out_type, [nc%space_dimid, nc%time_dimid], &
                                               nc%L_orbit_varid), &
                                  "netcdf_io_initialize_output nf90_def_var L_orbit_varid"  )
            call netcdf_io_check( nf90_def_var(nc%id, nc%L_spin_varname, nc%out_type, [nc%space_dimid, nc%time_dimid], &
                                               nc%L_spin_varid), &
                                  "netcdf_io_initialize_output nf90_def_var L_spin_varid"  )
            call netcdf_io_check( nf90_def_var(nc%id, nc%L_escape_varname, nc%out_type, [nc%space_dimid, nc%time_dimid], &
                                               nc%L_escape_varid), &
                                  "netcdf_io_initialize_output nf90_def_var L_escape_varid"  )
            call netcdf_io_check( nf90_def_var(nc%id, nc%E_collisions_varname, nc%out_type, nc%time_dimid, nc%E_collisions_varid), &
                                  "netcdf_io_initialize_output nf90_def_var E_collisions_varid"  )
            call netcdf_io_check( nf90_def_var(nc%id, nc%E_untracked_varname, nc%out_type, nc%time_dimid, nc%E_untracked_varid), &
                                  "netcdf_io_initialize_output nf90_def_var E_untracked_varid"  )
            call netcdf_io_check( nf90_def_var(nc%id, nc%GMescape_varname, nc%out_type, nc%time_dimid, nc%GMescape_varid), &
                                  "netcdf_io_initialize_output nf90_def_var GMescape_varid"  )
         end if

         call netcdf_io_check( nf90_def_var(nc%id, nc%j2rp2_varname, nc%out_type, nc%time_dimid, nc%j2rp2_varid), &
                                  "netcdf_io_initialize_output nf90_def_var j2rp2_varid"  )
         call netcdf_io_check( nf90_def_var(nc%id, nc%j4rp4_varname, nc%out_type, nc%time_dimid, nc%j4rp4_varid), &
                                  "netcdf_io_initialize_output nf90_def_var j4rp4_varid"  )

         ! Set fill mode to NaN for all variables
         call netcdf_io_check( nf90_inquire(nc%id, nVariables=nvar), "netcdf_io_initialize_output nf90_inquire nVariables" )
         do varid = 1, nvar
            call netcdf_io_check( nf90_inquire_variable(nc%id, varid, xtype=vartype, ndims=ndims), &
                                  "netcdf_io_initialize_output nf90_inquire_variable"  )
            select case(vartype)
            case(NF90_INT)
               if (varid == nc%status_varid) then
                  ! Be sure the status variable fill value is the INACTIVE symbolic value
                  call netcdf_io_check( nf90_def_var_fill(nc%id, varid, NO_FILL, INACTIVE), &
                                  "netcdf_io_netcdf_initialize_output nf90_def_var_fill status variable"  )
               else
                  call netcdf_io_check( nf90_def_var_fill(nc%id, varid, NO_FILL, NF90_FILL_INT), &
                                  "netcdf_io_initialize_output nf90_def_var_fill NF90_INT"  )
               end if
            case(NF90_FLOAT)
               call netcdf_io_check( nf90_def_var_fill(nc%id, varid, NO_FILL, sfill), &
                                  "netcdf_io_initialize_output nf90_def_var_fill NF90_FLOAT"  )
            case(NF90_DOUBLE)
               call netcdf_io_check( nf90_def_var_fill(nc%id, varid, NO_FILL, dfill), &
                                  "netcdf_io_initialize_output nf90_def_var_fill NF90_DOUBLE"  )
            case(NF90_CHAR)
               call netcdf_io_check( nf90_def_var_fill(nc%id, varid, NO_FILL, 0), &
                                  "netcdf_io_initialize_output nf90_def_var_fill NF90_CHAR"  )
            end select
         end do

         ! Set special fill mode for discard time so that we can make use of it for non-discarded bodies.
         if (param%lclose) then
            select case (vartype)
            case(NF90_FLOAT)
               call netcdf_io_check( nf90_def_var_fill(nc%id, nc%discard_time_varid, NO_FILL, huge(1.0_SP)), &
                                  "netcdf_io_initialize_output nf90_def_var_fill discard_time NF90_FLOAT"  )
            case(NF90_DOUBLE)
               call netcdf_io_check( nf90_def_var_fill(nc%id, nc%discard_time_varid, NO_FILL, huge(1.0_DP)), &
                                  "netcdf_io_initialize_output nf90_def_var_fill discard_time NF90_DOUBLE"  )
            end select
         end if

         ! Take the file out of define mode
         call netcdf_io_check( nf90_enddef(nc%id), "netcdf_io_initialize_output nf90_enddef"  )

         ! Add in the space dimension coordinates
         call netcdf_io_check( nf90_put_var(nc%id, nc%space_varid, nc%space_coords, start=[1], count=[NDIM]), &
                                  "netcdf_io_initialize_output nf90_put_var space"  )

      end associate
      return

      667 continue
      write(*,*) "Error creating NetCDF output file. " // trim(adjustl(errmsg))
      call base_util_exit(FAILURE,param%display_unit)
   end subroutine swiftest_io_netcdf_initialize_output


   module subroutine swiftest_io_netcdf_open(self, param, readonly)
      !! author: Carlisle A. Wishard, Dana Singh, and David A. Minton
      !!
      !! Opens a NetCDF file and does the variable inquiries to activate variable ids
      implicit none
      ! Arguments
      class(swiftest_netcdf_parameters), intent(inout) :: self     !! Parameters used to identify a particular NetCDF dataset
      class(swiftest_parameters),        intent(inout) :: param    !! Current run configuration parameters
      logical, optional,                 intent(in)    :: readonly !! Logical flag indicating that this should be open read only
      ! Internals
      integer(I4B) :: mode, status
      character(len=STRMAX) :: errmsg

      mode = NF90_WRITE
      if (present(readonly)) then
         if (readonly) mode = NF90_NOWRITE
      end if

      associate(nc => self)
         write(errmsg,*) "swiftest_io_netcdf_open nf90_open ",trim(adjustl(nc%file_name))
         call netcdf_io_check( nf90_open(nc%file_name, mode, nc%id), errmsg)
         self%lfile_is_open = .true.

         ! Dimensions
         call netcdf_io_check( nf90_inq_dimid(nc%id, nc%time_dimname, nc%time_dimid), &
                                  "swiftest_io_netcdf_open nf90_inq_dimid time_dimid"  )
         call netcdf_io_check( nf90_inquire_dimension(nc%id, nc%time_dimid, nc%time_dimname, len=nc%max_tslot), &
                                  "swiftest_io_netcdf_open nf90_inquire_dimension max_tslot"  )
         call netcdf_io_check( nf90_inq_dimid(nc%id, nc%space_dimname, nc%space_dimid), &
                                  "swiftest_io_netcdf_open nf90_inq_dimid space_dimid"  )
         call netcdf_io_check( nf90_inq_dimid(nc%id, nc%name_dimname, nc%name_dimid), &
                                  "swiftest_io_netcdf_open nf90_inq_dimid name_dimid"  )
         call netcdf_io_check( nf90_inq_dimid(nc%id, nc%str_dimname, nc%str_dimid), &
                                  "swiftest_io_netcdf_open nf90_inq_dimid str_dimid"  )

         ! Dimension coordinates
         call netcdf_io_check( nf90_inq_varid(nc%id, nc%time_dimname, nc%time_varid), &
                                  "swiftest_io_netcdf_open nf90_inq_varid time_varid" )
         call netcdf_io_check( nf90_inq_varid(nc%id, nc%space_dimname, nc%space_varid), &
                                  "swiftest_io_netcdf_open nf90_inq_varid space_varid" )
         call netcdf_io_check( nf90_inq_varid(nc%id, nc%name_dimname, nc%name_varid), &
                                  "swiftest_io_netcdf_open nf90_inq_varid name_varid" )

         ! Required Variables
         call netcdf_io_check( nf90_inq_varid(nc%id, nc%id_varname, nc%id_varid), &
                                  "swiftest_io_netcdf_open nf90_inq_varid name_varid" )
         call netcdf_io_check( nf90_inq_varid(nc%id, nc%Gmass_varname, nc%Gmass_varid), &
                                  "swiftest_io_netcdf_open nf90_inq_varid Gmass_varid" )

         if ((param%out_form == "XV") .or. (param%out_form == "XVEL")) then
            call netcdf_io_check( nf90_inq_varid(nc%id, nc%rh_varname, nc%rh_varid), &
                                  "swiftest_io_netcdf_open nf90_inq_varid rh_varid" )
            call netcdf_io_check( nf90_inq_varid(nc%id, nc%vh_varname, nc%vh_varid), &
                                  "swiftest_io_netcdf_open nf90_inq_varid vh_varid" )

            if (param%lgr) then
               !! check if pseudovelocity vectors exist in this file. If they are, set the correct flag so we know whe should not 
               !! do the conversion.
               status = nf90_inq_varid(nc%id, nc%gr_pseudo_vh_varname, nc%gr_pseudo_vh_varid)
               nc%lpseudo_vel_exists = (status == NF90_NOERR)
               if (param%lrestart .and. .not.nc%lpseudo_vel_exists) then
                  write(*,*) "Warning! Pseudovelocity not found in input file for GR enabled run. " &
                          // "If this is a restarted run, bit-identical trajectories are not guarunteed!"
               end if

            end if
         end if

         if ((param%out_form == "EL") .or. (param%out_form == "XVEL")) then
            call netcdf_io_check( nf90_inq_varid(nc%id, nc%a_varname, nc%a_varid), &
                                  "swiftest_io_netcdf_open nf90_inq_varid a_varid" )
            call netcdf_io_check( nf90_inq_varid(nc%id, nc%e_varname, nc%e_varid), &
                                  "swiftest_io_netcdf_open nf90_inq_varid e_varid" )
            call netcdf_io_check( nf90_inq_varid(nc%id, nc%inc_varname, nc%inc_varid), &
                                  "swiftest_io_netcdf_open nf90_inq_varid inc_varid" )
            call netcdf_io_check( nf90_inq_varid(nc%id, nc%capom_varname, nc%capom_varid), &
                                  "swiftest_io_netcdf_open nf90_inq_varid capom_varid" )
            call netcdf_io_check( nf90_inq_varid(nc%id, nc%omega_varname, nc%omega_varid), &
                                  "swiftest_io_netcdf_open nf90_inq_varid omega_varid" )
            call netcdf_io_check( nf90_inq_varid(nc%id, nc%capm_varname, nc%capm_varid), &
                                  "swiftest_io_netcdf_open nf90_inq_varid capm_varid" )
         end if

         if (param%lclose) then
            call netcdf_io_check( nf90_inq_varid(nc%id, nc%radius_varname, nc%radius_varid), &
                                  "swiftest_io_netcdf_open nf90_inq_varid radius_varid" )
         end if 

         if (param%lrotation) then
            call netcdf_io_check( nf90_inq_varid(nc%id, nc%Ip_varname, nc%Ip_varid), &
                                  "swiftest_io_netcdf_open nf90_inq_varid Ip_varid" )
            call netcdf_io_check( nf90_inq_varid(nc%id, nc%rot_varname, nc%rot_varid), &
                                  "swiftest_io_netcdf_open nf90_inq_varid rot_varid" )
         end if

         ! if (param%ltides) then
         !    call netcdf_io_check( nf90_inq_varid(nc%id, nc%k2_varname, nc%k2_varid), &
         !                         "swiftest_io_netcdf_open nf90_inq_varid k2_varid" )
         !    call netcdf_io_check( nf90_inq_varid(nc%id, nc%q_varname, nc%Q_varid), &
         !                         "swiftest_io_netcdf_open nf90_inq_varid Q_varid" )
         ! end if

         ! Optional variables The User Doesn't Need to Know About
         status = nf90_inq_varid(nc%id, nc%mass_varname, nc%mass_varid)
         status = nf90_inq_varid(nc%id, nc%rhill_varname, nc%rhill_varid)
         status = nf90_inq_varid(nc%id, nc%npl_varname, nc%npl_varid)
         status = nf90_inq_varid(nc%id, nc%status_varname, nc%status_varid)
         status = nf90_inq_varid(nc%id, nc%ntp_varname, nc%ntp_varid)
         status = nf90_inq_varid(nc%id, nc%j2rp2_varname, nc%j2rp2_varid)
         status = nf90_inq_varid(nc%id, nc%j4rp4_varname, nc%j4rp4_varid)
         status = nf90_inq_varid(nc%id, nc%ptype_varname, nc%ptype_varid)
         status = nf90_inq_varid(nc%id, nc%varpi_varname, nc%varpi_varid)
         status = nf90_inq_varid(nc%id, nc%lam_varname, nc%lam_varid)
         status = nf90_inq_varid(nc%id, nc%f_varname, nc%f_varid)
         status = nf90_inq_varid(nc%id, nc%cape_varname, nc%cape_varid)

         if (param%lmtiny_pl) status = nf90_inq_varid(nc%id, nc%nplm_varname, nc%nplm_varid)

         if (param%lclose) then
            status = nf90_inq_varid(nc%id, nc%origin_type_varname, nc%origin_type_varid)
            status = nf90_inq_varid(nc%id, nc%origin_time_varname, nc%origin_time_varid)
            status = nf90_inq_varid(nc%id, nc%origin_rh_varname, nc%origin_rh_varid)
            status = nf90_inq_varid(nc%id, nc%origin_vh_varname, nc%origin_vh_varid)
            status = nf90_inq_varid(nc%id, nc%collision_id_varname, nc%collision_id_varid)
            status = nf90_inq_varid(nc%id, nc%discard_time_varname, nc%discard_time_varid)
            status = nf90_inq_varid(nc%id, nc%discard_rh_varname, nc%discard_rh_varid)
            status = nf90_inq_varid(nc%id, nc%discard_vh_varname, nc%discard_vh_varid)
            status = nf90_inq_varid(nc%id, nc%discard_body_id_varname, nc%discard_body_id_varid)
         end if

         if (param%lenergy) then
            status = nf90_inq_varid(nc%id, nc%ke_orb_varname, nc%KE_orb_varid)
            status = nf90_inq_varid(nc%id, nc%ke_spin_varname, nc%KE_spin_varid)
            status = nf90_inq_varid(nc%id, nc%pe_varname, nc%PE_varid)
            status = nf90_inq_varid(nc%id, nc%be_varname, nc%BE_varid)
            status = nf90_inq_varid(nc%id, nc%te_varname, nc%TE_varid)
            status = nf90_inq_varid(nc%id, nc%L_orbit_varname, nc%L_orbit_varid)
            status = nf90_inq_varid(nc%id, nc%L_spin_varname, nc%L_spin_varid)
            status = nf90_inq_varid(nc%id, nc%L_escape_varname, nc%L_escape_varid)
            status = nf90_inq_varid(nc%id, nc%E_collisions_varname, nc%E_collisions_varid)
            status = nf90_inq_varid(nc%id, nc%E_untracked_varname, nc%E_untracked_varid)
            status = nf90_inq_varid(nc%id, nc%GMescape_varname, nc%GMescape_varid)
         end if

      end associate

      return
   end subroutine swiftest_io_netcdf_open


   module subroutine swiftest_io_netcdf_get_valid_masks(self, plmask, tpmask, plmmask, Gmtiny)
      !! author: David A. Minton
      !!
      !! Given an open NetCDF, returns logical masks indicating which bodies in the body arrays are active pl and tp type at the 
      !! current time. Uses the value of tslot stored in the NetCDF parameter object as the definition of current time
      use, intrinsic :: ieee_exceptions
      use, intrinsic :: ieee_arithmetic
      implicit none
      ! Arguments
      class(swiftest_netcdf_parameters),  intent(inout)          :: self    !! Parameters used to identify a particular NetCDF 
                                                                            !!   dataset
      logical, dimension(:), allocatable, intent(out)            :: plmask  !! Logical mask indicating which bodies are massive 
                                                                            !!   bodies
      logical, dimension(:), allocatable, intent(out)            :: tpmask  !! Logical mask indicating which bodies are test 
                                                                            !!   particles
      logical, dimension(:), allocatable, intent(out), optional  :: plmmask !! Logical mask indicating which bodies are fully 
                                                                            !!   interacting massive bodies
      real(DP),                           intent(in),  optional  :: Gmtiny  !! The cutoff G*mass between semi-interacting and fully
                                                                            !!   interacting massive bodies
      ! Internals
      real(DP), dimension(:), allocatable :: Gmass, a
      real(DP), dimension(:,:), allocatable :: rh
      integer(I4B), dimension(:), allocatable :: body_status 
      logical, dimension(:), allocatable :: lvalid
      integer(I4B) :: idmax, status
      logical, dimension(size(IEEE_ALL))      :: fpe_halting_modes

      call ieee_get_halting_mode(IEEE_ALL,fpe_halting_modes)  ! Save the current halting modes so we can turn them off temporarily
      call ieee_set_halting_mode(IEEE_ALL,.false.)

      call netcdf_io_check( nf90_inquire_dimension(self%id, self%name_dimid, len=idmax), &
                                  "swiftest_io_netcdf_get_valid_masks nf90_inquire_dimension name_dimid"  )

      allocate(Gmass(idmax))
      allocate(tpmask(idmax))
      allocate(plmask(idmax))
      allocate(lvalid(idmax))
      associate(tslot => self%tslot)

         call netcdf_io_check( nf90_get_var(self%id, self%Gmass_varid, Gmass, start=[1,tslot], count=[idmax,1]), &
                                  "swiftest_io_netcdf_get_valid_masks nf90_getvar Gmass_varid"  )

         status = nf90_inq_varid(self%id, self%status_varname, self%status_varid) 
         if (status == NF90_NOERR) then
            allocate(body_status(idmax))
            call netcdf_io_check( nf90_get_var(self%id, self%status_varid, body_status, start=[1, tslot], count=[idmax,1]), &
                                  "swiftest_io_netcdf_get_valid_masks nf90_getvar status_varid"  )
            lvalid(:) = body_status(:) /= INACTIVE
         else
            status = nf90_inq_varid(self%id, self%rh_varname, self%rh_varid) 
            if (status == NF90_NOERR) then
               allocate(rh(NDIM,idmax))
               call netcdf_io_check( nf90_get_var(self%id, self%rh_varid, rh, start=[1, 1, tslot], count=[NDIM,idmax,1]), &
                                  "swiftest_io_netcdf_get_valid_masks nf90_getvar rh_varid"  )
               lvalid(:) = ieee_is_normal(rh(1,:)) 
            else
               status = nf90_inq_varid(self%id, self%a_varname, self%a_varid) 
               if (status == NF90_NOERR) then
                  allocate(a(idmax))
                  call netcdf_io_check( nf90_get_var(self%id, self%a_varid, a, start=[1, tslot], count=[idmax,1]), &
                                  "swiftest_io_netcdf_get_valid_masks nf90_getvar a_varid"  )
                  lvalid(:) = ieee_is_normal(a(:))
               else
                  lvalid(:) = .false.
               end if
            end if
         end if

         plmask(:) = ieee_is_normal(Gmass(:))
         where(plmask(:)) plmask(:) = Gmass(:) > 0.0_DP
         tpmask(:) = .not. plmask(:)
         plmask(1) = .false. ! This is the central body

         ! Select only active bodies
         plmask(:) = plmask(:) .and. lvalid(:)
         tpmask(:) = tpmask(:) .and. lvalid(:)

         if (present(plmmask) .and. present(Gmtiny)) then
            allocate(plmmask, source=plmask)
            where(plmask(:))
               plmmask = Gmass(:) > Gmtiny
            endwhere
         end if

         call ieee_set_halting_mode(IEEE_ALL,fpe_halting_modes)

      end associate

      return
   end subroutine swiftest_io_netcdf_get_valid_masks


   module function swiftest_io_netcdf_read_frame_system(self, nc, param) result(ierr)
      !! author: The Purdue Swiftest Team - David A. Minton, Carlisle A. Wishard, Jennifer L.L. Pouplin, and Jacob R. Elliott
      !!
      !! Read a frame (header plus records for each massive body and active test particle) from an output binary file
      implicit none
      ! Arguments
      class(swiftest_nbody_system),      intent(inout) :: self  !! Swiftest nbody_system object
      class(swiftest_netcdf_parameters), intent(inout) :: nc    !! Parameters used to identify a particular NetCDF dataset
      class(swiftest_parameters),        intent(inout) :: param !! Current run configuration parameters 
      ! Return
      integer(I4B)                                :: ierr  !! Error code: returns 0 if the read is successful
      ! Internals
      integer(I4B)                              :: i, idmax, npl_check, ntp_check, str_max, status, npl, ntp
      real(DP), dimension(:), allocatable       :: rtemp
      real(DP), dimension(:,:), allocatable     :: vectemp
      integer(I4B), dimension(:), allocatable   :: itemp
      logical, dimension(:), allocatable        :: tpmask, plmask

      associate(cb => self%cb, pl => self%pl, tp => self%tp, tslot => nc%tslot)
         call nc%open(param, readonly=.true.)
         call nc%find_tslot(self%t, tslot)
         call netcdf_io_check( nf90_inquire_dimension(nc%id, nc%time_dimid, len=nc%max_tslot), &
                                  "netcdf_io_read_frame_system nf90_inquire_dimension time_dimid"  )
         tslot = min(tslot, nc%max_tslot)

         call self%read_hdr(nc, param)

         ! Save these values as variables as they get reset by the setup method
         npl = pl%nbody
         ntp = tp%nbody
         call pl%setup(npl, param)
         call tp%setup(ntp, param)
         call netcdf_io_check( nf90_inquire_dimension(nc%id, nc%name_dimid, len=idmax), &
                                  "netcdf_io_read_frame_system nf90_inquire_dimension name_dimid"  )
         allocate(rtemp(idmax))
         allocate(vectemp(NDIM,idmax))
         allocate(itemp(idmax))

         call netcdf_io_check( nf90_inquire_dimension(nc%id, nc%str_dimid, len=str_max), &
                                  "netcdf_io_read_frame_system nf90_inquire_dimension str_dimid"  )

         call nc%get_valid_masks(plmask, tpmask)

         ! Check to make sure the number of bodies is correct
         npl_check = count(plmask(:))
         ntp_check = count(tpmask(:))

         if (npl_check /= npl) then
            write(*,*) "Error reading in NetCDF file: The recorded value of npl does not match the number of active massive bodies"
            write(*,*) "Recorded: ",npl
            write(*,*) "Active  : ",npl_check
         end if

         if (ntp_check /= ntp) then
            write(*,*) "Error reading in NetCDF file: The recorded value of ntp does not match the number of active test particles"
            write(*,*) "Recorded: ",ntp
            write(*,*) "Active  : ",ntp_check
            call base_util_exit(FAILURE,param%display_unit)
         end if

         ! Now read in each variable and split the outputs by body type
         if ((param%in_form == "XV") .or. (param%in_form == "XVEL")) then
            call netcdf_io_check( nf90_get_var(nc%id, nc%rh_varid, vectemp, start=[1, 1, tslot], count=[NDIM,idmax,1]), &
                                  "netcdf_io_read_frame_system nf90_getvar rh_varid"  )
            do i = 1, NDIM
               if (npl > 0) pl%rh(i,:) = pack(vectemp(i,:), plmask(:))
               if (ntp > 0) tp%rh(i,:) = pack(vectemp(i,:), tpmask(:))
            end do

            if (param%lgr .and. nc%lpseudo_vel_exists) then
               call netcdf_io_check( nf90_get_var(nc%id, nc%gr_pseudo_vh_varid, vectemp, start=[1, 1, tslot], &
                                                  count=[NDIM,idmax,1]), &
                                  "netcdf_io_read_frame_system nf90_getvar gr_pseudo_vh_varid"  )
               do i = 1, NDIM
                  if (npl > 0) pl%vh(i,:) = pack(vectemp(i,:), plmask(:))
                  if (ntp > 0) tp%vh(i,:) = pack(vectemp(i,:), tpmask(:))
               end do
            else
               call netcdf_io_check( nf90_get_var(nc%id, nc%vh_varid, vectemp, start=[1, 1, tslot], count=[NDIM,idmax,1]), &
                                  "netcdf_io_read_frame_system nf90_getvar vh_varid"  )
               do i = 1, NDIM
                  if (npl > 0) pl%vh(i,:) = pack(vectemp(i,:), plmask(:))
                  if (ntp > 0) tp%vh(i,:) = pack(vectemp(i,:), tpmask(:))
               end do
            end if
         end if

         if ((param%in_form == "EL")  .or. (param%in_form == "XVEL")) then
            call netcdf_io_check( nf90_get_var(nc%id, nc%a_varid, rtemp, start=[1, tslot], count=[idmax,1]), &
                                  "netcdf_io_read_frame_system nf90_getvar a_varid"  )
            if (.not.allocated(pl%a)) allocate(pl%a(npl))
            if (.not.allocated(tp%a)) allocate(tp%a(ntp))
            if (npl > 0) pl%a(:) = pack(rtemp, plmask)
            if (ntp > 0) tp%a(:) = pack(rtemp, tpmask)

            call netcdf_io_check( nf90_get_var(nc%id, nc%e_varid, rtemp, start=[1, tslot], count=[idmax,1]), &
                                  "netcdf_io_read_frame_system nf90_getvar e_varid"  )
            if (.not.allocated(pl%e)) allocate(pl%e(npl))
            if (.not.allocated(tp%e)) allocate(tp%e(ntp))
            if (npl > 0) pl%e(:) = pack(rtemp, plmask)
            if (ntp > 0) tp%e(:) = pack(rtemp, tpmask)

            call netcdf_io_check( nf90_get_var(nc%id, nc%inc_varid, rtemp, start=[1, tslot], count=[idmax,1]), &
                                  "netcdf_io_read_frame_system nf90_getvar inc_varid"  )
            rtemp = rtemp * DEG2RAD
            if (.not.allocated(pl%inc)) allocate(pl%inc(npl))
            if (.not.allocated(tp%inc)) allocate(tp%inc(ntp))
            if (npl > 0) pl%inc(:) = pack(rtemp, plmask)
            if (ntp > 0) tp%inc(:) = pack(rtemp, tpmask)

            call netcdf_io_check( nf90_get_var(nc%id, nc%capom_varid, rtemp, start=[1, tslot], count=[idmax,1]), &
                                  "netcdf_io_read_frame_system nf90_getvar capom_varid"  )
            rtemp = rtemp * DEG2RAD
            if (.not.allocated(pl%capom)) allocate(pl%capom(npl))
            if (.not.allocated(tp%capom)) allocate(tp%capom(ntp))
            if (npl > 0) pl%capom(:) = pack(rtemp, plmask)
            if (ntp > 0) tp%capom(:) = pack(rtemp, tpmask)

            call netcdf_io_check( nf90_get_var(nc%id, nc%omega_varid, rtemp, start=[1, tslot], count=[idmax,1]), &
                                  "netcdf_io_read_frame_system nf90_getvar omega_varid"  )
            rtemp = rtemp * DEG2RAD
            if (.not.allocated(pl%omega)) allocate(pl%omega(npl))
            if (.not.allocated(tp%omega)) allocate(tp%omega(ntp))
            if (npl > 0) pl%omega(:) = pack(rtemp, plmask)
            if (ntp > 0) tp%omega(:) = pack(rtemp, tpmask)

            call netcdf_io_check( nf90_get_var(nc%id, nc%capm_varid, rtemp, start=[1, tslot], count=[idmax,1]), &
                                  "netcdf_io_read_frame_system nf90_getvar capm_varid"  )
            rtemp = rtemp * DEG2RAD
            if (.not.allocated(pl%capm)) allocate(pl%capm(npl))
            if (.not.allocated(tp%capm)) allocate(tp%capm(ntp))
            if (npl > 0) pl%capm(:) = pack(rtemp, plmask)
            if (ntp > 0) tp%capm(:) = pack(rtemp, tpmask)

         end if
      
         call netcdf_io_check( nf90_get_var(nc%id, nc%Gmass_varid, rtemp, start=[1, tslot], count=[idmax,1]), &
                                  "netcdf_io_read_frame_system nf90_getvar Gmass_varid"  )
         cb%Gmass = rtemp(1)
         cb%mass = cb%Gmass / param%GU

         ! Set initial central body mass for Helio bookkeeping
         cb%GM0 = cb%Gmass

         if (npl > 0) then
            pl%Gmass(:) = pack(rtemp, plmask)
            pl%mass(:) = pl%Gmass(:) / param%GU
            if (param%lmtiny_pl) pl%nplm = count(pack(rtemp,plmask) > param%GMTINY )

            status = nf90_get_var(nc%id, nc%rhill_varid, rtemp, start=[1, tslot], count=[idmax,1])
            if (status == NF90_NOERR) then
               pl%rhill(:) = pack(rtemp, plmask)
            end if
         end if

         if (param%lclose) then
            call netcdf_io_check( nf90_get_var(nc%id, nc%radius_varid, rtemp, start=[1, tslot], count=[idmax,1]), &
                                  "netcdf_io_read_frame_system nf90_getvar radius_varid"  )
            cb%radius = rtemp(1)

            if (npl > 0) pl%radius(:) = pack(rtemp, plmask)
         else
            cb%radius = param%rmin
            if (npl > 0) pl%radius(:) = 0.0_DP
         end if
         cb%R0 = cb%radius

         if (param%lrotation) then
            call netcdf_io_check( nf90_get_var(nc%id, nc%Ip_varid,  vectemp, start=[1, 1, tslot], count=[NDIM,idmax,1]), &
                                  "netcdf_io_read_frame_system nf90_getvar Ip_varid"  )
            cb%Ip(:) = vectemp(:,1)
            do i = 1, NDIM
               if (npl > 0) pl%Ip(i,:) = pack(vectemp(i,:), plmask(:))
            end do

            call netcdf_io_check( nf90_get_var(nc%id, nc%rot_varid, vectemp, start=[1, 1, tslot], count=[NDIM,idmax,1]), &
                                  "netcdf_io_read_frame_system nf90_getvar rot_varid"  )
            vectemp(:,:) = vectemp(:,:) * DEG2RAD
            cb%rot(:) = vectemp(:,1) 
            do i = 1, NDIM
               if (npl > 0) pl%rot(i,:) = pack(vectemp(i,:), plmask(:))
            end do

            ! Set initial central body angular momentum for bookkeeping
            cb%L0(:) = cb%Ip(3) * cb%mass * cb%R0**2 * cb%rot(:)         
         end if

         ! if (param%ltides) then
         !    call netcdf_io_check( nf90_get_var(nc%id, nc%k2_varid, rtemp, start=[1, tslot]), &
         !                        "netcdf_io_read_frame_system nf90_getvar k2_varid"  )
         !    cb%k2 = rtemp(1)
         !    if (npl > 0) pl%k2(:) = pack(rtemp, plmask)

         !    call netcdf_io_check( nf90_get_var(nc%id, nc%Q_varid,  rtemp,  start=[1, tslot]), &
         !                        "netcdf_io_read_frame_system nf90_getvar Q_varid"  )
         !    cb%Q = rtemp(1)
         !    if (npl > 0) pl%Q(:) = pack(rtemp, plmask)
         ! end if

         status = nf90_inq_varid(nc%id, nc%j2rp2_varname, nc%j2rp2_varid)
         if (status == NF90_NOERR) then
            call netcdf_io_check( nf90_get_var(nc%id, nc%j2rp2_varid, cb%j2rp2, start=[tslot]), &
                                  "netcdf_io_read_frame_system nf90_getvar j2rp2_varid"  )
         else 
            cb%j2rp2 = 0.0_DP
         end if

         status = nf90_inq_varid(nc%id, nc%j4rp4_varname, nc%j4rp4_varid)   
         if (status == NF90_NOERR) then      
            call netcdf_io_check( nf90_get_var(nc%id, nc%j4rp4_varid, cb%j4rp4, start=[tslot]), &
                                  "netcdf_io_read_frame_system nf90_getvar j4rp4_varid"  )
         else 
            cb%j4rp4 = 0.0_DP
         end if

         call self%read_particle_info(nc, param, plmask, tpmask) 

         if (param%in_form == "EL") then
            call pl%el2xv(cb)
            call tp%el2xv(cb)
         end if
         ! if this is a GR-enabled run, check to see if we got the pseudovelocities in. Otherwise, we'll need to generate them.
         if (param%lgr .and. .not.(nc%lpseudo_vel_exists)) then
            call pl%set_mu(cb)
            call tp%set_mu(cb)
            call pl%v2pv(param)
            call tp%v2pv(param)
         end if
         
      end associate

      call nc%close()

      ierr = 0
      return
   end function swiftest_io_netcdf_read_frame_system


   module subroutine swiftest_io_netcdf_read_hdr_system(self, nc, param) 
      !! author: David A. Minton
      !!
      !! Reads header information (variables that change with time, but not particle id). 
      !! This subroutine swiftest_significantly improves the output over the original binary file, allowing us to track energy, 
      !! momentum, and other quantities that previously were handled as separate output files.
      implicit none
      ! Arguments
      class(swiftest_nbody_system),      intent(inout) :: self  !! Swiftest nbody system object
      class(swiftest_netcdf_parameters), intent(inout) :: nc    !! Parameters used to for reading a NetCDF dataset to file
      class(swiftest_parameters),        intent(inout) :: param !! Current run configuration parameters
      ! Internals
      integer(I4B) :: status
      logical, dimension(:), allocatable        :: plmask, tpmask, plmmask

      associate(tslot => nc%tslot)
         call netcdf_io_check( nf90_get_var(nc%id, nc%time_varid, self%t, start=[tslot]), &
                                  "netcdf_io_read_hdr_system nf90_getvar time_varid"  )
         if (param%lmtiny_pl) then
            call nc%get_valid_masks(plmask, tpmask, plmmask, param%GMTINY)
         else
            call nc%get_valid_masks(plmask, tpmask)
         end if

         status = nf90_inq_varid(nc%id, nc%npl_varname, nc%npl_varid)
         if (status == NF90_NOERR) then
            call netcdf_io_check( nf90_get_var(nc%id, nc%npl_varid,  self%pl%nbody, start=[tslot]), &
                                  "netcdf_io_read_hdr_system nf90_getvar npl_varid"  )
         else
            self%pl%nbody = count(plmask(:))
         end if

         status = nf90_inq_varid(nc%id, nc%ntp_varname, nc%ntp_varid)
         if (status == NF90_NOERR) then
            call netcdf_io_check( nf90_get_var(nc%id, nc%ntp_varid,  self%tp%nbody, start=[tslot]), &
                                  "netcdf_io_read_hdr_system nf90_getvar ntp_varid"  )
         else
            self%tp%nbody = count(tpmask(:))
         end if

         if (param%lmtiny_pl) then
            status = nf90_inq_varid(nc%id, nc%nplm_varname, nc%nplm_varid)
            if (status == NF90_NOERR) then
               call netcdf_io_check( nf90_get_var(nc%id, nc%nplm_varid,  self%pl%nplm, start=[tslot]), &
                                  "netcdf_io_read_hdr_system nf90_getvar nplm_varid"  )
            else
               self%pl%nplm = count(plmmask(:))
            end if
         end if

         if (param%lenergy) then
            status = nf90_inq_varid(nc%id, nc%ke_orb_varname, nc%KE_orb_varid)
            if (status == NF90_NOERR) call netcdf_io_check( nf90_get_var(nc%id, nc%KE_orb_varid, self%ke_orbit, start=[tslot]), &
                                  "netcdf_io_read_hdr_system nf90_getvar KE_orb_varid"  )
            status = nf90_inq_varid(nc%id, nc%ke_spin_varname, nc%KE_spin_varid)
            if (status == NF90_NOERR) call netcdf_io_check( nf90_get_var(nc%id, nc%KE_spin_varid, self%ke_spin, start=[tslot]), &
                                  "netcdf_io_read_hdr_system nf90_getvar KE_spin_varid"  )
            status = nf90_inq_varid(nc%id, nc%pe_varname, nc%PE_varid)
            if (status == NF90_NOERR) call netcdf_io_check( nf90_get_var(nc%id, nc%PE_varid, self%pe, start=[tslot]), &
                                  "netcdf_io_read_hdr_system nf90_getvar PE_varid"  )
            status = nf90_inq_varid(nc%id, nc%be_varname, nc%BE_varid)
            if (status == NF90_NOERR) call netcdf_io_check( nf90_get_var(nc%id, nc%BE_varid, self%be, start=[tslot]), &
                                  "netcdf_io_read_hdr_system nf90_getvar BE_varid"  )
            status = nf90_inq_varid(nc%id, nc%te_varname, nc%TE_varid)
            if (status == NF90_NOERR) call netcdf_io_check( nf90_get_var(nc%id, nc%TE_varid, self%te, start=[tslot]), &
                                  "netcdf_io_read_hdr_system nf90_getvar TE_varid"  )
            status = nf90_inq_varid(nc%id, nc%L_orbit_varname, nc%L_orbit_varid)
            if (status == NF90_NOERR) call netcdf_io_check( nf90_get_var(nc%id, nc%L_orbit_varid, self%L_orbit(:), & 
                                                                         start=[1,tslot], count=[NDIM,1]), &
                                  "netcdf_io_read_hdr_system nf90_getvar L_orbit_varid"  )
            status = nf90_inq_varid(nc%id, nc%L_spin_varname, nc%L_spin_varid)
            if (status == NF90_NOERR) call netcdf_io_check( nf90_get_var(nc%id, nc%L_spin_varid, self%L_spin(:), start=[1,tslot], &
                                                            count=[NDIM,1]), &
                                  "netcdf_io_read_hdr_system nf90_getvar L_spin_varid"  )
            status = nf90_inq_varid(nc%id, nc%L_escape_varname, nc%L_escape_varid)
            if (status == NF90_NOERR) call netcdf_io_check( nf90_get_var(nc%id, nc%L_escape_varid, self%L_escape(:), &
                                                                         start=[1, tslot], count=[NDIM,1]), &
                                  "netcdf_io_read_hdr_system nf90_getvar L_escape_varid"  )
            status = nf90_inq_varid(nc%id, nc%E_collisions_varname, nc%E_collisions_varid)
            if (status == NF90_NOERR) call netcdf_io_check( nf90_get_var(nc%id, nc%E_collisions_varid, self%E_collisions, &
                                                                         start=[tslot]), &
                                  "netcdf_io_read_hdr_system nf90_getvar E_collisions_varid"  )
            status = nf90_inq_varid(nc%id, nc%E_untracked_varname, nc%E_untracked_varid)
            if (status == NF90_NOERR) call netcdf_io_check( nf90_get_var(nc%id, nc%E_untracked_varid, self%E_untracked, &
                                                                          start=[tslot]), &
                                  "netcdf_io_read_hdr_system nf90_getvar E_untracked_varid"  )
            status = nf90_inq_varid(nc%id, nc%GMescape_varname, nc%GMescape_varid)
            if (status == NF90_NOERR)  call netcdf_io_check( nf90_get_var(nc%id, nc%GMescape_varid, self%GMescape, start=[tslot]), &
                                  "netcdf_io_read_hdr_system nf90_getvar GMescape_varid"  )
         end if

      end associate

      return
   end subroutine swiftest_io_netcdf_read_hdr_system


   module subroutine swiftest_io_netcdf_read_particle_info_system(self, nc, param, plmask, tpmask)
      !! author: Carlisle A. Wishard, Dana Singh, and David A. Minton
      !!
      !! Reads particle information metadata from file
      implicit none
      ! Arguments
      class(swiftest_nbody_system),      intent(inout) :: self   !! Swiftest nbody system object
      class(swiftest_netcdf_parameters), intent(inout) :: nc     !! Parameters used to identify a particular NetCDF dataset
      class(swiftest_parameters),        intent(inout) :: param  !! Current run configuration parameters
      logical, dimension(:),             intent(in)    :: plmask !! Logical array indicating which index values belong to massive 
                                                                 !!     bodies
      logical, dimension(:),             intent(in)    :: tpmask !! Logical array indicating which index values belong to test 
                                                                 !!     particles

      ! Internals
      integer(I4B)                                :: i, idmax, status
      real(DP), dimension(:), allocatable         :: rtemp
      real(DP), dimension(:,:), allocatable       :: vectemp
      integer(I4B), dimension(:), allocatable     :: itemp
      character(len=NAMELEN), dimension(:), allocatable :: ctemp
      integer(I4B), dimension(:), allocatable     :: plind, tpind

      ! This string of spaces of length NAMELEN is used to clear out any old data left behind inside the string variables
      idmax = size(plmask)
      allocate(rtemp(idmax))
      allocate(vectemp(NDIM,idmax))
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

         call netcdf_io_check( nf90_get_var(nc%id, nc%id_varid, itemp), "netcdf_io_read_particle_info_system nf90_getvar id_varid")
         cb%id = itemp(1)
         pl%id(:) = pack(itemp, plmask)
         tp%id(:) = pack(itemp, tpmask)
         cb%id = 0
         pl%id(:) = pack([(i,i=0,idmax-1)],plmask)
         tp%id(:) = pack([(i,i=0,idmax-1)],tpmask)

         call netcdf_io_check( nf90_get_var(nc%id, nc%name_varid, ctemp, count=[NAMELEN, idmax]), &
                                  "netcdf_io_read_particle_info_system nf90_getvar name_varid"  )
         call cb%info%set_value(name=ctemp(1))
         do i = 1, npl
            call pl%info(i)%set_value(name=ctemp(plind(i)))
         end do
         do i = 1, ntp
            call tp%info(i)%set_value(name=ctemp(tpind(i)))
         end do

         status = nf90_get_var(nc%id, nc%ptype_varid, ctemp, count=[NAMELEN, idmax])
         if (status /= NF90_NOERR) then ! Set default particle types
            call cb%info%set_value(particle_type=CB_TYPE_NAME)

            ! Handle semi-interacting bodies in SyMBA
            do i = 1, npl
               if (param%lmtiny_pl .and. (pl%Gmass(i) < param%GMTINY)) then
                  call pl%info(i)%set_value(particle_type=PL_TINY_TYPE_NAME)
               else
                  call pl%info(i)%set_value(particle_type=PL_TYPE_NAME)
               end if
            end do
            do i = 1, ntp
               call tp%info(i)%set_value(particle_type=TP_TYPE_NAME)
            end do
         else ! Use particle types defined in input file
            call cb%info%set_value(particle_type=ctemp(1))
            do i = 1, npl
               call pl%info(i)%set_value(particle_type=ctemp(plind(i)))
            end do
            do i = 1, ntp
               call tp%info(i)%set_value(particle_type=ctemp(tpind(i)))
            end do
         end if

         call cb%info%set_value(status="ACTIVE")

         if (param%lclose) then
            status = nf90_inq_varid(nc%id, nc%origin_type_varname, nc%origin_type_varid)
            if (status == NF90_NOERR) then
               call netcdf_io_check( nf90_get_var(nc%id, nc%origin_type_varid, ctemp, count=[NAMELEN, idmax]), &
                                  "netcdf_io_read_particle_info_system nf90_getvar origin_type_varid"  )
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

            status = nf90_inq_varid(nc%id, nc%origin_time_varname, nc%origin_time_varid)
            if (status == NF90_NOERR) then
               call netcdf_io_check( nf90_get_var(nc%id, nc%origin_time_varid, rtemp), &
                                  "netcdf_io_read_particle_info_system nf90_getvar origin_time_varid"  )
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

            status = nf90_inq_varid(nc%id, nc%origin_rh_varname, nc%origin_rh_varid)
            if (status == NF90_NOERR) then
               call netcdf_io_check( nf90_get_var(nc%id, nc%origin_rh_varid, vectemp(:,:)), &
                                  "netcdf_io_read_particle_info_system nf90_getvar origin_rh_varid"  )
            else if ((param%out_form == "XV") .or. (param%out_form == "XVEL")) then
               call netcdf_io_check( nf90_get_var(nc%id, nc%rh_varid, vectemp(:,:)), &
                                  "netcdf_io_read_particle_info_system nf90_getvar rh_varid"  )
            else 
               vectemp(:,:) = 0._DP
            end if 

            do i = 1, npl
               call pl%info(i)%set_value(origin_rh=vectemp(:,plind(i)))
            end do
            do i = 1, ntp
               call tp%info(i)%set_value(origin_rh=vectemp(:,tpind(i)))
            end do

            status = nf90_inq_varid(nc%id, nc%origin_vh_varname, nc%origin_vh_varid)
            if (status == NF90_NOERR) then
               call netcdf_io_check( nf90_get_var(nc%id, nc%origin_vh_varid, vectemp(:,:)), &
                                  "netcdf_io_read_particle_info_system nf90_getvar origin_vh_varid"  )
            else if ((param%out_form == "XV") .or. (param%out_form == "XVEL")) then
               call netcdf_io_check( nf90_get_var(nc%id, nc%vh_varid, vectemp(:,:)), &
                                  "netcdf_io_read_particle_info_system nf90_getvar vh_varid"  )
            else
               vectemp(:,:) = 0._DP
            end if 
            
            do i = 1, npl
               call pl%info(i)%set_value(origin_vh=vectemp(:,plind(i)))
            end do
            do i = 1, ntp
               call tp%info(i)%set_value(origin_vh=vectemp(:,tpind(i)))
            end do

            status = nf90_inq_varid(nc%id, nc%collision_id_varname, nc%collision_id_varid)
            if (status == NF90_NOERR) then
               call netcdf_io_check( nf90_get_var(nc%id, nc%collision_id_varid, itemp), &
                                  "netcdf_io_read_particle_info_system nf90_getvar collision_id_varid"  )
               do i = 1, npl
                  call pl%info(i)%set_value(collision_id=itemp(plind(i)))
               end do
               do i = 1, ntp
                  call tp%info(i)%set_value(collision_id=itemp(tpind(i)))
               end do
            end if

            status = nf90_inq_varid(nc%id, nc%discard_time_varname, nc%discard_time_varid)
            if (status == NF90_NOERR) then
               call netcdf_io_check( nf90_get_var(nc%id, nc%discard_time_varid, rtemp), &
                                  "netcdf_io_read_particle_info_system nf90_getvar discard_time_varid"  )
               call cb%info%set_value(discard_time=rtemp(1))
               do i = 1, npl
                  call pl%info(i)%set_value(discard_time=rtemp(plind(i)))
               end do
               do i = 1, ntp
                  call tp%info(i)%set_value(discard_time=rtemp(tpind(i)))
               end do
            end if

            status = nf90_inq_varid(nc%id, nc%discard_rh_varname, nc%discard_rh_varid)
            if (status == NF90_NOERR) then
               call netcdf_io_check( nf90_get_var(nc%id, nc%discard_rh_varid, vectemp(:,:)), &
                                  "netcdf_io_read_particle_info_system nf90_getvar discard_rh_varid"  )
               do i = 1, npl
                  call pl%info(i)%set_value(discard_rh=vectemp(:,plind(i)))
               end do
               do i = 1, ntp
                  call tp%info(i)%set_value(discard_rh=vectemp(:,tpind(i)))
               end do
            end if

            status = nf90_inq_varid(nc%id, nc%discard_vh_varname, nc%discard_vh_varid)
            if (status == NF90_NOERR) then
               call netcdf_io_check( nf90_get_var(nc%id, nc%discard_vh_varid, vectemp(:,:)), &
                                  "netcdf_io_read_particle_info_system nf90_getvar discard_vh_varid"  )
               do i = 1, npl
                  call pl%info(i)%set_value(discard_vh=vectemp(:,plind(i)))
               end do
               do i = 1, ntp
                  call tp%info(i)%set_value(discard_vh=vectemp(:,tpind(i)))
               end do
            end if
         end if

      end associate

      return
   end subroutine swiftest_io_netcdf_read_particle_info_system


   module subroutine swiftest_io_netcdf_write_frame_body(self, nc, param)
      !! author: Carlisle A. Wishard, Dana Singh, and David A. Minton
      !!
      !! Write a frame of output of either test particle or massive body data to the binary output file
      !!    Note: If outputting to orbital elements, but sure that the conversion is done prior to calling this method
      implicit none
      ! Arguments
      class(swiftest_body),              intent(in)    :: self  !! Swiftest base object
      class(swiftest_netcdf_parameters), intent(inout) :: nc    !! Parameters used to for writing a NetCDF dataset to file
      class(swiftest_parameters),        intent(inout) :: param !! Current run configuration parameters 
      ! Internals
      integer(I4B)                              :: i, j, idslot, old_mode, tmp
      integer(I4B), dimension(:), allocatable   :: ind
      real(DP), dimension(NDIM)                 :: vh !! Temporary variable to store heliocentric velocity values when converting 
                                                      !! from pseudovelocity in GR-enabled runs
      real(DP)                                  :: a, e, inc, omega, capom, capm, varpi, lam, f, cape, capf
#ifdef COARRAY
      integer(I4B) :: ntp
      logical, dimension(:), allocatable        :: tpmask, plmask
#endif

      call self%write_info(nc, param)

      call netcdf_io_check( nf90_set_fill(nc%id, NF90_NOFILL, old_mode), "netcdf_io_write_frame_body nf90_set_fill"  )
      select type(self)
      class is (swiftest_body)
      select type (param)
      class is (swiftest_parameters)
         associate(n => self%nbody, tslot => nc%tslot)
            if (n == 0) return

            call util_sort(self%id(1:n), ind)

            do i = 1, n
               j = ind(i)
               call nc%find_idslot(self%id(j), idslot) 
               ! Convert from pseudovelocity to heliocentric without replacing the current value of pseudovelocity 
               if (param%lgr) call swiftest_gr_pseudovel2vel(param, self%mu(j), self%rh(:, j), self%vh(:, j), vh(:))

               if ((param%out_form == "XV") .or. (param%out_form == "XVEL")) then
                  call netcdf_io_check( nf90_put_var(nc%id, nc%rh_varid, self%rh(:, j), start=[1,idslot, tslot], count=[NDIM,1,1]),&
                                  "netcdf_io_write_frame_body nf90_put_var rh_varid"  )
                  if (param%lgr) then ! Convert from pseudovelocity to heliocentric without replacing the current value of 
                                      !  pseudovelocity
                     call netcdf_io_check( nf90_put_var(nc%id, nc%vh_varid, vh(:), start=[1,idslot, tslot], count=[NDIM,1,1]), &
                                  "netcdf_io_write_frame_body nf90_put_var vh_varid"  )
                     call netcdf_io_check( nf90_put_var(nc%id, nc%gr_pseudo_vh_varid, self%vh(:, j), start=[1,idslot, tslot], &
                                                        count=[NDIM,1,1]), &
                                  "netcdf_io_write_frame_body nf90_put_var gr_pseudo_vhx_varid"  )

                  else
                     call netcdf_io_check( nf90_put_var(nc%id, nc%vh_varid, self%vh(:, j), start=[1,idslot, tslot], &
                                                        count=[NDIM,1,1]), &
                                  "netcdf_io_write_frame_body nf90_put_var vh_varid"  )
                  end if
               end if

               if ((param%out_form == "EL") .or. (param%out_form == "XVEL")) then
                  if (param%lgr) then ! For GR-enabled runs, use the true value of velocity computed above
                     call swiftest_orbel_xv2el(self%mu(j), self%rh(1,j), self%rh(2,j), self%rh(3,j), &
                                       vh(1), vh(2), vh(3), &
                                       a, e, inc, capom, omega, capm, varpi, lam, f, cape, capf)
                  else !! For non-GR runs just convert from the velocity we have
                     call swiftest_orbel_xv2el(self%mu(j), self%rh(1,j), self%rh(2,j), self%rh(3,j), &
                                       self%vh(1,j), self%vh(2,j), self%vh(3,j), &
                                       a, e, inc, capom, omega, capm, varpi, lam, f, cape, capf)
                  end if
                  call netcdf_io_check( nf90_put_var(nc%id, nc%a_varid, a, start=[idslot, tslot]), &
                                  "netcdf_io_write_frame_body nf90_put_var body a_varid"  )
                  call netcdf_io_check( nf90_put_var(nc%id, nc%e_varid, e, start=[idslot, tslot]), &
                                  "netcdf_io_write_frame_body nf90_put_var body e_varid"  )
                  call netcdf_io_check( nf90_put_var(nc%id, nc%inc_varid, inc * RAD2DEG, start=[idslot, tslot]), &
                                  "netcdf_io_write_frame_body nf90_put_var body inc_varid"  )
                  call netcdf_io_check( nf90_put_var(nc%id, nc%capom_varid, capom * RAD2DEG, start=[idslot, tslot]), &
                                  "netcdf_io_write_frame_body nf90_put_var body capom_varid"  )
                  call netcdf_io_check( nf90_put_var(nc%id, nc%omega_varid, omega * RAD2DEG, start=[idslot, tslot]), &
                                  "netcdf_io_write_frame_body nf90_put_var body omega_varid"  )
                  call netcdf_io_check( nf90_put_var(nc%id, nc%capm_varid, capm * RAD2DEG, start=[idslot, tslot]), &
                                  "netcdf_io_write_frame_body nf90_put_var body capm_varid"  ) 
                  call netcdf_io_check( nf90_put_var(nc%id, nc%varpi_varid, varpi * RAD2DEG, start=[idslot, tslot]), &
                                  "netcdf_io_write_frame_body nf90_put_var body varpi_varid"  ) 
                  call netcdf_io_check( nf90_put_var(nc%id, nc%lam_varid, lam * RAD2DEG, start=[idslot, tslot]), &
                                  "netcdf_io_write_frame_body nf90_put_var body lam_varid"  ) 
                  call netcdf_io_check( nf90_put_var(nc%id, nc%f_varid, f * RAD2DEG, start=[idslot, tslot]), &
                                  "netcdf_io_write_frame_body nf90_put_var body f_varid"  ) 
                  if (e < 1.0_DP) then
                     call netcdf_io_check( nf90_put_var(nc%id, nc%cape_varid, cape * RAD2DEG, start=[idslot, tslot]), &
                                  "netcdf_io_write_frame_body nf90_put_var body cape_varid"  ) 
                  else if (e > 1.0_DP) then
                     call netcdf_io_check( nf90_put_var(nc%id, nc%cape_varid, capf * RAD2DEG, start=[idslot, tslot]), &
                                  "netcdf_io_write_frame_body nf90_put_var body (capf) cape_varid"  ) 
                  end if
               end if

               select type(self)  
               class is (swiftest_pl)  ! Additional output if the passed polymorphic object is a massive body
                  call netcdf_io_check( nf90_put_var(nc%id, nc%Gmass_varid, self%Gmass(j), start=[idslot, tslot]), &
                                  "netcdf_io_write_frame_body nf90_put_var body Gmass_varid"  )
                  call netcdf_io_check( nf90_put_var(nc%id, nc%mass_varid, self%mass(j), start=[idslot, tslot]), &
                                  "netcdf_io_write_frame_body nf90_put_var body mass_varid"  )
                  if (param%lrhill_present) then
                     call netcdf_io_check( nf90_put_var(nc%id, nc%rhill_varid, self%rhill(j), start=[idslot, tslot]), &
                                  "netcdf_io_write_frame_body nf90_put_var body rhill_varid"  )
                  end if
                  if (param%lclose) call netcdf_io_check( nf90_put_var(nc%id, nc%radius_varid, self%radius(j), &
                                                                       start=[idslot, tslot]), &
                                  "netcdf_io_write_frame_body nf90_put_var body radius_varid"  )
                  if (param%lrotation) then
                     call netcdf_io_check( nf90_put_var(nc%id, nc%Ip_varid, self%Ip(:, j), start=[1,idslot, tslot], &
                                                        count=[NDIM,1,1]), &
                                  "netcdf_io_write_frame_body nf90_put_var body Ip_varid"  )
                     call netcdf_io_check( nf90_put_var(nc%id, nc%rot_varid, self%rot(:, j) * RAD2DEG, start=[1,idslot, tslot], &
                                                        count=[NDIM,1,1]), &
                                  "netcdf_io_write_frame_body nf90_put_var body rotx_varid"  )
                  end if
                  ! if (param%ltides) then
                  !    call netcdf_io_check( nf90_put_var(nc%id, nc%k2_varid, self%k2(j), start=[idslot, tslot]), &
                  !                "netcdf_io_write_frame_body nf90_put_var body k2_varid"  )
                  !    call netcdf_io_check( nf90_put_var(nc%id, nc%Q_varid, self%Q(j), start=[idslot, tslot]), &
                  !                "netcdf_io_write_frame_body nf90_put_var body Q_varid"  )
                  ! end if
               class is (swiftest_tp)
                  call netcdf_io_check( nf90_put_var(nc%id, nc%Gmass_varid, 0.0_DP, start=[idslot, tslot]), &
                                  "netcdf_io_write_frame_body nf90_put_var body Gmass_varid"  )
                  call netcdf_io_check( nf90_put_var(nc%id, nc%mass_varid, 0.0_DP, start=[idslot, tslot]), &
                                  "netcdf_io_write_frame_body nf90_put_var body mass_varid"  )
                  if (param%lrhill_present) then
                     call netcdf_io_check( nf90_put_var(nc%id, nc%rhill_varid, 0.0_DP, start=[idslot, tslot]), &
                                  "netcdf_io_write_frame_body nf90_put_var body rhill_varid"  )
                  end if
                  if (param%lclose) call netcdf_io_check( nf90_put_var(nc%id, nc%radius_varid, 0.0_DP, start=[idslot, tslot]), &
                                  "netcdf_io_write_frame_body nf90_put_var body radius_varid"  )
                  if (param%lrotation) then
                     call netcdf_io_check( nf90_put_var(nc%id, nc%Ip_varid, [0.0_DP,0.0_DP,0.0_DP], start=[1,idslot, tslot], &
                                                        count=[NDIM,1,1]), &
                                  "netcdf_io_write_frame_body nf90_put_var body Ip_varid"  )
                     call netcdf_io_check( nf90_put_var(nc%id, nc%rot_varid, [0.0_DP,0.0_DP,0.0_DP], start=[1,idslot, tslot], &
                                                        count=[NDIM,1,1]), &
                                  "netcdf_io_write_frame_body nf90_put_var body rotx_varid"  )
                  end if
               end select
            end do
         end associate
      end select
      end select
#ifdef COARRAY
      select type(self)
      class is (swiftest_tp)
         call nc%get_valid_masks(plmask, tpmask)
         ntp = count(tpmask(:))  
         call netcdf_io_check( nf90_put_var(nc%id, nc%ntp_varid, ntp, start=[nc%tslot]), &
                                  "netcdf_io_write_frame_body nf90_put_var ntp_varid"  )
      end select
#endif   

      call netcdf_io_check( nf90_set_fill(nc%id, old_mode, tmp), &
                                  "netcdf_io_write_frame_body nf90_set_fill old_mode"  )

      return
   end subroutine swiftest_io_netcdf_write_frame_body


   module subroutine swiftest_io_netcdf_write_frame_cb(self, nc, param)
      !! author: Carlisle A. Wishard, Dana Singh, and David A. Minton
      !!
      !! Write a frame of output of the central body
      implicit none
      ! Arguments
      class(swiftest_cb),                intent(in)    :: self  !! Swiftest base object
      class(swiftest_netcdf_parameters), intent(inout) :: nc    !! Parameters used to for writing a NetCDF dataset to file
      class(swiftest_parameters),        intent(inout) :: param !! Current run configuration parameters 
      ! Internals
      integer(I4B)                              :: idslot, old_mode, tmp

      associate(tslot => nc%tslot)
         call self%write_info(nc, param)

         call netcdf_io_check( nf90_set_fill(nc%id, NF90_NOFILL, old_mode), &
                                  "swiftest_io_netcdf_write_frame_cb nf90_set_fill"  )

         call nc%find_idslot(self%id, idslot) 
         call netcdf_io_check( nf90_put_var(nc%id, nc%id_varid, self%id, start=[idslot]), &
                                  "swiftest_io_netcdf_write_frame_cb nf90_put_var cb id_varid"  )
         call netcdf_io_check( nf90_put_var(nc%id, nc%status_varid, ACTIVE, start=[idslot, tslot]), &
                                  "swiftest_io_netcdf_write_frame_cb nf90_put_var cb id_varid"  )

         call netcdf_io_check( nf90_put_var(nc%id, nc%Gmass_varid, self%Gmass, start=[idslot, tslot]), &
                                  "swiftest_io_netcdf_write_frame_cb nf90_put_var cb Gmass_varid"  )
         call netcdf_io_check( nf90_put_var(nc%id, nc%mass_varid, self%mass, start=[idslot, tslot]), &
                                  "swiftest_io_netcdf_write_frame_cb nf90_put_var cb mass_varid"  )
         if (param%lclose) call netcdf_io_check( nf90_put_var(nc%id, nc%radius_varid, self%radius, start=[idslot, tslot]), &
                                  "swiftest_io_netcdf_write_frame_cb nf90_put_var cb radius_varid"  )
         call netcdf_io_check( nf90_put_var(nc%id, nc%j2rp2_varid, self%j2rp2, start=[tslot]), &
                                  "swiftest_io_netcdf_write_frame_cb nf90_put_var cb j2rp2_varid" )
         call netcdf_io_check( nf90_put_var(nc%id, nc%j4rp4_varid, self%j4rp4, start=[tslot]), &
                                  "swiftest_io_netcdf_write_frame_cb nf90_put_var cb j4rp4_varid" )
         if (param%lrotation) then
            call netcdf_io_check( nf90_put_var(nc%id, nc%Ip_varid, self%Ip(:), start=[1, idslot, tslot], count=[NDIM,1,1]), &
                                  "swiftest_io_netcdf_write_frame_cb nf90_put_var cb Ip_varid"  )
            call netcdf_io_check( nf90_put_var(nc%id, nc%rot_varid, self%rot(:) * RAD2DEG, start=[1, idslot, tslot], &
                                               count=[NDIM,1,1]), &
                                  "swiftest_io_netcdf_write_frame_cby nf90_put_var cb rot_varid"  )
         end if

         call netcdf_io_check( nf90_set_fill(nc%id, old_mode, tmp), &
                                  "swiftest_io_netcdf_write_frame_cb nf90_set_fill old_mode"  )
      end associate

      return
   end subroutine swiftest_io_netcdf_write_frame_cb


   module subroutine swiftest_io_netcdf_write_frame_system(self, nc, param)
      !! author: The Purdue Swiftest Team - David A. Minton, Carlisle A. Wishard, Jennifer L.L. Pouplin, and Jacob R. Elliott
      !!
      !! Write a frame (header plus records for each massive body and active test particle) to a output binary file
      implicit none
      ! Arguments
      class(swiftest_nbody_system),      intent(inout) :: self  !! Swiftest nbody_system object
      class(swiftest_netcdf_parameters), intent(inout) :: nc    !! Parameters used to for writing a NetCDF dataset to file
      class(swiftest_parameters),        intent(inout) :: param !! Current run configuration parameters 

      call self%write_hdr(nc, param)
#ifdef COARRAY
      if (this_image() == 1) then
#endif
         call self%cb%write_frame(nc, param)
         call self%pl%write_frame(nc, param)
#ifdef COARRAY
      end if ! this_image() == 1
#endif
      call self%tp%write_frame(nc, param)

      return
   end subroutine swiftest_io_netcdf_write_frame_system


   module subroutine swiftest_io_netcdf_write_hdr_system(self, nc, param) 
      !! author: David A. Minton
      !!
      !! Writes header information (variables that change with time, but not particle id). 
      !! This subroutine swiftest_significantly improves the output over the original binary file, allowing us to track energy,
      !! momentum, and other quantities that previously were handled as separate output files.
      implicit none
      ! Arguments
      class(swiftest_nbody_system),      intent(in)    :: self  !! Swiftest nbody system object
      class(swiftest_netcdf_parameters), intent(inout) :: nc    !! Parameters used to for writing a NetCDF dataset to file
      class(swiftest_parameters),        intent(inout) :: param !! Current run configuration parameters
      ! Internals
      integer(I4B) :: tslot

      call nc%find_tslot(self%t, tslot)
      call netcdf_io_check( nf90_put_var(nc%id, nc%time_varid, self%t, start=[tslot]), &
                                  "netcdf_io_write_hdr_system nf90_put_var time_varid"  )
      call netcdf_io_check( nf90_put_var(nc%id, nc%npl_varid, self%pl%nbody, start=[tslot]), &
                                  "netcdf_io_write_hdr_system nf90_put_var npl_varid"  )
#ifndef COARRAY
      call netcdf_io_check( nf90_put_var(nc%id, nc%ntp_varid, self%tp%nbody, start=[tslot]), &
                                  "netcdf_io_write_hdr_system nf90_put_var ntp_varid"  )
#endif
      if (param%lmtiny_pl) call netcdf_io_check( nf90_put_var(nc%id, nc%nplm_varid, self%pl%nplm, start=[tslot]), &
                                  "netcdf_io_write_hdr_system nf90_put_var nplm_varid"  )

      if (param%lenergy) then
         call netcdf_io_check( nf90_put_var(nc%id, nc%KE_orb_varid, self%ke_orbit, start=[tslot]), &
                                  "netcdf_io_write_hdr_system nf90_put_var KE_orb_varid"  )
         call netcdf_io_check( nf90_put_var(nc%id, nc%KE_spin_varid, self%ke_spin, start=[tslot]), &
                                  "netcdf_io_write_hdr_system nf90_put_var KE_spin_varid"  )
         call netcdf_io_check( nf90_put_var(nc%id, nc%PE_varid, self%pe, start=[tslot]), &
                                  "netcdf_io_write_hdr_system nf90_put_var PE_varid"  )
         call netcdf_io_check( nf90_put_var(nc%id, nc%BE_varid, self%be, start=[tslot]), &
                                  "netcdf_io_write_hdr_system nf90_put_var BE_varid"  )
         call netcdf_io_check( nf90_put_var(nc%id, nc%TE_varid, self%te, start=[tslot]), &
                                  "netcdf_io_write_hdr_system nf90_put_var TE_varid"  )
         call netcdf_io_check( nf90_put_var(nc%id, nc%L_orbit_varid, self%L_orbit(:), start=[1,tslot], count=[NDIM,1]), &
                                  "netcdf_io_write_hdr_system nf90_put_var L_orbit_varid"  )
         call netcdf_io_check( nf90_put_var(nc%id, nc%L_spin_varid, self%L_spin(:), start=[1,tslot], count=[NDIM,1]), &
                                  "netcdf_io_write_hdr_system nf90_put_var L_spin_varid"  )
         call netcdf_io_check( nf90_put_var(nc%id, nc%L_escape_varid, self%L_escape(:), start=[1,tslot], count=[NDIM,1]), &
                                  "netcdf_io_write_hdr_system nf90_put_var L_escape_varid"  )
         call netcdf_io_check( nf90_put_var(nc%id, nc%E_collisions_varid, self%E_collisions, start=[tslot]), &
                                  "netcdf_io_write_hdr_system nf90_put_var E_collisions_varid"  )
         call netcdf_io_check( nf90_put_var(nc%id, nc%E_untracked_varid, self%E_untracked, start=[tslot]), &
                                  "netcdf_io_write_hdr_system nf90_put_var E_untracked_varid"  )
         call netcdf_io_check( nf90_put_var(nc%id, nc%GMescape_varid, self%GMescape, start=[tslot]), &
                                  "netcdf_io_write_hdr_system nf90_put_var GMescape_varid"  )
      end if

      return
   end subroutine swiftest_io_netcdf_write_hdr_system


   module subroutine swiftest_io_netcdf_write_info_body(self, nc, param)
      !! author: Carlisle A. Wishard, Dana Singh, and David A. Minton
      !!
      !! Write all current particle information metadata to file
      implicit none
      ! Arguments
      class(swiftest_body),              intent(in)    :: self  !! Swiftest particle object
      class(swiftest_netcdf_parameters), intent(inout) :: nc      !! Parameters used to identify a particular NetCDF dataset
      class(swiftest_parameters),        intent(inout) :: param !! Current run configuration parameters
      ! Internals
      integer(I4B)                              :: i, j, idslot, old_mode, tmp
      integer(I4B), dimension(:), allocatable   :: ind
      character(len=NAMELEN) :: charstring
      character(len=NAMELEN), dimension(self%nbody) :: origin_type

      call netcdf_io_check( nf90_set_fill(nc%id, NF90_NOFILL, old_mode), &
                                  "netcdf_io_write_info_body nf90_set_fill NF90_NOFILL"  )

      select type(self)
         class is (swiftest_body)
         associate(n => self%nbody, tslot => nc%tslot)
            if (n == 0) return
            call util_sort(self%id(1:n), ind)
            call nc%get_idvals()

            do i = 1, n
               j = ind(i)
               call nc%find_idslot(self%id(j), idslot) 
               call netcdf_io_check( nf90_put_var(nc%id, nc%id_varid, self%id(j), start=[idslot]), &
                                  "netcdf_io_write_info_body nf90_put_var id_varid"  )
               call netcdf_io_check( nf90_put_var(nc%id, nc%status_varid, self%status(j), start=[idslot,tslot]), &
                                  "netcdf_io_write_info_body nf90_put_var status_varid"  )

               charstring = trim(adjustl(self%info(j)%name))
               call netcdf_io_check( nf90_put_var(nc%id, nc%name_varid, charstring, start=[1, idslot], count=[NAMELEN, 1]), &
                                  "netcdf_io_write_info_body nf90_put_var name_varid"  )

               charstring = trim(adjustl(self%info(j)%particle_type))
               call netcdf_io_check( nf90_put_var(nc%id, nc%ptype_varid, charstring, start=[1, idslot], count=[NAMELEN, 1]), &
                                  "netcdf_io_write_info_body nf90_put_var particle_type_varid"  )

               if (param%lclose) then
                  charstring = trim(adjustl(self%info(j)%origin_type))
                  origin_type(i) = charstring
                  call netcdf_io_check( nf90_put_var(nc%id, nc%origin_type_varid, charstring, start=[1, idslot], &
                                                     count=[NAMELEN, 1]), &
                                  "netcdf_io_write_info_body nf90_put_var origin_type_varid"  )
                  call netcdf_io_check( nf90_put_var(nc%id, nc%origin_time_varid,  self%info(j)%origin_time,  start=[idslot]), &
                                  "netcdf_io_write_info_body nf90_put_var origin_time_varid"  )
                  call netcdf_io_check( nf90_put_var(nc%id, nc%origin_rh_varid,    self%info(j)%origin_rh(:), start=[1,idslot], &
                                                     count=[NDIM,1]), &
                                  "netcdf_io_write_info_body nf90_put_var origin_rh_varid"  )
                  call netcdf_io_check( nf90_put_var(nc%id, nc%origin_vh_varid,    self%info(j)%origin_vh(:), start=[1,idslot], &
                                                     count=[NDIM,1]), &
                                  "netcdf_io_write_info_body nf90_put_var origin_vh_varid"  )
   
                  call netcdf_io_check( nf90_put_var(nc%id, nc%collision_id_varid, self%info(j)%collision_id, start=[idslot]), &
                                  "netcdf_io_write_info_body nf90_put_var collision_id_varid"  )
                  call netcdf_io_check( nf90_put_var(nc%id, nc%discard_time_varid, self%info(j)%discard_time, start=[idslot]), &
                                  "netcdf_io_write_info_body nf90_put_var discard_time_varid"  )
                  call netcdf_io_check( nf90_put_var(nc%id, nc%discard_rh_varid,   self%info(j)%discard_rh(:), start=[1,idslot], &
                                                     count=[NDIM,1]), &
                                  "netcdf_io_write_info_body nf90_put_var discard_rh_varid"  )
                  call netcdf_io_check( nf90_put_var(nc%id, nc%discard_vh_varid,   self%info(j)%discard_vh(:), start=[1,idslot], &
                                                     count=[NDIM,1]), &
                                  "netcdf_io_write_info_body nf90_put_var discard_vh_varid"  )
               end if

            end do
         end associate
      end select

      call netcdf_io_check( nf90_set_fill(nc%id, old_mode, tmp) )
      return
   end subroutine swiftest_io_netcdf_write_info_body


   module subroutine swiftest_io_netcdf_write_info_cb(self, nc, param)
      !! author: Carlisle A. Wishard, Dana Singh, and David A. Minton
      !!
      !! Write the central body info to file
      implicit none
      class(swiftest_cb),               intent(in)    :: self  !! Swiftest particle object
      class(swiftest_netcdf_parameters), intent(inout) :: nc      !! Parameters used to identify a particular NetCDF dataset
      class(swiftest_parameters),           intent(inout) :: param !! Current run configuration parameters
      ! Internals
      integer(I4B)                              :: idslot, old_mode, tmp
      character(len=NAMELEN) :: charstring

      ! This string of spaces of length NAMELEN is used to clear out any old data left behind inside the string variables
      call netcdf_io_check( nf90_set_fill(nc%id, NF90_NOFILL, old_mode), &
                                  "netcdf_io_write_info_cb nf90_set_fill NF90_NOFILL"  )

      call nc%find_idslot(self%id, idslot) 

      call netcdf_io_check( nf90_put_var(nc%id, nc%id_varid, self%id, start=[idslot]), &
                                  "netcdf_io_write_info_cb nf90_put_var id_varid"  )

      charstring = trim(adjustl(self%info%name))
      call netcdf_io_check( nf90_put_var(nc%id, nc%name_varid, charstring, start=[1, idslot], count=[NAMELEN, 1]), &
                                  "netcdf_io_write_info_cb nf90_put_var name_varid"  )

      charstring = trim(adjustl(self%info%particle_type))
      call netcdf_io_check( nf90_put_var(nc%id, nc%ptype_varid, charstring, start=[1, idslot], count=[NAMELEN, 1]), &
                                  "netcdf_io_write_info_cb nf90_put_var ptype_varid"  )

      if (param%lclose) then
         charstring = trim(adjustl(self%info%origin_type))
         call netcdf_io_check( nf90_put_var(nc%id, nc%origin_type_varid, charstring, start=[1, idslot], count=[NAMELEN, 1]), &
                                  "netcdf_io_write_info_body nf90_put_var cb origin_type_varid"  )

         call netcdf_io_check( nf90_put_var(nc%id, nc%origin_time_varid, self%info%origin_time, start=[idslot]), &
                                  "netcdf_io_write_info_body nf90_put_var cb origin_time_varid"  )
         call netcdf_io_check( nf90_put_var(nc%id, nc%origin_rh_varid, self%info%origin_rh(:), start=[1, idslot], count=[NDIM,1]), &
                                  "netcdf_io_write_info_body nf90_put_var cb origin_rh_varid"  )
         call netcdf_io_check( nf90_put_var(nc%id, nc%origin_vh_varid, self%info%origin_vh(:), start=[1, idslot], count=[NDIM,1]), &
                                  "netcdf_io_write_info_body nf90_put_var cb origin_vh_varid"  )

         call netcdf_io_check( nf90_put_var(nc%id, nc%collision_id_varid, self%info%collision_id, start=[idslot]), &
                                  "netcdf_io_write_info_body nf90_put_var cb collision_id_varid"  )
         call netcdf_io_check( nf90_put_var(nc%id, nc%discard_time_varid, self%info%discard_time, start=[idslot]), &
                                  "netcdf_io_write_info_body nf90_put_var cb discard_time_varid"  )
         call netcdf_io_check( nf90_put_var(nc%id, nc%discard_rh_varid, self%info%discard_rh(:), start=[1, idslot], &
                                            count=[NDIM,1]), &
                                  "netcdf_io_write_info_body nf90_put_var cb discard_rh_varid"  )
         call netcdf_io_check( nf90_put_var(nc%id, nc%discard_vh_varid, self%info%discard_vh(:), start=[1, idslot], &
                                            count=[NDIM,1]), &
                                  "netcdf_io_write_info_body nf90_put_var cb discard_vh_varid"  )
      end if
      call netcdf_io_check( nf90_set_fill(nc%id, old_mode, tmp) )

      return
   end subroutine swiftest_io_netcdf_write_info_cb


   module subroutine swiftest_io_remove_nul_char(string)
      !! author: David A. Minton
      !!
      !! Remove spaces and trailing NUL characters that are introduced by NetCDF strings
      implicit none
      ! Arguments
      character(len=*), intent(inout) :: string !! String to make upper case
      ! Internals
      integer(I4B) :: pos
  
      pos = index(string, achar(0)) - 1
      if (pos > 0) then
         string = trim(adjustl(string(1:pos)))
      else
         string = trim(adjustl(string))
      end if

      return
   end subroutine swiftest_io_remove_nul_char


   module subroutine swiftest_io_param_reader(self, unit, iotype, v_list, iostat, iomsg) 
      !! author: The Purdue Swiftest Team - David A. Minton, Carlisle A. Wishard, Jennifer L.L. Pouplin, and Jacob R. Elliott
      !!
      !! Read in parameters for the integration
      !! as the newline characters are ignored in the input file when compiled in ifort.
      !!
      !! Adapted from David E. Kaufmann's Swifter routine io_init_param.f90
      !! Adapted from Martin Duncan's Swift routine io_init_param.f
      implicit none
      ! Arguments
      class(swiftest_parameters), intent(inout) :: self       !! Collection of parameters
      integer(I4B),               intent(in)    :: unit       !! File unit number
      character(len=*),           intent(in)    :: iotype     !! Dummy argument passed to the  input/output procedure contains the 
                                                              !!    text from the char-literal-constant, prefixed with DT. 
                                                              !!    If you do not include a char-literal-constant, the iotype 
                                                              !!    argument contains only DT.
      character(len=*),           intent(in)    :: v_list(:)  !! The first element passes the integrator code to the reader
      integer(I4B),               intent(out)   :: iostat     !! IO status code
      character(len=*),           intent(inout) :: iomsg      !! Message to pass if iostat /= 0
      ! Internals
      logical                        :: tstart_set = .false.              !! Is the final time set in the input file?
      logical                        :: tstop_set = .false.               !! Is the final time set in the input file?
      logical                        :: dt_set = .false.                  !! Is the step size set in the input file?
      integer(I4B)                   :: ilength, ifirst, ilast, i         !! Variables used to parse input file
      character(STRMAX)              :: line                              !! Line of the input file
      character(len=:), allocatable  :: line_trim,param_name, param_value !! Strings used to parse the param file
      character(*),parameter         :: linefmt = '(A)'                   !! Format code for simple text string
      integer(I4B)                   :: nseeds, nseeds_from_file
      logical                        :: seed_set = .false.      !! Is the random seed set in the input file?
      real(DP)                       :: tratio, y
#ifdef COARRAY
      type(swiftest_parameters), codimension[*], save :: coparam
     
   if (this_image() == 1) then
      coparam = self
      associate(param => coparam) 
#else
      associate(param => self) 
#endif
         ! Parse the file line by line, extracting tokens then matching them up with known parameters if possible
         call random_seed(size = nseeds)
         if (allocated(param%seed)) deallocate(param%seed)
         allocate(param%seed(nseeds))
         open(unit = unit, file = param%param_file_name, status = 'old', err = 667, iomsg = iomsg)
         do
            read(unit = unit, fmt = linefmt, end = 1, err = 667, iomsg = iomsg) line
            line_trim = trim(adjustl(line))
            ilength = len(line_trim)
            if ((ilength /= 0)) then 
               ifirst = 1
               ! Read the pair of tokens. The first one is the parameter name, the second is the value.
               param_name = swiftest_io_get_token(line_trim, ifirst, ilast, iostat)
               if (param_name == '') cycle ! No parameter name (usually because this line is commented out)
               call swiftest_io_toupper(param_name)
               ifirst = ilast + 1
               param_value = swiftest_io_get_token(line_trim, ifirst, ilast, iostat)
               select case (param_name)
               case ("T0")
                  read(param_value, *, err = 667, iomsg = iomsg) param%t0
               case ("TSTART")
                  read(param_value, *, err = 667, iomsg = iomsg) param%tstart
                  tstart_set = .true.                  
               case ("TSTOP")
                  read(param_value, *, err = 667, iomsg = iomsg) param%tstop
                  tstop_set = .true.
               case ("DT")
                  read(param_value, *, err = 667, iomsg = iomsg) param%dt
                  dt_set = .true.
               case ("CB_IN")
                  param%incbfile = param_value
               case ("PL_IN")
                  param%inplfile = param_value
               case ("TP_IN")
                  param%intpfile = param_value
               case ("NC_IN")
                  param%nc_in = param_value
               case ("IN_TYPE")
                  call swiftest_io_toupper(param_value)
                  param%in_type = param_value
               case ("IN_FORM")
                  call swiftest_io_toupper(param_value)
                  param%in_form = param_value
               case ("ISTEP_OUT")
                  read(param_value, *) param%istep_out
               case ("NSTEP_OUT")
                  read(param_value, *) param%nstep_out
               case ("BIN_OUT")
                  param%outfile = param_value
               case ("OUT_TYPE")
                  call swiftest_io_toupper(param_value)
                  param%out_type = param_value
               case ("OUT_FORM")
                  call swiftest_io_toupper(param_value)
                  param%out_form = param_value
               case ("OUT_STAT")
                  call swiftest_io_toupper(param_value)
                  param%out_stat = param_value
               case ("DUMP_CADENCE")
                  read(param_value, *, err = 667, iomsg = iomsg) param%dump_cadence
               case ("CHK_CLOSE")
                  call swiftest_io_toupper(param_value)
                  if (param_value == "YES" .or. param_value == 'T') param%lclose = .true.
               case ("CHK_RMIN")
                  read(param_value, *, err = 667, iomsg = iomsg) param%rmin
               case ("CHK_RMAX")
                  read(param_value, *, err = 667, iomsg = iomsg) param%rmax
               case ("CHK_EJECT")
                  read(param_value, *, err = 667, iomsg = iomsg) param%rmaxu
               case ("CHK_QMIN")
                  read(param_value, *, err = 667, iomsg = iomsg) param%qmin
               case ("CHK_QMIN_COORD")
                  call swiftest_io_toupper(param_value)
                  param%qmin_coord = param_value
               case ("CHK_QMIN_RANGE")
                  read(param_value, *, err = 667, iomsg = iomsg) param%qmin_alo
                  ifirst = ilast + 2
                  param_value = swiftest_io_get_token(line, ifirst, ilast, iostat)
                  read(param_value, *, err = 667, iomsg = iomsg) param%qmin_ahi
               case ("EXTRA_FORCE")
                  call swiftest_io_toupper(param_value)
                  if (param_value == "YES" .or. param_value == 'T') param%lextra_force = .true.
               case ("BIG_DISCARD")
                  call swiftest_io_toupper(param_value)
                  if (param_value == "YES" .or. param_value == 'T' ) param%lbig_discard = .true.
               case ("RHILL_PRESENT")
                  call swiftest_io_toupper(param_value)
                  if (param_value == "YES" .or. param_value == 'T' ) param%lrhill_present = .true.
               case ("MU2KG")
                  read(param_value, *, err = 667, iomsg = iomsg) param%MU2KG
               case ("TU2S")
                  read(param_value, *, err = 667, iomsg = iomsg) param%TU2S
               case ("DU2M")
                  read(param_value, *, err = 667, iomsg = iomsg) param%DU2M
               case ("ENERGY")
                  call swiftest_io_toupper(param_value)
                  if (param_value == "YES" .or. param_value == 'T') param%lenergy = .true.
               case ("GR")
                  call swiftest_io_toupper(param_value)
                  if (param_value == "YES" .or. param_value == 'T') param%lgr = .true. 
               case ("ROTATION")
                  call swiftest_io_toupper(param_value)
                  if (param_value == "YES" .or. param_value == 'T') param%lrotation = .true. 
               case ("TIDES")
                  call swiftest_io_toupper(param_value)
                  if (param_value == "YES" .or. param_value == 'T') param%ltides = .true. 
               case ("INTERACTION_LOOPS")
                  call swiftest_io_toupper(param_value)
                  param%interaction_loops = param_value
               case ("ENCOUNTER_CHECK_PLPL")
                  call swiftest_io_toupper(param_value)
                  param%encounter_check_plpl = param_value
               case ("ENCOUNTER_CHECK_PLTP")
                  call swiftest_io_toupper(param_value)
                  param%encounter_check_pltp = param_value
               case ("ENCOUNTER_CHECK")
                  call swiftest_io_toupper(param_value)
                  param%encounter_check_plpl = param_value
                  param%encounter_check_pltp = param_value
               case ("FIRSTKICK")
                  call swiftest_io_toupper(param_value)
                  if (param_value == "NO" .or. param_value == 'F') param%lfirstkick = .false. 
               case ("FIRSTENERGY")
                  call swiftest_io_toupper(param_value)
                  if (param_value == "NO" .or. param_value == 'F') param%lfirstenergy = .false. 
               case("EORBIT_ORIG")
                  read(param_value, *, err = 667, iomsg = iomsg) param%E_orbit_orig 
               case("GMTOT_ORIG")
                  read(param_value, *, err = 667, iomsg = iomsg) param%GMtot_orig 
               case("LTOT_ORIG")
                  read(param_value, *, err = 667, iomsg = iomsg) param%L_total_orig(1)
                  do i = 2, NDIM
                     ifirst = ilast + 2
                     param_value = swiftest_io_get_token(line, ifirst, ilast, iostat) 
                     read(param_value, *, err = 667, iomsg = iomsg) param%L_total_orig(i)
                  end do
               case("LORBIT_ORIG")
                  read(param_value, *, err = 667, iomsg = iomsg) param%L_orbit_orig(1)
                  do i = 2, NDIM
                     ifirst = ilast + 2
                     param_value = swiftest_io_get_token(line, ifirst, ilast, iostat) 
                     read(param_value, *, err = 667, iomsg = iomsg) param%L_orbit_orig(i)
                  end do
               case("LSPIN_ORIG")
                  read(param_value, *, err = 667, iomsg = iomsg) param%L_spin_orig(1)
                  do i = 2, NDIM
                     ifirst = ilast + 2
                     param_value = swiftest_io_get_token(line, ifirst, ilast, iostat) 
                     read(param_value, *, err = 667, iomsg = iomsg) param%L_spin_orig(i)
                  end do
               case("LESCAPE")
                  read(param_value, *, err = 667, iomsg = iomsg) param%L_escape(1)
                  do i = 2, NDIM
                     ifirst = ilast + 2
                     param_value = swiftest_io_get_token(line, ifirst, ilast, iostat) 
                     read(param_value, *, err = 667, iomsg = iomsg) param%L_escape(i)
                  end do
               case("GMESCAPE")
                  read(param_value, *, err = 667, iomsg = iomsg) param%GMescape 
               case("ECOLLISIONS")
                  read(param_value, *, err = 667, iomsg = iomsg) param%E_collisions
               case("EUNTRACKED")
                  read(param_value, *, err = 667, iomsg = iomsg) param%E_untracked
               case ("COLLISION_MODEL")
                  call swiftest_io_toupper(param_value)
                  read(param_value, *) param%collision_model
               case ("GMTINY")
                  read(param_value, *) param%GMTINY
               case ("MIN_GMFRAG")
                  read(param_value, *) param%min_GMfrag
               case ("NFRAG_REDUCTION")
                  read(param_value, *) param%nfrag_reduction
               case ("ENCOUNTER_SAVE")
                  call swiftest_io_toupper(param_value)
                  read(param_value, *) param%encounter_save
               case ("COARRAY")
                  call swiftest_io_toupper(param_value)
                  if (param_value == "YES" .or. param_value == 'T') param%lcoarray = .true. 
               case("SEED")
                  read(param_value, *) nseeds_from_file
                  ! Because the number of seeds can vary between compilers/systems, we need to make sure we can handle cases in 
                  ! which the input file has a different number of seeds than the current nbody_system. If the number of seeds in 
                  ! the file is smaller than required, we will use them as a source to fill in the missing elements.
                  ! If the number of seeds in the file is larger than required, we will truncate the seed array.
                  if (nseeds_from_file > nseeds) then
                     nseeds = nseeds_from_file
                     deallocate(param%seed)
                     allocate(param%seed(nseeds))
                     do i = 1, nseeds
                        ifirst = ilast + 2
                        param_value = swiftest_io_get_token(line, ifirst, ilast, iostat) 
                        read(param_value, *) param%seed(i)
                     end do
                  else ! Seed array in file is too small
                     do i = 1, nseeds_from_file
                        ifirst = ilast + 2
                        param_value = swiftest_io_get_token(line, ifirst, ilast, iostat) 
                        read(param_value, *) param%seed(i)
                     end do
                     param%seed(nseeds_from_file+1:nseeds) = [(param%seed(1) - param%seed(nseeds_from_file) + i, &
                                                               i=nseeds_from_file+1, nseeds)]
                  end if
                  seed_set = .true.
               case ("RESTART")
                  if (param_value == "NO" .or. param_value == 'F') then
                     param%lrestart = .false. 
                  else if (param_value == "YES" .or. param_value == 'T') then
                     param%lrestart = .true.
                  end if 
               ! Ignore SyMBA-specific, not-yet-implemented, or obsolete input parameters
               case ("NPLMAX", "NTPMAX", "YARKOVSKY", "YORP")
               case default
                  write(*,*) "Ignoring unknown parameter -> ",param_name
               end select
            end if
         end do
         1 continue
         close(unit)
         iostat = 0

         ! Do basic sanity checks on the input values
         if ((.not. tstart_set) .or. (.not. tstop_set) .or. (.not. dt_set)) then
            write(iomsg,*) 'Valid simulation time not set'
            iostat = -1
            return
         end if
         if (param%dt <= 0.0_DP) then
            write(iomsg,*) 'Invalid timestep: '
            iostat = -1
            return
         end if
         if (param%inplfile == "") then
            write(iomsg,*) 'No valid massive body file in input file'
            iostat = -1
            return
         end if
         if ((param%in_type /= "ASCII") .and. (param%in_type /= "NETCDF_FLOAT") .and. (param%in_type /= "NETCDF_DOUBLE"))  then
            write(iomsg,*) 'Invalid input file type:',trim(adjustl(param%in_type))
            iostat = -1
            return
         end if
         if (param%istep_out <= 0) then
            write(iomsg,*) 'Invalid ISTEP_OUT. Must be a positive integer'
            iostat = -1
            return
         end if
         if (param%nstep_out <= 0) then
            param%nstep_out = int((param%tstop - param%t0) / (param%istep_out * param%dt))
            param%fstep_out = 1.0_DP ! Linear output time
            param%ltstretch = .false.
         else
            param%fstep_out = 1._DP
            tratio = (param%TSTOP - param%T0) / (param%istep_out * param%dt)
            if (int(tratio) == param%nstep_out) then
               param%ltstretch = .false.
            else
               param%ltstretch = .true.
               y = time_stretcher(param%fstep_out) 
               call solve_roots(time_stretcher,param%fstep_out)
            end if
         end if
         if (param%dump_cadence < 0) then
            write(iomsg,*) 'Invalid DUMP_CADENCE. Must be a positive integer or 0.'
            iostat = -1
            return
         end if
         if ((param%istep_out > 0) .and. (param%outfile == "")) then
            write(iomsg,*) 'Invalid outfile'
            iostat = -1
            return
         end if
         param%lrestart = (param%out_stat == "APPEND")
         if (param%outfile /= "") then
            if ((param%out_type /= "NETCDF_FLOAT") .and. (param%out_type /= "NETCDF_DOUBLE")) then
               write(iomsg,*) 'Invalid out_type: ',trim(adjustl(param%out_type))
               iostat = -1
               return
            end if
            if ((param%out_form /= "EL") .and. (param%out_form /= "XV") .and. (param%out_form /= "XVEL")) then
               write(iomsg,*) 'Invalid out_form: ',trim(adjustl(param%out_form))
               iostat = -1
               return
            end if
            if ((param%out_stat /= "NEW") .and. (param%out_stat /= "REPLACE") .and. (param%out_stat /= "APPEND")  &
          .and. (param%out_stat /= "UNKNOWN")) then
               write(iomsg,*) 'Invalid out_stat: ',trim(adjustl(param%out_stat))
               iostat = -1
               return
            end if
         end if
         if (param%qmin > 0.0_DP) then
            if ((param%qmin_coord /= "HELIO") .and. (param%qmin_coord /= "BARY")) then
               write(iomsg,*) 'Invalid qmin_coord: ',trim(adjustl(param%qmin_coord))
               iostat = -1
               return
            end if
            if ((param%qmin_alo <= 0.0_DP) .or. (param%qmin_ahi <= 0.0_DP)) then
               write(iomsg,*) 'Invalid qmin vals'
               iostat = -1
               return
            end if
         end if
         if (param%ltides .and. .not. param%lrotation) then
            write(iomsg,*) 'Tides require rotation to be turned on'
            iostat = -1
            return
         end if

         if ((param%MU2KG < 0.0_DP) .or. (param%TU2S < 0.0_DP) .or. (param%DU2M < 0.0_DP)) then
            write(iomsg,*) 'Invalid unit conversion factor'
            iostat = -1
            return
         end if

         ! Calculate the G for the nbody_system units
         param%GU = GC / (param%DU2M**3 / (param%MU2KG * param%TU2S**2))


         if ((param%encounter_save /= "NONE")       .and. &
             (param%encounter_save /= "TRAJECTORY") .and. &
             (param%encounter_save /= "CLOSEST")    .and. &
             (param%encounter_save /= "BOTH")) then
            write(iomsg,*) 'Invalid encounter_save parameter: ',trim(adjustl(param%out_type))
            write(iomsg,*) 'Valid options are NONE, TRAJECTORY, CLOSEST, or BOTH'
            iostat = -1
            return
         end if

         param%lenc_save_trajectory = (param%encounter_save == "TRAJECTORY") .or. (param%encounter_save == "BOTH")
         param%lenc_save_closest = (param%encounter_save == "CLOSEST") .or. (param%encounter_save == "BOTH")

         if ((param%integrator == INT_RMVS) .or. (param%integrator == INT_SYMBA)) then
            if (.not.param%lclose) then
               write(iomsg,*) 'This integrator requires CHK_CLOSE to be enabled.'
               iostat = -1
               return
            end if
         end if

         param%lmtiny_pl = (param%integrator == INT_SYMBA) 

         if (param%lmtiny_pl .and. param%GMTINY < 0.0_DP) then
            write(iomsg,*) "GMTINY invalid or not set: ", param%GMTINY
            iostat = -1
            return
         end if

         if ((param%collision_model /= "MERGE")       .and. &
             (param%collision_model /= "BOUNCE")    .and. &
             (param%collision_model /= "FRAGGLE")) then
            write(iomsg,*) 'Invalid collision_model parameter: ',trim(adjustl(param%out_type))
            write(iomsg,*) 'Valid options are MERGE, BOUNCE, or FRAGGLE'
            iostat = -1
            return
         end if

         if (param%collision_model == "FRAGGLE" ) then
            if (seed_set) then
               call random_seed(put = param%seed)
            else
               call random_seed(get = param%seed)
            end if
            if (param%min_GMfrag < 0.0_DP) param%min_GMfrag = param%GMTINY
            if (param%nfrag_reduction < 1.0_DP) then
               write(iomsg,*) "Warning: NFRAG_REDUCTION value invalid. Setting to 1.0" 
               param%nfrag_reduction = 1.0_DP
            end if
         end if
   
         ! Determine if the GR flag is set correctly for this integrator
         select case(param%integrator)
         case(INT_WHM, INT_RMVS, INT_HELIO, INT_SYMBA)
         case default   
            if (param%lgr) write(iomsg, *) 'GR is not yet implemented for this integrator. This parameter will be ignored.'
            param%lgr = .false.
         end select

         if (param%lgr) then
            ! Calculate the inverse speed of light in the nbody_system units
            param%inv_c2 = einsteinC * param%TU2S / param%DU2M
            param%inv_c2 = (param%inv_c2)**(-2)
         end if

         select case(trim(adjustl(param%interaction_loops)))
         case("TRIANGULAR")
            param%lflatten_interactions = .false.
         case("FLAT")
            param%lflatten_interactions = .true.
         case default
            write(*,*) "Unknown value for parameter INTERACTION_LOOPS: -> ",trim(adjustl(param%interaction_loops))
            write(*,*) "Must be one of the following: TRIANGULAR or FLAT"
            write(*,*) "Using default value of TRIANGULAR"
            param%interaction_loops = "TRIANGULAR"
            param%lflatten_interactions = .false.
         end select

         select case(trim(adjustl(param%encounter_check_plpl)))
         case("TRIANGULAR")
            param%lencounter_sas_plpl = .false.
         case("SORTSWEEP")
            param%lencounter_sas_plpl = .true.
         case default
            write(*,*) "Unknown value for parameter ENCOUNTER_CHECK_PLPL: -> ",trim(adjustl(param%encounter_check_plpl))
            write(*,*) "Must be one of the following: TRIANGULAR or SORTSWEEP"
            write(*,*) "Using default value of TRIANGULAR"
            param%encounter_check_plpl = "TRIANGULAR"
            param%lencounter_sas_plpl = .false.
         end select

         select case(trim(adjustl(param%encounter_check_pltp)))
         case("TRIANGULAR")
            param%lencounter_sas_pltp = .false.
         case("SORTSWEEP")
            param%lencounter_sas_pltp = .true.
         case default
            write(*,*) "Unknown value for parameter ENCOUNTER_CHECK_PLTP: -> ",trim(adjustl(param%encounter_check_pltp))
            write(*,*) "Must be one of the following: TRIANGULAR or SORTSWEEP"
            write(*,*) "Using default value of TRIANGULAR"
            param%encounter_check_pltp = "TRIANGULAR"
            param%lencounter_sas_pltp = .false.
         end select


         if (param%lcoarray) then
#ifdef COARRAY
            if (num_images() == 1) then
               write(iomsg, *) "Only one Coarray image detected. Coarrays will not be used."
               param%lcoarray = .false.
            end if

            select case(param%integrator)
            case(INT_WHM, INT_RMVS, INT_HELIO)
            case default   
               write(iomsg, *) "Coarray-based parallelization of test particles are not compatible with this integrator. " &
                            // "This parameter will be ignored."
               param%lcoarray = .false.
            end select
#else
            write(iomsg,*) "Coarray capability not detected. Swiftest must be compiled with Coarrays enabled. to use this feature."
            param%lcoarray = .false.
#endif
         end if

         iostat = 0

      end associate

#ifdef COARRAY
   end if ! this_image() == 1
      call coparam%coclone()
#endif
      select type(param => self)
      type is (swiftest_parameters)
#ifdef COARRAY
         param = coparam
#endif
         call param%set_display(param%display_style)

         if (.not.param%lrestart) then
#ifdef COARRAY
            if (this_image() == 1 .or. param%log_output) then
#endif
               call param%writer(unit = param%display_unit, iotype = "none", v_list = [0], iostat = iostat, iomsg = iomsg)
               if (param%log_output) flush(param%display_unit) 
#ifdef COARRAY
            end if !(this_image() == 1)
            write(COLLISION_LOG_OUT,'("collision_coimage",I0.3,".log")') this_image()
#endif
            ! A minimal log of collision outcomes is stored in the following log file
            ! More complete data on collisions is stored in the NetCDF output files
            call swiftest_io_log_start(param, COLLISION_LOG_OUT, "Collision logfile")
         end if
         ! Print the contents of the parameter file to standard output
      end select

      return 
      667 continue
      write(*,*) "Error reading param file: ", trim(adjustl(iomsg))

      contains
         function time_stretcher(fstep_out) result(ans)
            !! author: David A. Minton
            !!
            !! Equation for the time stretchinf function. Solving the roots of this equation gives the time stretching factor for 
            !! non-linear file output cadence.
            implicit none
            ! Arguments
            real(DP), intent(in) :: fstep_out
            ! Result
            real(DP)             :: ans

            if (abs(fstep_out-1.0_DP) < epsilon(1.0_DP)) then
               ans = self%nstep_out - tratio
            else
               ans = (1.0_DP - fstep_out**(self%nstep_out))/ (1.0_DP - fstep_out) - tratio
            end if

            return
         end function time_stretcher

   end subroutine swiftest_io_param_reader


   module subroutine swiftest_io_param_writer(self, unit, iotype, v_list, iostat, iomsg) 
      !! author: David A. Minton
      !!
      !! Dump integration parameters to file
      !!
      !! Adapted from David E. Kaufmann's Swifter routine io_param_restart.f90
      !! Adapted from Martin Duncan's Swift routine io_param_restart.f
      implicit none
      ! Arguments
      class(swiftest_parameters),intent(in)     :: self       !! Collection of parameters
      integer, intent(in)                       :: unit       !! File unit number
      character(len=*), intent(in)              :: iotype     !! Dummy argument passed to the  input/output procedure contains the 
                                                              !! text from the char-literal-constant, prefixed with DT. 
                                                              !!    If you do not include a char-literal-constant, the iotype 
                                                              !!    argument contains only DT.
      integer, intent(in)                       :: v_list(:)  !! Not used in this procedure
      integer, intent(out)                      :: iostat     !! IO status code
      character(len=*), intent(inout)           :: iomsg      !! Message to pass if iostat /= 0
      ! Internals
      integer(I4B)                   :: nseeds

      associate(param => self)
         call io_param_writer_one("T0", param%t0, unit)
         call io_param_writer_one("TSTART", param%tstart, unit)
         call io_param_writer_one("TSTOP", param%tstop, unit)
         call io_param_writer_one("DT", param%dt, unit)
         call io_param_writer_one("IN_TYPE", param%in_type, unit)
         if (param%in_type == "ASCII") then
            call io_param_writer_one("CB_IN", param%incbfile, unit)
            call io_param_writer_one("PL_IN", param%inplfile, unit)
            call io_param_writer_one("TP_IN", param%intpfile, unit)
         else 
            call io_param_writer_one("NC_IN", param%nc_in, unit)
         end if

         call io_param_writer_one("IN_FORM", param%in_form, unit)
         if (param%dump_cadence > 0) call io_param_writer_one("DUMP_CADENCE",param%dump_cadence, unit)
         if (param%istep_out > 0) then
            call io_param_writer_one("ISTEP_OUT", param%istep_out, unit)
            call io_param_writer_one("BIN_OUT", param%outfile, unit)
            call io_param_writer_one("OUT_TYPE", param%out_type, unit)
            call io_param_writer_one("OUT_FORM", param%out_form, unit)
            call io_param_writer_one("OUT_STAT", "APPEND", unit) 
         end if
         call io_param_writer_one("CHK_RMIN", param%rmin, unit)
         call io_param_writer_one("CHK_RMAX", param%rmax, unit)
         call io_param_writer_one("CHK_EJECT", param%rmaxu, unit)
         call io_param_writer_one("CHK_QMIN", param%qmin, unit)
         if (param%qmin >= 0.0_DP) then
            call io_param_writer_one("CHK_QMIN_COORD", param%qmin_coord, unit)
            call io_param_writer_one("CHK_QMIN_RANGE", [param%qmin_alo, param%qmin_ahi], unit)
         end if
         call io_param_writer_one("MU2KG", param%MU2KG, unit)
         call io_param_writer_one("TU2S", param%TU2S , unit)
         call io_param_writer_one("DU2M", param%DU2M, unit)
         call io_param_writer_one("RHILL_PRESENT", param%lrhill_present, unit)
         call io_param_writer_one("EXTRA_FORCE", param%lextra_force, unit)
         call io_param_writer_one("CHK_CLOSE", param%lclose, unit)
         call io_param_writer_one("ENERGY", param%lenergy, unit)
         call io_param_writer_one("GR", param%lgr, unit)
         call io_param_writer_one("ROTATION", param%lrotation, unit)
         call io_param_writer_one("TIDES", param%ltides, unit)
         call io_param_writer_one("INTERACTION_LOOPS", param%interaction_loops, unit)
         call io_param_writer_one("ENCOUNTER_CHECK_PLPL", param%encounter_check_plpl, unit)
         call io_param_writer_one("ENCOUNTER_CHECK_PLTP", param%encounter_check_pltp, unit)
         call io_param_writer_one("ENCOUNTER_SAVE", param%encounter_save, unit)
         call io_param_writer_one("COARRAY", param%lcoarray, unit)

         if (param%lenergy) then
            call io_param_writer_one("FIRSTENERGY", param%lfirstenergy, unit)
         end if
         call io_param_writer_one("FIRSTKICK",param%lfirstkick, unit)

         if (param%GMTINY >= 0.0_DP) call io_param_writer_one("GMTINY",param%GMTINY, unit)
         if (param%min_GMfrag >= 0.0_DP) call io_param_writer_one("MIN_GMFRAG",param%min_GMfrag, unit)
         call io_param_writer_one("COLLISION_MODEL",param%collision_model, unit)
         if (param%collision_model == "FRAGGLE" ) then
            call io_param_writer_one("NFRAG_REDUCTION",param%nfrag_reduction, unit)
            nseeds = size(param%seed)
            call random_seed(get = param%seed)
            call io_param_writer_one("SEED", [nseeds, param%seed(:)], unit)
         end if
   
         iostat = 0
         iomsg = "UDIO not implemented"
      end associate

      return
   end subroutine swiftest_io_param_writer


   module subroutine swiftest_io_param_writer_one_char(param_name, param_value, unit)
      !! author: David A. Minton
      !!
      !! Writes a single parameter name/value pair to a file unit. 
      !! This version is for character param_value type
      implicit none
      ! Arguments
      character(len=*), intent(in) :: param_name  !! Name of parameter to print
      character(len=*), intent(in) :: param_value !! Value of parameter to print
      integer(I4B),     intent(in) :: unit        !! Open file unit number to print parameter to
      ! Internals
      character(len=NAMELEN) :: param_name_fixed_width !! Parameter label converted to fixed-width string
      character(len=STRMAX)  :: iomsg             !! Message to pass if iostat /= 0

      write(param_name_fixed_width, *) param_name
      write(unit, *, err = 667, iomsg = iomsg) adjustl(param_name_fixed_width) // " " // trim(adjustl(param_value)) 

      return
      667 continue
      write(*,*) 'Error writing parameter: ',trim(adjustl(iomsg))
   end subroutine swiftest_io_param_writer_one_char


   module subroutine swiftest_io_param_writer_one_DP(param_name, param_value, unit)
      !! author: David A. Minton
      !!
      !! Writes a single parameter name/value pair to a file unit. 
      !! This version is for real(DP) param_value type
      implicit none
      ! Arguments
      character(len=*), intent(in)    :: param_name  !! Name of parameter to print
      real(DP),         intent(in)    :: param_value !! Value of parameter to print
      integer(I4B),     intent(in)    :: unit        !! Open file unit number to print parameter to
      ! Internals
      character(len=STRMAX) :: param_value_string   !! Parameter value converted to a string
      character(*),parameter :: Rfmt  = '(ES25.17)' !! Format label for real values 

      write(param_value_string,Rfmt) param_value
      call io_param_writer_one(param_name, param_value_string, unit)

      return
   end subroutine swiftest_io_param_writer_one_DP


   module subroutine swiftest_io_param_writer_one_DParr(param_name, param_value, unit)
      !! author: David A. Minton
      !!
      !! Writes a single parameter name/value pair to a file unit. 
      !! This version is for real(DP) arrays () param_value type
      implicit none
      ! Arguments
      character(len=*),       intent(in) :: param_name  !! Name of parameter to print
      real(DP), dimension(:), intent(in) :: param_value !! Value of parameter to print
      integer(I4B),           intent(in) :: unit        !! Open file unit number to print parameter to
      ! Internals
      character(len=STRMAX) :: param_value_string   !! Parameter value converted to a string
      character(*),parameter :: Rfmt  = '(ES25.17)' !! Format label for real values 
      character(len=25) :: arr_val
      integer(I4B) :: i, narr

      narr = size(param_value)
      do i = 1, narr
         write(arr_val, Rfmt) param_value(i)
         if (i == 1) then
            write(param_value_string, *) arr_val
         else
            param_value_string = trim(adjustl(param_value_string)) // " " // arr_val
         end if
      end do

      call io_param_writer_one(param_name, param_value_string, unit)

      return
   end subroutine swiftest_io_param_writer_one_DParr


   module subroutine swiftest_io_param_writer_one_I4B(param_name, param_value, unit)
      !! author: David A. Minton
      !!
      !! Writes a single parameter name/value pair to a file unit. 
      !! This version is for integer(I4B) param_value type
      implicit none
      ! Arguments
      character(len=*), intent(in) :: param_name  !! Name of parameter to print
      integer(I4B),     intent(in) :: param_value !! Value of parameter to print
      integer(I4B),     intent(in) :: unit        !! Open file unit number to print parameter to
      ! Internals
      character(len=STRMAX) :: param_value_string   !! Parameter value converted to a string
      character(*),parameter :: Ifmt  = '(I0)'      !! Format label for integer values

      write(param_value_string,Ifmt) param_value
      call io_param_writer_one(param_name, param_value_string, unit)

      return
   end subroutine swiftest_io_param_writer_one_I4B


   module subroutine swiftest_io_param_writer_one_I8B(param_name, param_value, unit)
      !! author: David A. Minton
      !!
      !! Writes a single parameter name/value pair to a file unit. 
      !! This version is for integer(I8B) param_value type
      implicit none
      ! Arguments
      character(len=*), intent(in)    :: param_name  !! Name of parameter to print
      integer(I8B),     intent(in)    :: param_value !! Value of parameter to print
      integer(I4B),     intent(in)    :: unit        !! Open file unit number to print parameter to
      ! Internals
      character(len=STRMAX) :: param_value_string   !! Parameter value converted to a string
      character(*),parameter :: Ifmt  = '(I0)'      !! Format label for integer values

      write(param_value_string,Ifmt) param_value
      call io_param_writer_one(param_name, param_value_string, unit)

      return
   end subroutine swiftest_io_param_writer_one_I8B


   module subroutine swiftest_io_param_writer_one_I4Barr(param_name, param_value, unit)
      !! author: David A. Minton
      !!
      !! Writes a single parameter name/value pair to a file unit. 
      !! This version is for integer(I4B) arrays param_value type
      implicit none
      ! Arguments
      character(len=*),           intent(in) :: param_name  !! Name of parameter to print
      integer(I4B), dimension(:), intent(in) :: param_value !! Value of parameter to print
      integer(I4B),               intent(in) :: unit        !! Open file unit number to print parameter to
      ! Internals
      character(len=STRMAX) :: param_value_string !! Parameter value converted to a string
      character(*),parameter :: Ifmt  = '(I0)'    !! Format label for integer values
      character(len=25) :: arr_val
      integer(I4B) :: i, narr

      narr = size(param_value)
      do i = 1, narr
         write(arr_val, Ifmt) param_value(i)
         if (i == 1) then
            write(param_value_string, *) trim(adjustl(arr_val))
         else
            param_value_string = trim(adjustl(param_value_string)) // " " // trim(adjustl(arr_val))
         end if
      end do

      call io_param_writer_one(param_name, param_value_string, unit)

      return
   end subroutine swiftest_io_param_writer_one_I4Barr


   module subroutine swiftest_io_param_writer_one_logical(param_name, param_value, unit)
      !! author: David A. Minton
      !!
      !! Writes a single parameter name/value pair to a file unit. 
      !! This version is for logical param_value type
      implicit none
      ! Arguments
      character(len=*), intent(in) :: param_name  !! Name of parameter to print
      logical,          intent(in) :: param_value !! Value of parameter to print
      integer(I4B),     intent(in) :: unit        !! Open file unit number to print parameter to
      ! Internals
      character(len=STRMAX) :: param_value_string   !! Parameter value converted to a string
      character(*),parameter :: Lfmt  = '(L1)'         !! Format label for logical values 

      write(param_value_string,Lfmt) param_value
      call io_param_writer_one(param_name, param_value_string, unit)

      return
   end subroutine swiftest_io_param_writer_one_logical

#ifdef QUADPREC
   module subroutine swiftest_io_param_writer_one_QP(param_name, param_value, unit)
      !! author: David A. Minton
      !!
      !! Writes a single parameter name/value pair to a file unit. 
      !! This version is for real(QP) param_value type
      implicit none
      ! Arguments
      character(len=*), intent(in) :: param_name  !! Name of parameter to print
      real(QP),         intent(in) :: param_value !! Value of parameter to print
      integer(I4B),     intent(in) :: unit        !! Open file unit number to print parameter to
      ! Internals
      character(len=STRMAX) :: param_value_string   !! Parameter value converted to a string
      character(*),parameter :: Rfmt  = '(ES25.17)' !! Format label for real values 

      write(param_value_string,Rfmt) param_value
      call io_param_writer_one(param_name, param_value_string, unit)

      return
   end subroutine swiftest_io_param_writer_one_QP
#endif

   module subroutine swiftest_io_read_in_body(self, param) 
      !! author: The Purdue Swiftest Team - David A. Minton, Carlisle A. Wishard, Jennifer L.L. Pouplin, and Jacob R. Elliott
      !!
      !! Read in either test particle or massive body data 
      !!
      !! Adapted from David E. Kaufmann's Swifter routine swiftest_init_pl.f90 and swiftest_init_tp.f90
      !! Adapted from Martin Duncan's Swift routine swiftest_init_pl.f and swiftest_init_tp.f
      implicit none
      ! Arguments
      class(swiftest_body),       intent(inout) :: self   !! Swiftest particle object
      class(swiftest_parameters), intent(inout) :: param !! Current run configuration parameters
      ! Internals
      integer(I4B)                  :: iu = LUN
      integer(I4B)                  :: i, nbody
      character(len=:), allocatable :: infile
      character(STRMAX)             :: errmsg
      ! Internals
      integer(I4B)                                :: ierr  !! Error code: returns 0 if the read is successful

      ! Select the appropriate polymorphic class (test particle or massive body)
      if (param%in_type /= "ASCII") return ! Not for NetCDF

      select type(self)
      class is (swiftest_pl)
         infile = param%inplfile
      class is (swiftest_tp)
         infile = param%intpfile
      end select

      open(unit = iu, file = infile, status = 'old', form = 'FORMATTED', err = 667, iomsg = errmsg)
      read(iu, *, err = 667, iomsg = errmsg) nbody

      call self%setup(nbody, param)
      ierr = 0
      if (nbody > 0) then
         ierr = self%read_frame(iu, param)
         self%status(:) = ACTIVE
         self%lmask(:) = .true.
         do i = 1, nbody
            call self%info(i)%set_value(status="ACTIVE")
         end do
      end if
      close(iu, err = 667, iomsg = errmsg)

      if (ierr == 0) return

      667 continue
      write(*,*) 'Error reading in initial conditions file: ',trim(adjustl(errmsg))
      return
   end subroutine swiftest_io_read_in_body


   module subroutine swiftest_io_read_in_cb(self, param) 
      !! author: David A. Minton
      !!
      !! Reads in central body data 
      !!
      !! Adapted from David E. Kaufmann's Swifter routine swiftest_init_pl.f90
      !! Adapted from Martin Duncan's Swift routine swiftest_init_pl.f
      implicit none
      ! Arguments
      class(swiftest_cb),         intent(inout) :: self
      class(swiftest_parameters), intent(inout) :: param
      ! Internals
      integer(I4B)            :: iu = LUN
      character(len=STRMAX)   :: errmsg
      integer(I4B)            :: ierr
      character(len=NAMELEN)  :: name

      if (param%in_type /= "ASCII") return ! Not for NetCDF

      self%id = 0
      open(unit = iu, file = param%incbfile, status = 'old', form = 'FORMATTED', err = 667, iomsg = errmsg)
      read(iu, *, err = 667, iomsg = errmsg) name
      call self%info%set_value(name=name)
      read(iu, *, err = 667, iomsg = errmsg) self%Gmass
      self%mass = real(self%Gmass / param%GU, kind=DP)
      read(iu, *, err = 667, iomsg = errmsg) self%radius
      read(iu, *, err = 667, iomsg = errmsg) self%j2rp2
      read(iu, *, err = 667, iomsg = errmsg) self%j4rp4
      if (param%lrotation) then
         read(iu, *, err = 667, iomsg = errmsg) self%Ip(1), self%Ip(2), self%Ip(3)
         read(iu, *, err = 667, iomsg = errmsg) self%rot(1), self%rot(2), self%rot(3)
         self%rot(:) = self%rot(:) * DEG2RAD
      end if
      ierr = 0
      close(iu, err = 667, iomsg = errmsg)

      if (ierr == 0) then
   
         if (param%rmin < 0.0) param%rmin = self%radius
         
         self%GM0 = self%Gmass
         self%dGM = 0.0_DP
         self%R0 = self%radius
         if (param%lrotation) then
            self%L0(:) = self%Ip(3) * self%mass * self%radius**2 * self%rot(:)
            self%dL(:) = 0.0_DP
         end if
      end if
      return

      667 continue
      write(*,*) "Error reading central body file: " // trim(adjustl(errmsg))
      call base_util_exit(FAILURE,param%display_unit)
   end subroutine swiftest_io_read_in_cb


   module subroutine swiftest_io_read_in_system(self, nc, param)
      !! author: David A. Minton and Carlisle A. Wishard
      !!
      !! Reads in the nbody_system from input files
      implicit none
      ! Arguments
      class(swiftest_nbody_system),      intent(inout) :: self
      class(swiftest_netcdf_parameters), intent(inout) :: nc     !! Parameters used to identify a particular NetCDF dataset
      class(swiftest_parameters),        intent(inout) :: param
      ! Internals
      integer(I4B) :: ierr, i
      class(swiftest_parameters), allocatable :: tmp_param

      if (param%in_type == "ASCII") then
         call self%cb%read_in(param)
         call self%pl%read_in(param)
         if (self%pl%nbody > 0) self%pl%id(:) = [(i, i = 1, self%pl%nbody)]
         call self%tp%read_in(param)
         if (self%tp%nbody > 0) self%tp%id(:) = [(i, i = self%pl%nbody + 1, self%pl%nbody + 1 + self%tp%nbody)]
         self%maxid = self%pl%nbody + self%tp%nbody
         ! Copy over param file variable inputs
         self%E_orbit_orig = param%E_orbit_orig
         self%GMtot_orig = param%GMtot_orig
         self%L_total_orig(:) = param%L_total_orig(:)
         self%L_orbit_orig(:) = param%L_orbit_orig(:)
         self%L_spin_orig(:) = param%L_spin_orig(:)
         self%L_escape(:) = param%L_escape(:)
         self%E_collisions = param%E_collisions
         self%E_untracked = param%E_untracked
      else
         allocate(tmp_param, source=param)
         nc%file_name = param%nc_in
         tmp_param%out_form = param%in_form
         if (.not. param%lrestart) then
            ! Turn off energy computation so we don't have to feed it into the initial conditions
            tmp_param%lenergy = .false.
         end if
         ierr = self%read_frame(nc, tmp_param)
         deallocate(tmp_param)
         if (ierr /=0) call base_util_exit(FAILURE,param%display_unit)
      end if

      param%loblatecb = ((abs(self%cb%j2rp2) > 0.0_DP) .or. (abs(self%cb%j4rp4) > 0.0_DP))
      if (.not.param%loblatecb) then
         if (allocated(self%pl%aobl)) deallocate(self%pl%aobl)
         if (allocated(self%tp%aobl)) deallocate(self%tp%aobl)
      else
         if (self%pl%nbody > 0) then
            if (.not. allocated(self%pl%aobl)) allocate(self%pl%aobl(NDIM,self%pl%nbody))
            self%pl%aobl(:,:) = 0.0_DP
         end if
         if (self%tp%nbody > 0) then
            if (.not. allocated(self%tp%aobl)) allocate(self%tp%aobl(NDIM,self%tp%nbody))
            self%tp%aobl(:,:) = 0.0_DP
         end if
      end if

      return
   end subroutine swiftest_io_read_in_system


   module function swiftest_io_read_frame_body(self, iu, param) result(ierr)
      !! author: David A. Minton
      !!
      !! Reads a frame of output of either test particle or massive body data from a binary output file
      !!
      !! Adapted from David E. Kaufmann's Swifter routine  io_read_frame.f90
      !! Adapted from Hal Levison's Swift routine io_read_frame.f
      implicit none
      ! Arguments
      class(swiftest_body),       intent(inout) :: self   !! Swiftest particle object
      integer(I4B),               intent(inout) :: iu     !! Unit number for the output file to write frame to
      class(swiftest_parameters), intent(inout) :: param  !! Current run configuration parameters 
      ! Result
      integer(I4B)                              :: ierr  !! Error code: returns 0 if the read is successful
      ! Internals
      character(len=STRMAX)   :: errmsg
      character(len=NAMELEN), dimension(self%nbody)  :: name
      integer(I4B) :: i
      real(QP)                      :: val

      if (self%nbody == 0) return

      if ((param%in_form /= "EL") .and. (param%in_form /= "XV")) then
         write(errmsg, *) trim(adjustl(param%in_form)) // " is not a recognized format code for input files."
         goto 667
      end if

      associate(n => self%nbody)

         if (param%in_form == "EL") then
            if (.not.allocated(self%a))     allocate(self%a(n))
            if (.not.allocated(self%e))     allocate(self%e(n))
            if (.not.allocated(self%inc))   allocate(self%inc(n))
            if (.not.allocated(self%capom)) allocate(self%capom(n))
            if (.not.allocated(self%omega)) allocate(self%omega(n))
            if (.not.allocated(self%capm))  allocate(self%capm(n))
         end if

         select case(param%in_type)
         case ("ASCII")
            do i = 1, n
               select type(self)
               class is (swiftest_pl)
                  if (param%lrhill_present) then
                     read(iu, *, err = 667, iomsg = errmsg) name(i), val, self%rhill(i)
                  else
                     read(iu, *, err = 667, iomsg = errmsg) name(i), val
                  end if
                  self%Gmass(i) = real(val, kind=DP)
                  self%mass(i) = real(val / param%GU, kind=DP)
                  if (param%lclose) read(iu, *, err = 667, iomsg = errmsg) self%radius(i)
               class is (swiftest_tp)
                  read(iu, *, err = 667, iomsg = errmsg) name(i)
               end select
               call self%info(i)%set_value(name=name(i))

               select case(param%in_form)
               case ("XV")
                  read(iu, *, err = 667, iomsg = errmsg) self%rh(1, i), self%rh(2, i), self%rh(3, i)
                  read(iu, *, err = 667, iomsg = errmsg) self%vh(1, i), self%vh(2, i), self%vh(3, i)
               case ("EL")
                  read(iu, *, err = 667, iomsg = errmsg) self%a(i), self%e(i), self%inc(i)
                  read(iu, *, err = 667, iomsg = errmsg) self%capom(i), self%omega(i), self%capm(i)
               end select

               select type (self)
               class is (swiftest_pl)
                  if (param%lrotation) then
                     read(iu, *, err = 667, iomsg = errmsg) self%Ip(1, i), self%Ip(2, i), self%Ip(3, i)
                     read(iu, *, err = 667, iomsg = errmsg) self%rot(1, i), self%rot(2, i), self%rot(3, i)
                     self%rot(:,i) = self%rot(:,i) * DEG2RAD 
                  end if
                  ! if (param%ltides) then
                  !    read(iu, *, err = 667, iomsg = errmsg) self%k2(i)
                  !    read(iu, *, err = 667, iomsg = errmsg) self%Q(i)
                  ! end if
               end select
            end do
         end select

         if (param%in_form == "EL") then
            self%inc(1:n)   = self%inc(1:n) * DEG2RAD
            self%capom(1:n) = self%capom(1:n) * DEG2RAD
            self%omega(1:n) = self%omega(1:n) * DEG2RAD
            self%capm(1:n)  = self%capm(1:n) * DEG2RAD
         end if
      end associate

      ierr = 0
      return

      667 continue
      select type (self)
      class is (swiftest_pl)
         write(*,*) "Error reading massive body file: " // trim(adjustl(errmsg))
      class is (swiftest_tp)
         write(*,*) "Error reading test particle file: " // trim(adjustl(errmsg))
      class default
         write(*,*) "Error reading body file: " // trim(adjustl(errmsg))
      end select
      call base_util_exit(FAILURE,param%display_unit)
   end function swiftest_io_read_frame_body


   module subroutine swiftest_io_read_in_param(self, param_file_name) 
      !! author: David A. Minton
      !!
      !! Read in parameters for the integration
      !!
      !! Adapted from David E. Kaufmann's Swifter routine io_init_param.f90
      !! Adapted from Martin Duncan's Swift routine io_init_param.f
      implicit none
      ! Arguments
      class(swiftest_parameters),intent(inout) :: self             !! Current run configuration parameters
      character(len=*), intent(in)             :: param_file_name !! Parameter input file name (i.e. param.in)
      ! Internals
      integer(I4B)      :: ierr = 0 !! Input error code
      character(STRMAX) :: errmsg   !! Error message in UDIO procedure

      ! Read in name of parameter file
      self%param_file_name = trim(adjustl(param_file_name))

      !! todo: Currently this procedure does not work in user-defined derived-type input mode 
      !!    as the newline characters are ignored in the input file when compiled in ifort.

      !read(LUN,'(DT)', iostat= ierr, iomsg = errmsg) self
      call self%reader(LUN, iotype= "none", v_list=[""], iostat = ierr, iomsg = errmsg)
      if (ierr == 0) return

      write(self%display_unit,*) "Error reading parameter file: " // trim(adjustl(errmsg))
      call base_util_exit(FAILURE)
   end subroutine swiftest_io_read_in_param


   module subroutine swiftest_io_set_display_param(self, display_style)
      !! author: David A. Minton
      !!
      !! Sets the display style parameters. If display is "STANDARD" then output goes to stdout. If display is "COMPACT" 
      !! then it is redirected to a log file and a progress-bar is used for stdout
      implicit none
      ! Arguments
      class(swiftest_parameters), intent(inout) :: self            !! Current run configuration parameters
      character(*),               intent(in)    :: display_style   !! Style of the output display 
      ! Internals
      character(STRMAX) :: errmsg
      logical           :: fileExists   

      select case(display_style)
      case ('STANDARD')
         self%display_unit = OUTPUT_UNIT !! stdout from iso_fortran_env
         self%log_output = .false.
      case ('COMPACT', 'PROGRESS')
#ifdef COARRAY
         if (self%lcoarray) then
            write(SWIFTEST_LOG_FILE,'("swiftest_coimage",I0.3,".log")') this_image()
         else
            write(SWIFTEST_LOG_FILE,'("swiftest.log")')
         end if 
#endif
         inquire(file=SWIFTEST_LOG_FILE, exist=fileExists)
         if (self%lrestart.and.fileExists) then
            open(unit=SWIFTEST_LOG_OUT, file=SWIFTEST_LOG_FILE, status="OLD", position="APPEND", err = 667, iomsg = errmsg)
         else
            open(unit=SWIFTEST_LOG_OUT, file=SWIFTEST_LOG_FILE, status="REPLACE", err = 667, iomsg = errmsg)
         end if
         self%display_unit = SWIFTEST_LOG_OUT 
         self%log_output = .true.
      case default
         write(*,*) display_style, " is an unknown display style"
         call base_util_exit(USAGE)
      end select

      self%display_style = display_style

      return

      667 continue
      write(*,*) "Error opening swiftest log file: " // trim(adjustl(errmsg))
      call base_util_exit(FAILURE,self%display_unit)
   end subroutine swiftest_io_set_display_param


   module subroutine swiftest_io_toupper(string)
      !! author: David A. Minton
      !!
      !! Convert string to uppercase
      !!
      !! Adapted from David E. Kaufmann's Swifter routine: util_toupper.f90
      implicit none
      ! Arguments
      character(*), intent(inout) :: string !! String to make upper case
      ! Internals
      integer(I4B) :: i, length, idx
   
      length = len(string)
      do i = 1, length
         idx = iachar(string(i:i))
         if ((idx >= lowercase_begin) .and. (idx <= lowercase_end)) then
            idx = idx + uppercase_offset
            string(i:i) = achar(idx)
         end if
      end do
   
      return
   end subroutine swiftest_io_toupper


   module subroutine swiftest_io_initialize_output_file_system(self, nc, param)
      !! author: The Purdue Swiftest Team - David A. Minton, Carlisle A. Wishard, Jennifer L.L. Pouplin, and Jacob R. Elliott
      !!
      !! Write a frame (header plus records for each massive body and active test particle) to output binary file
      !! There is no direct file output from this subroutine
      !!
      !! Adapted from David E. Kaufmann's Swifter routine  io_write_frame.f90
      !! Adapted from Hal Levison's Swift routine io_write_frame.f
      implicit none
      ! Arguments
      class(swiftest_nbody_system),      intent(inout) :: self   !! Swiftest nbody_system object
      class(swiftest_netcdf_parameters), intent(inout) :: nc     !! Parameters used to identify a particular NetCDF dataset
      class(swiftest_parameters),        intent(inout) :: param !! Current run configuration parameters 
      ! Internals

      character(len=2*STRMAX)          :: errmsg
      logical                          :: fileExists

      associate (pl => self%pl, tp => self%tp, npl => self%pl%nbody, ntp => self%tp%nbody, lfirst => self%lfirst_io)
         nc%file_name = param%outfile
         if (lfirst) then
            inquire(file=param%outfile, exist=fileExists)
#ifdef COARRAY
            if (this_image() /= 1) param%out_stat = 'APPEND'
#endif
            
            select case(param%out_stat)
            case('APPEND')
               if (.not.fileExists) then
                  errmsg = trim(adjustl(param%outfile)) // " not found! You must specify OUT_STAT = NEW, REPLACE, or UNKNOWN"
                  goto 667
               end if
               call nc%open(param)
            case('NEW')
               if (fileExists) then
                  errmsg = trim(adjustl(param%outfile))// " already exists! You must specify OUT_STAT = APPEND, REPLACE, or UNKNOWN"
                  goto 667
               end if
               call nc%initialize(param)
            case('REPLACE', 'UNKNOWN')
               call nc%initialize(param)
            end select

            lfirst = .false.
         end if

      end associate

      return

      667 continue
      write(*,*) "Error writing nbody_system frame: " // trim(adjustl(errmsg))
      call base_util_exit(FAILURE,param%display_unit)
   end subroutine swiftest_io_initialize_output_file_system

end submodule s_swiftest_io
