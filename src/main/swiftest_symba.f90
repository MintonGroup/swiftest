program swiftest_symba
   !! author: The Purdue Swiftest Team -  David A. Minton, Carlisle A. Wishard, Jennifer L.L. Pouplin, and Jacob R. Elliott
   !!
   !! Driver program for the Symplectic Massive Body Algorithm
   !!
   !! Adapted from Swifter by David E. Kaufmanna swiftert_symba.f90
   !! Adapted from Hal Levison and Martin Duncan's Swift program swift_symba5.f
   !! Reference: Duncan, M. J., Levison, H. F. & Lee, M. H. 1998. Astron. J., 116, 2067.
   use swiftest

   !> The following are temporary until the conversion to the new module structure is complete
   use module_interfaces
   implicit none

   ! Arguments
   type(user_input_parameters)  :: param    ! derived type containing user-defined parameters
   integer(I4B)      :: nplmax         ! maximum number of planets
   integer(I4B)      :: ntpmax         ! maximum number of test particles
   integer(I4B)      :: istep_out      ! time steps between binary outputs
   integer(I4B)      :: istep_dump     ! time steps between dumps
   real(DP)          :: t0             ! integration start time
   real(DP)          :: tstop          ! integration stop time
   real(DP)          :: dt             ! time step
   real(DP)          :: j2rp2          ! j2*r^2 term for central body
   real(DP)          :: j4rp4          ! j4*r^4 term for central body
   real(DP)          :: rmin           ! minimum heliocentric radius for test particle
   real(DP)          :: rmax           ! maximum heliocentric radius for test particle
   real(DP)          :: rmaxu          ! maximum unbound heliocentric radius for test particle
   real(DP)          :: qmin           ! minimum pericenter distance for test particle
   real(DP)          :: qmin_alo       ! minimum semimajor axis for qmin
   real(DP)          :: qmin_ahi       ! maximum semimajor axis for qmin
   character(strmax) :: qmin_coord     ! coordinate frame to use for qmin
   character(strmax) :: encounter_file ! name of output file for encounters
   character(strmax) :: inplfile       ! name of input file for planets
   character(strmax) :: intpfile       ! name of input file for test particles
   character(strmax) :: in_type        ! format of input data files
   character(strmax) :: outfile        ! name of output binary file
   character(strmax) :: out_type       ! binary format of output file
   character(strmax) :: out_form       ! data to write to output file
   character(strmax) :: out_stat       ! open status for output binary file

   ! Internals
   logical                       :: lfrag_add
   integer(I4B)                  :: npl, ntp, ntp0, nsppl, nsptp, iout, idump, iloop
   integer(I4B)                  :: nplplenc, npltpenc, nmergeadd, nmergesub, fragmax
   real(DP)                      :: t, tfrac, tbase, mtiny, ke, pe, te, eoffset
   real(DP), dimension(ndim)     :: htot
   character(strmax)             :: inparfile
   type(symba_pl)                :: symba_plA
   type(symba_tp)                :: symba_tpA
   type(swiftest_tp)             :: discard_tpA
   type(swiftest_pl)             :: discard_plA
   type(symba_plplenc)           :: plplenc_list
   type(symba_pltpenc)           :: pltpenc_list
   type(symba_merger)            :: mergeadd_list, mergesub_list
   integer(I4B), parameter       :: egyiu = 72
   real(DP)                      :: start, finish

   ! Executable code
   call cpu_time(start)
   call util_version
   nthreads = 1                        
   write(*, 100, advance = "no") "Enter name of parameter data file: "
   read(*, 100) inparfile
   100 format(a)
   inparfile = trim(adjustl(inparfile))
   ! read in the param.in file and get simulation parameters
   call param%read_from_file(inparfile)
   param%lmtiny = .true. ! Turn this on for SyMBA

   ! temporary until the conversion to the derived type argument list is complete
   nplmax = param%nplmax
   ntpmax = param%ntpmax
   t0 = param%t0
   tstop = param%tstop
   dt = param%dt
   inplfile = param%inplfile
   intpfile = param%intpfile
   in_type = param%in_type
   istep_out = param%istep_out
   outfile = param%outfile
   out_type = param%out_type
   out_form = param%out_form
   out_stat = param%out_stat
   istep_dump = param%istep_dump
   j2rp2 = param%j2rp2
   j4rp4 = param%j4rp4
   rmin = param%rmin
   rmax = param%rmax
   rmaxu = param%rmaxu
   qmin = param%qmin
   qmin_coord = param%qmin_coord
   qmin_alo = param%qmin_alo
   qmin_ahi = param%qmin_ahi
   encounter_file = param%encounter_file
   mtiny = param%mtiny
   !^^^^^^^^^^^^^^^^^^^^^^^^^
   if (.not. param%lrhill_present) then
      write(*, *) "Swiftest error:"
      write(*, *) "   Integrator SyMBA requires massive body Hill sphere radii on input"
      call util_exit(failure)
   end if
   ! read in the total number of bodies from the input files

   call symba_plA%read_from_file(param)
   call symba_tpA%read_from_file(param)

   ! Save central body mass in vector form so that elemental functions can be evaluated with it
   call symba_tpA%set_mu(symba_plA%mass(1))
   call symba_plA%set_mu(symba_plA%mass(1))

   !Temporary until everything get switched over
   npl = symba_plA%nbody
   ntp = symba_tpA%nbody

   ! create arrays of data structures big enough to store the number of bodies we are adding
   call mergeadd_list%alloc(10*npl)!DM: Why 10*npl?
   call mergesub_list%alloc(npl)
   call plplenc_list%alloc(10*npl)!DM: See ^
   call pltpenc_list%alloc(ntp)!DM: See ^

   ! reads in initial conditions of all massive bodies from input file
   ! reorder by mass 
   call symba_reorder_pl(npl, symba_plA)
   call util_valid(npl, ntp, symba_plA%helio%swiftest, symba_tpA%helio%swiftest)
   ntp0 = ntp
   t = t0
   tbase = t0
   iloop = 0
   iout = istep_out
   idump = istep_dump
   nmergeadd = 0
   nmergesub = 0
   nsppl = 0
   nsptp = 0
   eoffset = 0.0_DP
   fragmax = 0 
   if (istep_out > 0) then
      call io_write_frame(t, symba_plA%helio%swiftest, symba_tpA%helio%swiftest, outfile, out_type, out_form, out_stat)
   end if
   if (out_stat == "old") then
      open(unit = egyiu, file = energy_file, form = "formatted", status = "old", action = "write", position = "append")
   else 
      open(unit = egyiu, file = energy_file, form = "formatted", status = "replace", action = "write")
   end if
   300 format(7(1x, e23.16))
   write(*, *) " *************** Main Loop *************** "
   if (param%lenergy) then 
      call symba_energy(npl, symba_plA%helio%swiftest, j2rp2, j4rp4, ke, pe, te, htot)
      write(egyiu,300) t, ke, pe, te, htot
   end if
   call symba_energy(npl, symba_plA%helio%swiftest, j2rp2, j4rp4, ke, pe, te, htot)
   do while ((t < tstop) .and. ((ntp0 == 0) .or. (ntp > 0)))
      call symba_step(t, dt, param,npl,ntp,symba_plA, symba_tpA,nplplenc, npltpenc,&
            plplenc_list, pltpenc_list, nmergeadd, nmergesub, mergeadd_list, mergesub_list, &
            eoffset, fragmax)
      iloop = iloop + 1
      if (iloop == loopmax) then
          tbase = tbase + iloop*dt
          iloop = 0
      end if
      t = tbase + iloop*dt
      ldiscard = .false. 
      ldiscard_tp = .false.
      lfrag_add = .false.
      call symba_discard_merge_pl(t, npl, symba_plA, nplplenc, plplenc_list)                                  ! check this 
      call symba_discard_pl(t, npl, nplmax, nsppl, symba_plA, rmin, rmax, rmaxu, qmin, qmin_coord, qmin_alo,    &    ! check this 
            qmin_ahi, j2rp2, j4rp4, eoffset)
      call symba_discard_tp(t, npl, ntp, nsptp, symba_plA, symba_tpA, dt, rmin, rmax, rmaxu, qmin, qmin_coord, &    ! check this 
            qmin_alo, qmin_ahi, param%lclose, param%lrhill_present)
      if ((ldiscard .eqv. .true.) .or. (ldiscard_tp .eqv. .true.) .or. (lfrag_add .eqv. .true.)) then
         call symba_rearray(npl, ntp, nsppl, nsptp, symba_plA, symba_tpA, nmergeadd, mergeadd_list, discard_plA, &
            discard_tpA,param)
         if ((ldiscard .eqv. .true.) .or. (ldiscard_tp .eqv. .true.)) then
            call io_discard_write_symba(t, mtiny, npl, ntp, nsppl, nsptp, nmergeadd, symba_plA, &
               discard_plA, discard_tpA, mergeadd_list, mergesub_list, discard_file, param%lbig_discard) 
            nmergeadd = 0
            nmergesub = 0
            nsppl = 0
            nsptp = 0
         end if 
         if (param%lenergy) then 
            call symba_energy(npl, symba_plA%helio%swiftest, j2rp2, j4rp4, ke, pe, te, htot)
            write(egyiu,300) t, ke, pe, te, htot
         end if
      end if
      if (istep_out > 0) then
         iout = iout - 1
         if (iout == 0) then
            call io_write_frame(t, symba_plA%helio%swiftest, symba_tpA%helio%swiftest, outfile, out_type, out_form, out_stat)
            iout = istep_out
            if (param%lenergy) then 
               call symba_energy(npl, symba_plA%helio%swiftest, j2rp2, j4rp4, ke, pe, te, htot)
               write(egyiu,300) t, ke, pe, te, htot
            end if
            call symba_energy(npl, symba_plA%helio%swiftest, j2rp2, j4rp4, ke, pe, te, htot)
         end if
      end if
      if (istep_dump > 0) then
         idump = idump - 1
         if (idump == 0) then
            tfrac = (t - t0)/(tstop - t0)
            write(*, 200) t, tfrac, npl, ntp
200         format(" Time = ", es12.5, "; fraction done = ", f5.3, "; number of active pl, tp = ", i5, ", ", i5)
            call param%dump_to_file(t)
            call io_dump_pl(npl, symba_plA%helio%swiftest, param%lclose, param%lrhill_present)
            if (ntp > 0) call io_dump_tp(ntp, symba_tpA%helio%swiftest)
            idump = istep_dump
         end if
      end if
      if (allocated(plplenc_list%status)) then
         plplenc_list%status(:) = 0
         plplenc_list%lvdotr(:) = .false.
         plplenc_list%level(:) = 0
         plplenc_list%index1(:) = 0
         plplenc_list%index2(:) = 0
         plplenc_list%enc_child(:) = 0 
         plplenc_list%enc_parent(:) = 0
      end if

      if (allocated(pltpenc_list%status)) then
         pltpenc_list%status(:) = 0
         pltpenc_list%lvdotr(:) = .false.
         pltpenc_list%level(:) = 0
         pltpenc_list%indexpl(:) = 0
         pltpenc_list%indextp(:) = 0
      end if

      if (allocated(mergeadd_list%name)) then
         mergeadd_list%name(:) = 0
         mergeadd_list%index_ps(:) = 0
         mergeadd_list%status(:) = 0
         mergeadd_list%ncomp(:) = 0
         mergeadd_list%xh(:,:) = 0
         mergeadd_list%vh(:,:) = 0
         mergeadd_list%mass(:) = 0
         mergeadd_list%radius(:) = 0
      end if

      if (allocated(mergesub_list%name)) then
         mergesub_list%name(:) = 0
         mergesub_list%index_ps(:) = 0
         mergesub_list%status(:) = 0
         mergesub_list%ncomp(:) = 0
         mergesub_list%xh(:,:) = 0
         mergesub_list%vh(:,:) = 0
         mergesub_list%mass(:) = 0
         mergesub_list%radius(:) = 0
      end if

      if (allocated(discard_plA%name)) call swiftest_pl_deallocate(discard_plA)
      if (allocated(discard_tpA%name)) call swiftest_tp_deallocate(discard_tpA)

   end do
   call param%dump_to_file(t)
   call io_dump_pl(npl, symba_plA%helio%swiftest, param%lclose, param%lrhill_present)
   if (ntp > 0) call io_dump_tp(ntp, symba_tpA%helio%swiftest)
   if (param%lenergy) then 
      call symba_energy(npl, symba_plA%helio%swiftest, j2rp2, j4rp4, ke, pe, te, htot)
      write(egyiu,300) t, ke, pe, te, htot
      close(egyiu)
   end if

   call symba_pl_deallocate(symba_plA)
   call symba_merger_deallocate(mergeadd_list)
   call symba_merger_deallocate(mergesub_list)
   call symba_plplenc_deallocate(plplenc_list)
   if (ntp > 0) then
      call symba_tp_deallocate(symba_tpA)
      call symba_pltpenc_deallocate(pltpenc_list)
   end if
   call cpu_time(finish)
   write(*,*) 'Time: ', finish - start
   call util_exit(SUCCESS)

   stop

end program swiftest_symba
