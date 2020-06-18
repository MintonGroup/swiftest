program swiftest_symba
   !! author: The Purdue Swiftest Team - David A. Minton, Carlisle A. Wishard, Jennifer L.L. Pouplin, and Jacob R. Elliott
   !!
   !! Driver program for the Symplectic Massive Body Algorithm
   !!
   !! Adapted from Swifter by David E. Kaufmanna swiftert_symba.f90
   !! Adapted from Hal Levison and Martin Duncan's Swift program swift_symba5.f
   !! Reference: Duncan, M. J., Levison, H. F. & Lee, M. H. 1998. Astron. J., 116, 2067.
   use swiftest

   implicit none

   ! Arguments
   type(swiftest_configuration)  :: config ! derived type containing user-defined parameters

   ! Internals
   logical                   :: lfrag_add
   integer(I4B)              :: iout, idump, iloop
   real(DP)                  :: t, tfrac, tbase, mtiny, ke, pe, te, eoffset
   real(DP), dimension(NDIM) :: htot
   character(STRMAX)         :: inparfile
   type(symba_pl)            :: symba_plA
   type(symba_tp)            :: symba_tpA
   type(symba_tp)            :: discard_tpA
   type(symba_pl)            :: discard_plA
   type(symba_plplenc)       :: plplenc_list
   type(symba_pltpenc)       :: pltpenc_list
   type(symba_merger)        :: mergeadd_list, mergesub_list
   integer(I4B), parameter   :: egyiu = 72
   real(DP)                  :: start, finish

   ! Executable code
   call cpu_time(start)
   call util_version
   nthreads = 1                        
   write(*, 100, advance = "no") "Enter name of parameter data file: "
   read(*, 100) inparfile
   100 format(a)
   inparfile = trim(adjustl(inparfile))
   ! read in the param.in file and get simulation parameters
   call config%read_from_file(inparfile)
   config%lmtiny = .true. ! Turn this on for SyMBA

   !^^^^^^^^^^^^^^^^^^^^^^^^^
   if (.not. config%lrhill_present) then
      write(*, *) "Swiftest error:"
      write(*, *) "   Integrator SyMBA requires massive body Hill sphere radii on input"
      call util_exit(failure)
   end if
   ! read in the total number of bodies from the input files

   call symba_plA%read_from_file(config)
   call symba_tpA%read_from_file(config)

   ! Save central body mass in vector form so that elemental functions can be evaluated with it
   call symba_tpA%set_vec(symba_plA%mass(1),config%dt)
   call symba_plA%set_vec(symba_plA%mass(1),config%dt)

   ! Save system mass to both objects
   call symba_plA%set_msys(symba_plA)
   call symba_tpA%set_msys(symba_plA)

   ! create arrays of data structures big enough to store the number of bodies we are adding
   call mergeadd_list%alloc(10*npl)!DM: Why 10*npl?
   call mergesub_list%alloc(npl)
   call plplenc_list%alloc(10*npl)!DM: See ^
   call pltpenc_list%alloc(ntp)!DM: See ^

   ! reads in initial conditions of all massive bodies from input file
   ! reorder by mass 
   call symba_reorder_pl(npl, symba_plA)
   call util_valid(npl, ntp, symba_plA, symba_tpA)
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
      call io_write_frame(t, symba_plA, symba_tpA, outfile, out_type, out_form, out_stat)
   end if
   if (out_stat == "old") then
      open(unit = egyiu, file = energy_file, form = "formatted", status = "old", action = "write", position = "append")
   else 
      open(unit = egyiu, file = energy_file, form = "formatted", status = "replace", action = "write")
   end if
   300 format(7(1x, e23.16))
   write(*, *) " *************** Main Loop *************** "
   if (config%lenergy) then 
      call symba_energy(npl, symba_plA, config%j2rp2, config%j4rp4, ke, pe, te, htot)
      write(egyiu,300) t, ke, pe, te, htot
   end if
   call symba_energy(npl, symba_plA, config%j2rp2, config%j4rp4, ke, pe, te, htot)
   do while ((t < config%tstop) .and. ((ntp0 == 0) .or. (ntp > 0)))
      call symba_step(t, dt, symba_plA, symba_tpA, plplenc_list, pltpenc_list, mergeadd_list, mergesub_list, eoffset, config)
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
      call symba_discard_pl(t, npl, config%nplmax, nsppl, symba_plA, rmin, rmax, config%rmaxu, qmin, qmin_coord, qmin_alo,    &    ! check this 
            qmin_ahi, config%j2rp2, config%j4rp4, eoffset)
      call symba_discard_tp(t, npl, ntp, nsptp, symba_plA, symba_tpA, dt, rmin, rmax, config%rmaxu, qmin, qmin_coord, &    ! check this 
            qmin_alo, qmin_ahi, config%lclose, config%lrhill_present)
      if ((ldiscard .eqv. .true.) .or. (ldiscard_tp .eqv. .true.) .or. (lfrag_add .eqv. .true.)) then
         call symba_rearray(npl, ntp, nsppl, nsptp, symba_plA, symba_tpA, nmergeadd, mergeadd_list, discard_plA, &
            discard_tpA,param)
         if ((lspill_list .eqv. .true.) .or. (ldiscard_tp .eqv. .true.)) then
            call io_discard_write_symba(t, mtiny, npl, ntp, nsppl, nsptp, nmergeadd, symba_plA, &
               discard_plA, discard_tpA, mergeadd_list, mergesub_list, discard_file, config%lbig_discard) 
            nmergeadd = 0
            nmergesub = 0
            nsppl = 0
            nsptp = 0
         end if 
         if (config%lenergy) then 
            call symba_energy(npl, symba_plA, config%j2rp2, config%j4rp4, ke, pe, te, htot)
            write(egyiu,300) t, ke, pe, te, htot
         end if
      end if
      if (istep_out > 0) then
         iout = iout - 1
         if (iout == 0) then
            call io_write_frame(t, symba_plA, symba_tpA, outfile, out_type, out_form, out_stat)
            iout = istep_out
            if (config%lenergy) then 
               call symba_energy(npl, symba_plA, config%j2rp2, config%j4rp4, ke, pe, te, htot)
               write(egyiu,300) t, ke, pe, te, htot
            end if
            call symba_energy(npl, symba_plA, config%j2rp2, config%j4rp4, ke, pe, te, htot)
         end if
      end if
      if (istep_dump > 0) then
         idump = idump - 1
         if (idump == 0) then
            tfrac = (t - t0)/(config%tstop - t0)
            write(*, 200) t, tfrac, npl, ntp
200         format(" Time = ", es12.5, "; fraction done = ", f5.3, "; number of active pl, tp = ", i5, ", ", i5)
            call config%dump_to_file(t)
            call io_dump_pl(npl, symba_plA, config%lclose, config%lrhill_present)
            if (ntp > 0) call io_dump_tp(ntp, symba_tpA)
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

      if (allocated(discard_plA%name)) call swiftest_deallocate_pl(discard_plA)
      if (allocated(discard_tpA%name)) call swiftest_deallocate_tp(discard_tpA)

   end do
   call config%dump_to_file(t)
   call io_dump_pl(npl, symba_plA, config%lclose, config%lrhill_present)
   if (ntp > 0) call io_dump_tp(ntp, symba_tpA)
   if (config%lenergy) then 
      call symba_energy(npl, symba_plA, config%j2rp2, config%j4rp4, ke, pe, te, htot)
      write(egyiu,300) t, ke, pe, te, htot
      close(egyiu)
   end if

   call symba_deallocate_pl(symba_plA)
   call symba_deallocate_merger(mergeadd_list)
   call symba_deallocate_merger(mergesub_list)
   call symba_deallocate_plplenc(plplenc_list)
   if (ntp > 0) then
      call symba_deallocate_tp(symba_tpA)
      call symba_deallocate_pltpenc(pltpenc_list)
   end if
   call cpu_time(finish)
   write(*,*) 'Time: ', finish - start
   call util_exit(SUCCESS)

   stop

end program swiftest_symba
