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

   type(swiftest_configuration)  :: config !! derived type containing user-defined parameters
   logical                   :: lfrag_add
   integer(I4B)              :: iout, idump, iloop
   real(DP)                  :: t, tfrac, tbase
   real(DP), dimension(NDIM) :: htot
   character(STRMAX)         :: config_file_name
   type(symba_pl)            :: symba_plA
   type(symba_tp)            :: symba_tpA
   type(symba_tp)            :: discard_tpA
   type(symba_pl)            :: discard_plA
   type(symba_plplenc)       :: plplenc_list
   type(symba_pltpenc)       :: pltpenc_list
   type(symba_merger)        :: mergeadd_list, mergesub_list
   integer(I4B), parameter   :: egyiu = 72
   real(DP)                  :: start_cpu_time, finish_cpu_time

   call cpu_time(start_cpu_time)
   call util_version ! Splash screen
   config_file_name = io_read_config_file_name()
   call config%read_from_file(config_file_name, integrator = SYMBA)

   call symba_set_initial_conditions(symba_plA, symba_tpA, config)

   t = t0
   tbase = t0
   iloop = 0
   iout = config%istep_out
   idump = config%istep_dump
   eoffset = 0.0_DP
   if (istep_out > 0) then
      call io_write_frame(t, symba_plA, symba_tpA, config)
   end if
   ! TODO: Replace with subroutine call
   !if (out_stat == "old") then
   !   open(unit = egyiu, file = ENERGY_FILE, form = "formatted", status = "old", action = "write", position = "append")
   !else 
   !   open(unit = egyiu, file = ENERGY_FILE, form = "formatted", status = "replace", action = "write")
   !end if
   !300 format(7(1x, e23.16))


   write(*, *) " *************** Main Loop *************** "
   call symba_plA%calc_conserved()
   !if (config%lenergy) write(egyiu,300) t, ke, pe, te, htot

   do while (t < config%tstop)
      call symba_step(t, dt, symba_plA, symba_tpA, plplenc_list, pltpenc_list, mergeadd_list, mergesub_list, config)
      iloop = iloop + 1
      if (iloop == LOOPMAX) then
          tbase = tbase + iloop*dt
          iloop = 0
      end if
      t = tbase + iloop * dt
      call symba_discard_merge_pl(t, symba_plA, plplenc_list)                                  
      call symba_plA%discard(t, dt, eoffset, config)
      call symba_tpA%discard(t, dt, symba_plA, config) 
      if ((symba_plA%ldiscard) .or. (symba_tpA%ldiscard) .or. (lfrag_add .eqv. .true.)) then
         call symba_rearray(npl, ntp, nsppl, nsptp, symba_plA, symba_tpA, nmergeadd, mergeadd_list, discard_plA, &
            discard_tpA,param)
         !if ((lspill_list) .or. (ldiscard_tpt
            !call io_discard_write_symba(t, mtiny, npl, ntp, nsppl, nsptp, nmergeadd, symba_plA, &
            !   discard_plA, discard_tpA, mergeadd_list, mergesub_list, discard_file, config%lbig_discard) 
            !nmergeadd = 0
            !nmergesub = 0
            !nsppl = 0
            !nsptp = 0
         end if 
         !if (config%lenergy) then 
         !   call symba_plA%calc_conserved()
         !   write(egyiu,300) t, ke, pe, te, htot
         !end if
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

      if (allocated(discard_plA%name)) call nbody_deallocate_pl(discard_plA)
      if (allocated(discard_tpA%name)) call nbody_deallocate_tp(discard_tpA)

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
   call cpu_time(finish_cpu_time)
   write(*,*) 'Time: ', finish_cpu_time - start_cpu_time
   call util_exit(SUCCESS)

   stop

end program swiftest_symba
