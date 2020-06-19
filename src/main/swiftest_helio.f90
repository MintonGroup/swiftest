program swifter_helio
   !! author: The Purdue Swiftest Team - David A. Minton, Carlisle A. Wishard, Jennifer L.L. Pouplin, and Jacob R. Elliott
   !!
   !! Driver program for Democratic Heliocentric Method
   !!
   !! Adapted from Swifter by David E. Kaufmanna swiftert_helio.f90
   !! Adapted from Hal Levison and Martin Duncan's Swift program swift_helio.f
   !! Reference: Duncan, M. J., Levison, H. F. & Lee, M. H. 1998. Astron. J., 116, 2067.
   use swiftest
   implicit none

   type(swiftest_configuration)  :: config !! Object containing user-defined configuration parameters
   integer(I4B)                  :: narg   !! Number of command line arguments passed
   integer(I4B)                  :: ierr   !! I/O error code 
   logical                       :: lfirst !! Flag indicating that this is the first time through the main loop
   integer(I4B)                  :: iout   
   integer(I4B)                  :: idump
   integer(I4B)                  :: iloop  
   real(DP)                      :: t
   real(DP)                      :: tfrac
   real(DP)                      :: tbase
   character(len=:),allocatable  :: config_file_name
   type(helio_pl)                :: helio_pla !! 
   type(helio_tp)                :: helio_tpa
   real(DP)                      :: start_cpu_time, finish_cpu_time


   call cpu_time(start_cpu_time)
   call util_version ! Splash screen
   config_file_name = io_read_config_file_name()
   call config%read_from_file(config_file_name, integrator = HELIO)


   call set_point(helio_pla)
   if (ntp > 0) then
      allocate(helio_tpa(ntpmax))
      call set_point(helio_tpa)
   end if
   call helio_setup(npl, ntp, helio_pla, helio_tpa, helio_pl1p, helio_tp1p, swifter_pl1p, swifter_tp1p)
   call io_init_pl(inplfile, in_type, lclose, lrhill_present, npl, swifter_pl1p)
   call io_init_tp(intpfile, in_type, ntp, swifter_tp1p)
   call util_valid(npl, ntp, swifter_pl1p, swifter_tp1p)
   lfirst = .true.
   ntp0 = ntp
   t = t0
   tbase = t0
   iloop = 0
   iout = istep_out
   idump = istep_dump
   nsp = 0
   eoffset = 0.0_DP
   nullify(helio_tpd1p)
   if (istep_out > 0) call io_write_frame(t, npl, ntp, swifter_pl1p, swifter_tp1p, outfile, out_type, out_form, out_stat)
   write(*, *) " *************** main loop *************** "
   do while ((t < tstop) .and. ((ntp0 == 0) .or. (ntp > 0)))
      call helio_step(lfirst, lextra_force, t, npl, nplmax, ntp, ntpmax, helio_pl1p, helio_tp1p, j2rp2, j4rp4, dt)
      iloop = iloop + 1
      if (iloop == LOOPMAX) then
         tbase = tbase + iloop*dt
         iloop = 0
      end if
      t = tbase + iloop*dt
      call helio_discard(t, npl, ntp, nsp, helio_pl1p, helio_tp1p, helio_tpd1p, dt, rmin, rmax, rmaxu, qmin, qmin_coord,    &
         qmin_alo, qmin_ahi, lclose, lrhill_present)
      if (nsp > 0) then
         swifter_tp1p => helio_tp1p%swifter
         swifter_tpd1p => helio_tpd1p%swifter
         call io_discard_write(t, npl, nsp, swifter_pl1p, swifter_tpd1p, discard_file, lbig_discard)
         nsp = 0
         nullify(helio_tpd1p)
      end if
      if (istep_out > 0) then
         iout = iout - 1
         if (iout == 0) then
            call io_write_frame(t, npl, ntp, swifter_pl1p, swifter_tp1p, outfile, out_type, out_form, out_stat)
            iout = istep_out
         end if
      end if
      if (istep_dump > 0) then
         idump = idump - 1
         if (idump == 0) then
            tfrac = (t - t0)/(tstop - t0)
            write(*, 200) t, tfrac, npl, ntp
 200          format(" time = ", es12.5, "; fraction done = ", f5.3, "; number of active pl, tp = ", i5, ", ", i5)
            call io_dump_param(nplmax, ntpmax, ntp, t, tstop, dt, in_type, istep_out, outfile, out_type, out_form,      &
               istep_dump, j2rp2, j4rp4, lclose, rmin, rmax, rmaxu, qmin, qmin_coord, qmin_alo, qmin_ahi,         &
               encounter_file, lextra_force, lbig_discard, lrhill_present)
            call io_dump_pl(npl, swifter_pl1p, lclose, lrhill_present)
            if (ntp > 0) call io_dump_tp(ntp, swifter_tp1p)
            idump = istep_dump
         end if
      end if
   end do
   call io_dump_param(nplmax, ntpmax, ntp, t, tstop, dt, in_type, istep_out, outfile, out_type, out_form, istep_dump, j2rp2,    &
      j4rp4, lclose, rmin, rmax, rmaxu, qmin, qmin_coord, qmin_alo, qmin_ahi, encounter_file, lextra_force, lbig_discard,   &
      lrhill_present)
   call io_dump_pl(npl, swifter_pl1p, lclose, lrhill_present)
   if (ntp > 0) call io_dump_tp(ntp, swifter_tp1p)
   if (allocated(helio_pla)) deallocate(helio_pla)
   if (allocated(helio_tpa)) deallocate(helio_tpa)


   call cpu_time(finish_cpu_time)
   write(*,*) 'Time: ', finish_cpu_time - start_cpu_time
   call util_exit(SUCCESS)

   stop

end program swifter_helio
