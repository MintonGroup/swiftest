program swifter_helio
   !! author: The Purdue Swiftest Team - David A. Minton, Carlisle A. Wishard, Jennifer L.L. Pouplin, and Jacob R. Elliott
   !!
   !! Driver program for Democratic Heliocentric Method
   !!
   !! Adapted from Swifter by David E. Kaufmanna swiftert_helio.f90
   !! Adapted from Hal Levison and Martin Duncan's Swift program swift_helio.f
   !! Reference: Duncan, M. J., Levison, H. F. & Lee, M. H. 1998. Astron. J., 116, 2067.

! modules
   use swiftest
   implicit none

! arguments
   logical     :: lclose       ! check for planet-test particle encounters
   logical     :: lextra_force   ! use user-supplied force routines
   logical     :: lbig_discard   ! dump planet data with discards
   logical     :: lrhill_present ! hill's sphere radius present
   integer(I4B)    :: nplmax       ! maximum number of planets
   integer(I4B)    :: ntpmax       ! maximum number of test particles
   integer(I4B)    :: istep_out    ! time steps between binary outputs
   integer(I4B)    :: istep_dump   ! time steps between dumps
   real(DP)      :: t0         ! integration start time
   real(DP)      :: tstop      ! integration stop time
   real(DP)      :: dt         ! time step
   real(DP)      :: j2rp2      ! j2*r^2 term for central body
   real(DP)      :: j4rp4      ! j4*r^4 term for central body
   real(DP)      :: rmin       ! minimum heliocentric radius for test particle
   real(DP)      :: rmax       ! maximum heliocentric radius for test particle
   real(DP)      :: rmaxu      ! maximum unbound heliocentric radius for test particle
   real(DP)      :: qmin       ! minimum pericenter distance for test particle
   real(DP)      :: qmin_alo     ! minimum semimajor axis for qmin
   real(DP)      :: qmin_ahi     ! maximum semimajor axis for qmin
   character(STRMAX) :: qmin_coord   ! coordinate frame to use for qmin
   character(STRMAX) :: encounter_file ! name of output file for encounters
   character(STRMAX) :: inplfile     ! name of input file for planets
   character(STRMAX) :: intpfile     ! name of input file for test particles
   character(STRMAX) :: in_type      ! format of input data files
   character(STRMAX) :: outfile      ! name of output binary file
   character(STRMAX) :: out_type     ! binary format of output file
   character(STRMAX) :: out_form     ! data to write to output file
   character(STRMAX) :: out_stat     ! open status for output binary file

! internals
   logical                         :: lfirst
   integer(I4B)                        :: npl, ntp, ntp0, nsp, iout, idump, iloop
   real(DP)                          :: t, tfrac, tbase, eoffset
   character(STRMAX)                     :: inparfile
   type(swifter_pl), pointer               :: swifter_pl1p
   type(swifter_tp), pointer               :: swifter_tp1p, swifter_tpd1p
   type(helio_pl), dimension(:), allocatable, target :: helio_pla
   type(helio_tp), dimension(:), allocatable, target :: helio_tpa
   type(helio_pl), pointer                 :: helio_pl1p
   type(helio_tp), pointer                 :: helio_tp1p, helio_tpd1p

! executable code
   call util_version
   write(*, 100, advance = "no") "enter name of parameter data file: "
   read(*, 100) inparfile
 100 format(a)
   inparfile = trim(adjustl(inparfile))
   call io_init_param(inparfile, nplmax, ntpmax, t0, tstop, dt, inplfile, intpfile, in_type, istep_out, outfile, out_type,    &
      out_form, out_stat, istep_dump, j2rp2, j4rp4, lclose, rmin, rmax, rmaxu, qmin, qmin_coord, qmin_alo, qmin_ahi,      &
      encounter_file, lextra_force, lbig_discard, lrhill_present)
   call io_getn(inplfile, intpfile, in_type, npl, nplmax, ntp, ntpmax)
   allocate(helio_pla(nplmax))
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
   call util_exit(success)

   stop

end program swifter_helio
