module module_interfaces

     implicit none

     interface
          subroutine coord_h2b(npl, swiftest_pla, msys)
               use swiftest_globals
               use swiftest_data_structures
               implicit none
               integer(I4B), intent(in)         :: npl
               real(DP), intent(out)            :: msys
               type(swiftest_pl), intent(inout) :: swiftest_pla
          end subroutine coord_h2b
     end interface

     interface
          subroutine coord_h2b_tp(ntp, swiftest_tpa, swiftest_pla)
               use swiftest_globals
               use swiftest_data_structures
               implicit none
               integer(I4B), intent(in)         :: ntp
               type(swiftest_tp), intent(inout) :: swiftest_tpa
               type(swiftest_pl), intent(inout) :: swiftest_pla
          end subroutine coord_h2b_tp
     end interface

     interface
          subroutine coord_vb2vh(npl, swiftest_pla)
               use swiftest_globals
               use swiftest_data_structures
               implicit none
               integer(I4B), intent(in)         :: npl
               type(swiftest_pl), intent(inout) :: swiftest_pla
          end subroutine coord_vb2vh
     end interface

     interface
          subroutine coord_vb2vh_tp(ntp, swiftest_tpa, vs)
               use swiftest_globals
               use swiftest_data_structures
               implicit none
               integer(I4B), intent(in)              :: ntp
               real(DP), dimension(ndim), intent(in) :: vs
               type(swiftest_tp), intent(inout)      :: swiftest_tpa
          end subroutine coord_vb2vh_tp
     end interface

     interface
          subroutine coord_vh2vb(npl, swiftest_pla, msys)
               use swiftest_globals
               use swiftest_data_structures
               implicit none
               integer(I4B), intent(in)        :: npl
               real(DP), intent(out)           :: msys
               type(swiftest_pl),intent(inout) :: swiftest_pla
          end subroutine coord_vh2vb
     end interface

     interface
          subroutine coord_vh2vb_tp(ntp, swiftest_tpa, vs)
               use swiftest_globals
               use swiftest_data_structures
               implicit none
               integer(I4B), intent(in)              :: ntp
               real(DP), dimension(ndim), intent(in) :: vs
               type(swiftest_tp), intent(inout)      :: swiftest_tpa
          end subroutine coord_vh2vb_tp
     end interface

     interface
          subroutine discard(t, dt, npl, ntp, swiftest_pla, swiftest_tpa, rmin, rmax, rmaxu, qmin,  &
               qmin_alo, qmin_ahi, qmin_coord, lclose, lrhill_present)
               use swiftest_globals
               use swiftest_data_structures
               implicit none
               logical(lgt), intent(in)  :: lclose, lrhill_present
               integer(I4B), intent(in)  :: npl, ntp
               real(DP), intent(in)      :: t, dt, rmin, rmax, rmaxu, qmin, qmin_alo, qmin_ahi
               character(*), intent(in)  :: qmin_coord
               type(swiftest_pl), intent(inout) :: swiftest_pla
               type(swiftest_tp), intent(inout) :: swiftest_tpa
          end subroutine discard
     end interface

     interface
          subroutine discard_peri(t, npl, ntp, swiftest_pla, swiftest_tpa, msys, qmin, qmin_alo, & 
               qmin_ahi, qmin_coord, lrhill_present)
               use swiftest_globals
               use swiftest_data_structures
               implicit none
               logical(lgt), intent(in)  :: lrhill_present
               integer(I4B), intent(in)  :: npl, ntp
               real(DP), intent(in)      :: t, msys, qmin, qmin_alo, qmin_ahi
               character(*), intent(in)  :: qmin_coord
               type(swiftest_pl), intent(inout) :: swiftest_pla
               type(swiftest_tp), intent(inout) :: swiftest_tpa
          end subroutine discard_peri
     end interface

     interface
          subroutine discard_pl_close(dx, dv, dt, r2crit, iflag, r2min)
               use swiftest_globals
               implicit none
               integer(I4B), intent(out)             :: iflag
               real(DP), intent(in)                  :: dt, r2crit
               real(DP), dimension(ndim), intent(in) :: dx, dv
               real(DP), intent(out)                 :: r2min
          end subroutine discard_pl_close
     end interface

     interface
          subroutine discard_pl(t, dt, npl, ntp, swiftest_pla, swiftest_tpa)
               use swiftest_globals
               use swiftest_data_structures
               implicit none
               integer(I4B), intent(in)  :: npl, ntp
               real(DP), intent(in)      :: t, dt
               type(swiftest_pl), intent(inout) :: swiftest_pla
               type(swiftest_tp), intent(inout) :: swiftest_tpa
          end subroutine discard_pl
     end interface

     interface
          subroutine discard_sun(t, ntp, msys, swifter_tpa, rmin, rmax, rmaxu)
               use swiftest_globals
               use swiftest_data_structures
               implicit none
               integer(I4B), intent(in)  :: ntp
               real(DP), intent(in)      :: t, msys, rmin, rmax, rmaxu
               type(swiftest_tp), intent(inout) :: swifter_tpa
          end subroutine discard_sun
     end interface

     interface
          subroutine drift_dan(mu, x0, v0, dt0, iflag)
               use swiftest_globals
               implicit none
               integer(I4B), intent(out)                :: iflag
               real(DP), intent(in)                     :: mu, dt0
               real(DP), dimension(ndim), intent(inout) :: x0, v0
          end subroutine drift_dan
     end interface

     interface
          subroutine drift_kepmd(dm, es, ec, x, s, c)
               use swiftest_globals
               implicit none
               real(DP), intent(in)  :: dm, es, ec
               real(DP), intent(out) :: x, s, c
          end subroutine drift_kepmd
     end interface

     interface
          subroutine drift_kepu(dt,r0,mu,alpha,u,fp,c1,c2,c3,iflag)
               use swiftest_globals
               implicit none
               integer(I4B), intent(out) :: iflag
               real(DP), intent(in)      :: dt, r0, mu, alpha, u
               real(DP), intent(out)     :: fp, c1, c2, c3
          end subroutine drift_kepu
     end interface

     interface
          subroutine drift_kepu_fchk(dt, r0, mu, alpha, u, s, f)
               use swiftest_globals
               implicit none
               real(DP), intent(in)  :: dt, r0, mu, alpha, u, s
               real(DP), intent(out) :: f
          end subroutine drift_kepu_fchk
     end interface

     interface
          subroutine drift_kepu_guess(dt, r0, mu, alpha, u, s)
               use swiftest_globals
               implicit none
               real(DP), intent(in)  :: dt, r0, mu, alpha, u
               real(DP), intent(out) :: s
          end subroutine drift_kepu_guess
     end interface

     interface
          subroutine drift_kepu_lag(s, dt, r0, mu, alpha, u, fp, c1, c2, c3, iflag)
               use swiftest_globals
               implicit none
               integer(I4B), intent(out) :: iflag
               real(DP), intent(in)      :: dt, r0, mu, alpha, u
               real(DP), intent(inout)   :: s
               real(DP), intent(out)     :: fp, c1, c2, c3
          end subroutine drift_kepu_lag
     end interface

     interface
          subroutine drift_kepu_new(s, dt, r0, mu, alpha, u, fp, c1, c2, c3, iflag)
               use swiftest_globals
               implicit none
               integer(I4B), intent(out) :: iflag
               real(DP), intent(in)      :: dt, r0, mu, alpha, u
               real(DP), intent(inout)   :: s
               real(DP), intent(out)     :: fp, c1, c2, c3
          end subroutine drift_kepu_new
     end interface

     interface
          subroutine drift_kepu_p3solve(dt, r0, mu, alpha, u, s, iflag)
               use swiftest_globals
               implicit none
               integer(I4B), intent(out) :: iflag
               real(DP), intent(in)      :: dt, r0, mu, alpha, u
               real(DP), intent(out)     :: s
          end subroutine drift_kepu_p3solve
     end interface

     interface
          subroutine drift_kepu_stumpff(x, c0, c1, c2, c3)
               use swiftest_globals
               implicit none
               real(DP), intent(inout) :: x
               real(DP), intent(out)   :: c0, c1, c2, c3
          end subroutine drift_kepu_stumpff
     end interface

     interface
          subroutine drift_one(mu, x, v, dt, iflag)
               use swiftest_globals
               implicit none
               integer(I4B), intent(out)                :: iflag
               real(DP), intent(in)                     :: mu, dt
               real(DP), dimension(ndim), intent(inout) :: x, v
          end subroutine drift_one
     end interface



     interface
          subroutine io_discard_write_symba(t, mtiny, npl, ntp, nsppl, nsptp, nmergeadd, nmergesub, symba_pla, & 
               discard_pla, discard_tpa, mergeadd_list, mergesub_list, fname, lbig_discard)
               use swiftest_globals
               use swiftest_data_structures
               use symba
               implicit none
               logical(lgt), intent(in)                       :: lbig_discard
               integer(I4B), intent(in)                       :: npl, ntp, nsppl, nsptp, nmergeadd, nmergesub
               real(DP), intent(in)                           :: t, mtiny
               character(*), intent(in)                       :: fname
               type(symba_pl), intent(inout)                  :: symba_pla
               type(swiftest_tp), intent(inout)               :: discard_tpa
               type(swiftest_pl), intent(inout)               :: discard_pla
               type(symba_merger), intent(inout)              :: mergeadd_list, mergesub_list
          end subroutine io_discard_write_symba
     end interface

     interface
          subroutine io_dump_param(nplmax, ntpmax, ntp, t, tstop, dt, in_type, istep_out, outfile, out_type, out_form,            &
               istep_dump, j2rp2, j4rp4, lclose, rmin, rmax, rmaxu, qmin, qmin_coord, qmin_alo, qmin_ahi, encounter_file,         &
               lextra_force, lbig_discard, lrhill_present, mtiny, lpython)
               use swiftest_globals
               implicit none
               logical(lgt), intent(in) :: lclose, lextra_force, lbig_discard, lrhill_present, lpython
               integer(I4B), intent(in) :: nplmax, ntpmax, ntp, istep_out, istep_dump
               real(DP), intent(in)     :: t, tstop, dt, j2rp2, j4rp4, rmin, rmax, rmaxu, qmin, qmin_alo, qmin_ahi, mtiny
               character(*), intent(in) :: qmin_coord, encounter_file, in_type, outfile, out_type, out_form
          end subroutine io_dump_param
     end interface

     interface
          subroutine io_dump_pl(npl, swiftest_pla, lclose, lrhill_present)
               use swiftest_globals
               use swiftest_data_structures
               implicit none
               logical(lgt), intent(in)        :: lclose, lrhill_present
               integer(I4B), intent(in)        :: npl
               type(swiftest_pl), intent(inout):: swiftest_pla
          end subroutine io_dump_pl
     end interface

     interface
          subroutine io_dump_tp(ntp, swiftest_tpa)
               use swiftest_globals
               use swiftest_data_structures
               implicit none
               integer(I4B), intent(in)        :: ntp
               type(swiftest_tp), intent(inout):: swiftest_tpa
          end subroutine io_dump_tp
     end interface

     interface
          subroutine io_getn(inplfile, intpfile, in_type, npl, nplmax, ntp, ntpmax)
               use swiftest_globals
               implicit none
               integer(I4B), intent(inout) :: nplmax, ntpmax
               integer(I4B), intent(out)   :: npl, ntp
               character(*), intent(in)    :: inplfile, intpfile, in_type
          end subroutine io_getn
     end interface

     interface
          subroutine io_get_token(buffer, ilength, ifirst, ilast, ierr)
               use swiftest_globals
               implicit none
               integer(I4B), intent(in)    :: ilength
               integer(I4B), intent(inout) :: ifirst
               integer(I4B), intent(out)   :: ilast, ierr
               character(*), intent(in)    :: buffer
          end subroutine io_get_token
     end interface

     interface
          subroutine io_init_param(inparfile, nplmax, ntpmax, t0, tstop, dt, inplfile, intpfile, in_type, istep_out, outfile,     &
               out_type, out_form, out_stat, istep_dump, j2rp2, j4rp4, lclose, rmin, rmax, rmaxu, qmin, qmin_coord, qmin_alo,     &
               qmin_ahi, encounter_file, lextra_force, lbig_discard, lrhill_present, mtiny, lpython, lenergy)
               use swiftest_globals
               implicit none
               logical(lgt), intent(out) :: lclose, lextra_force, lbig_discard, lrhill_present, lpython, lenergy
               integer(I4B), intent(out) :: nplmax, ntpmax, istep_out, istep_dump
               real(DP), intent(out)     :: t0, tstop, dt, j2rp2, j4rp4, rmin, rmax, rmaxu, qmin, qmin_alo, qmin_ahi
               real(DP), intent(out), optional :: mtiny
               character(*), intent(in)  :: inparfile
               character(*), intent(out) :: qmin_coord, encounter_file, inplfile, intpfile, in_type, outfile, out_type, out_form, &
                                            out_stat
          end subroutine io_init_param
     end interface

     interface
          subroutine io_init_pl(inplfile, in_type, lclose, lrhill_present, npl, symba_pla)
               use swiftest_globals
               use swiftest_data_structures
               use symba
               use helio
               implicit none
               logical(lgt), intent(in)         :: lclose, lrhill_present
               integer(I4B), intent(in)         :: npl
               character(*), intent(in)         :: inplfile, in_type
               type(symba_pl), intent(inout)    :: symba_pla
          end subroutine io_init_pl
     end interface

     interface
          subroutine io_init_tp(intpfile, in_type, ntp, symba_tpa)
               use swiftest_globals
               use swiftest_data_structures
               use symba
               use helio
               implicit none
               integer(I4B), intent(in)         :: ntp
               character(*), intent(in)         :: intpfile, in_type
               type(symba_tp), intent(inout)    :: symba_tpa
          end subroutine io_init_tp
     end interface

     interface
          subroutine io_open(iu, fname, fopenstat, fmt, ierr)
               use swiftest_globals
               implicit none
               integer(I4B), intent(in)  :: iu
               integer(I4B), intent(out) :: ierr
               character(*), intent(in)  :: fname, fopenstat, fmt
          end subroutine io_open
     end interface

     interface
          function io_read_encounter(t, name1, name2, mass1, mass2, xh1, xh2, vh1, vh2, encounter_file, out_type)
               use swiftest_globals
               implicit none
               integer(I4B)                           :: io_read_encounter
               integer(I4B), intent(out)              :: name1, name2
               real(DP), intent(out)                  :: t, mass1, mass2
               real(DP), dimension(ndim), intent(out) :: xh1, xh2, vh1, vh2
               character(*), intent(in)               :: encounter_file,out_type
          end function io_read_encounter
     end interface

     interface
          function io_read_hdr(iu, t, npl, ntp, iout_form, out_type)
               use swiftest_globals
               implicit none
               integer(I4B)               :: io_read_hdr
               integer(I4B), intent(in)   :: iu
               integer(I4B), intent(out)  :: npl, ntp, iout_form
               real(DP), intent(out)      :: t
               character(*), intent(in)   :: out_type
          end function io_read_hdr
     end interface

     interface
          function io_read_line(iu, name, d1, d2, d3, d4, d5, d6, out_type, mass, radius)
               use swiftest_globals
               implicit none
               integer(I4B)                    :: io_read_line
               integer(I4B), intent(in)        :: iu
               integer(I4B), intent(out)       :: name
               real(DP), intent(out)           :: d1, d2, d3, d4, d5, d6
               real(DP), optional, intent(out) :: mass, radius
               character(*), intent(in)        :: out_type
          end function io_read_line
     end interface

     interface
          subroutine io_write_encounter(t, name1, name2, mass1, mass2, radius1, radius2, &
               xh1, xh2, vh1, vh2, encounter_file, out_type)
               use swiftest_globals
               implicit none
               integer(I4B), intent(in)              :: name1, name2
               real(DP), intent(in)                  :: t, mass1, mass2, radius1, radius2
               real(DP), dimension(ndim), intent(in) :: xh1, xh2, vh1, vh2
               character(*), intent(in)              :: encounter_file, out_type
          end subroutine io_write_encounter
     end interface

     interface
          subroutine io_write_frame(t, npl, ntp, swiftest_pla, swiftest_tpa, outfile, &
               out_type, out_form, out_stat)
               use swiftest_globals
               use swiftest_data_structures
               implicit none
               integer(I4B), intent(in)  :: npl, ntp
               real(DP), intent(in)      :: t
               character(*), intent(in)  :: outfile, out_type, out_form, out_stat
               type(swiftest_pl), intent(inout) :: swiftest_pla
               type(swiftest_tp), intent(inout) :: swiftest_tpa
          end subroutine io_write_frame
     end interface

     interface
          subroutine io_write_hdr(iu, t, npl, ntp, iout_form, out_type)
               use swiftest_globals
               implicit none
               integer(I4B), intent(in) :: iu, npl, ntp, iout_form
               real(DP), intent(in)     :: t
               character(*), intent(in) :: out_type
          end subroutine io_write_hdr
     end interface

     interface
          subroutine io_write_line(iu, name, d1, d2, d3, d4, d5, d6, out_type, mass, radius)
               use swiftest_globals
               implicit none
               integer(I4B), intent(in)       :: iu, name
               real(DP), intent(in)           :: d1, d2, d3, d4, d5, d6
               real(DP), optional, intent(in) :: mass, radius
               character(*), intent(in)       :: out_type
          end subroutine io_write_line
     end interface

     interface
          subroutine obl_acc(npl, swiftest_pla, j2rp2, j4rp4, xh, irh, aobl)
               use swiftest_globals
               use swiftest_data_structures
               implicit none
               integer(I4B), intent(in)                    :: npl
               real(DP), intent(in)                        :: j2rp2, j4rp4
               real(DP), dimension(npl), intent(in)        :: irh
               real(DP), dimension(ndim, npl), intent(in)  :: xh
               real(DP), dimension(ndim, npl), intent(out) :: aobl
               type(swiftest_pl), intent(inout)            :: swiftest_pla
          end subroutine obl_acc
     end interface

     interface
          subroutine obl_acc_tp(ntp, xht, j2rp2, j4rp4, irht, aoblt, msun)
               use swiftest_globals
               implicit none
               integer(I4B), intent(in)                    :: ntp
               real(DP), intent(in)                        :: j2rp2, j4rp4, msun
               real(DP), dimension(ntp), intent(in)        :: irht
               real(DP), dimension(ndim, ntp), intent(in)  :: xht
               real(DP), dimension(ndim, ntp), intent(out) :: aoblt
          end subroutine obl_acc_tp
     end interface

     interface
          subroutine obl_pot(npl, swiftest_pla, j2rp2, j4rp4, xh, irh, oblpot)
               use swiftest_globals
               use swiftest_data_structures
               implicit none
               integer(I4B), intent(in)                   :: npl
               real(DP), intent(in)                       :: j2rp2, j4rp4
               real(DP), intent(out)                      :: oblpot
               real(DP), dimension(npl), intent(in)       :: irh
               real(DP), dimension(ndim, npl), intent(in) :: xh
               type(swiftest_pl), intent(inout)           :: swiftest_pla
          end subroutine obl_pot
     end interface

     interface
          subroutine orbel_scget(angle, sx, cx)
               use swiftest_globals
               implicit none
               real(DP), intent(in)  :: angle
               real(DP), intent(out) :: sx, cx
          end subroutine orbel_scget
     end interface

     interface
          subroutine orbel_xv2aeq(x, v, mu, a, e, q)
               use swiftest_globals
               implicit none
               real(DP), intent(in)                  :: mu
               real(DP), dimension(ndim), intent(in) :: x, v
               real(DP), intent(out)                 :: a, e, q
          end subroutine orbel_xv2aeq
     end interface

     interface
          subroutine orbel_xv2aqt(x, v, mu, a, q, capm, tperi)
               use swiftest_globals
               implicit none
               real(DP), intent(in)                  :: mu
               real(DP), dimension(ndim), intent(in) :: x, v
               real(DP), intent(out)                 :: a, q, capm, tperi
          end subroutine orbel_xv2aqt
     end interface

     interface
          subroutine orbel_xv2el(x, v, mu, a, e, inc, capom, omega, capm)
               use swiftest_globals
               implicit none
               real(DP), intent(in)                  :: mu
               real(DP), dimension(ndim), intent(in) :: x, v
               real(DP), intent(out)                 :: a, e, inc, capom, omega, capm
          end subroutine orbel_xv2el
     end interface

     interface
          subroutine python_io_write_frame_pl(t, symba_pla, npl, out_stat)
               use swiftest_globals
               use swiftest_data_structures
               use helio
               use symba
               implicit none
               real(DP), intent(in)      :: t
               type(symba_pl),intent(in) :: symba_pla
               integer, intent(in)       :: npl
               character(*), intent(in)  :: out_stat
          end subroutine python_io_write_frame_pl
     end interface

     interface
          subroutine python_io_write_frame_tp(t, symba_tpa, ntp, out_stat)
               use swiftest_globals
               use swiftest_data_structures
               use helio
               use symba
               implicit none
               real(DP), intent(in)      :: t
               type(symba_tp),intent(in) :: symba_tpa
               integer, intent(in)       :: ntp
               character(*), intent(in)  :: out_stat
          end subroutine python_io_write_frame_tp
     end interface



     interface
          subroutine rmvs_chk_ind(xr, vr, dt, r2crit, iflag)
               use swiftest_globals
               implicit none
               real(DP), intent(in)                  :: dt, r2crit
               real(DP), dimension(ndim), intent(in) :: xr, vr
               integer(I4B), intent(out)             :: iflag
          end subroutine rmvs_chk_ind
     end interface

     interface
          subroutine symba_casedisruption (t, dt, index_enc, nmergeadd, nmergesub, mergeadd_list, mergesub_list, eoffset, vbs, & 
          encounter_file, out_type, npl, ntp, symba_pla, symba_tpa, nplplenc, npltpenc, pltpenc_list, plplenc_list, &
          nplmax, ntpmax, fragmax, mres, rres, array_index1_child, array_index2_child, m1, m2, rad1, rad2, x1, x2, v1, v2)
          use swiftest_globals
          use symba
          implicit none
          integer(I4B), intent(in)                         :: index_enc, nplmax, ntpmax
          integer(I4B), intent(inout)                      :: npl, ntp, nmergeadd, nmergesub, nplplenc, npltpenc, fragmax
          real(DP), intent(in)                             :: t, dt
          real(DP), intent(inout)                          :: eoffset, m1, m2, rad1, rad2
          real(DP), dimension(3), intent(inout)            :: mres, rres
          real(DP), dimension(ndim), intent(in)            :: vbs
          real(DP), dimension(ndim), intent(inout)         :: x1, x2, v1, v2
          character(*), intent(in)                         :: encounter_file, out_type
          type(symba_plplenc), intent(inout)               :: plplenc_list
          type(symba_pltpenc), intent(inout)               :: pltpenc_list
          type(symba_merger), intent(inout)                :: mergeadd_list, mergesub_list
          type(symba_pl), intent(inout)                    :: symba_pla
          type(symba_tp), intent(inout)                    :: symba_tpa
          integer(I4B), dimension(npl), intent(inout)      :: array_index1_child, array_index2_child

          end subroutine symba_casedisruption
     end interface

     interface
          subroutine symba_casehitandrun (t, dt, index_enc, nmergeadd, nmergesub, mergeadd_list, mergesub_list, eoffset, vbs, & 
          encounter_file, out_type, npl, ntp, symba_pla, symba_tpa, nplplenc, npltpenc, pltpenc_list, plplenc_list, &
          nplmax, ntpmax, fragmax, mres, rres, array_index1_child, array_index2_child, m1, m2, rad1, rad2, x1, x2, v1, v2)
          use swiftest_globals
          use swiftest_data_structures
          use helio
          use symba
          implicit none
          integer(I4B), intent(in)                         :: index_enc, nplmax, ntpmax
          integer(I4B), intent(inout)                      :: npl, ntp, nmergeadd, nmergesub, nplplenc, npltpenc, fragmax
          real(DP), intent(in)                             :: t, dt
          real(DP), intent(inout)                          :: eoffset, m1, m2, rad1, rad2
          real(DP), dimension(3), intent(inout)            :: mres, rres
          real(DP), dimension(ndim), intent(in)            :: vbs
          real(DP), dimension(ndim), intent(inout)         :: x1, x2, v1, v2
          character(*), intent(in)                         :: encounter_file, out_type
          type(symba_plplenc), intent(inout)               :: plplenc_list
          type(symba_pltpenc), intent(inout)               :: pltpenc_list
          type(symba_merger), intent(inout)                :: mergeadd_list, mergesub_list
          type(symba_pl), intent(inout)                    :: symba_pla
          type(symba_tp), intent(inout)                    :: symba_tpa
          integer(I4B), dimension(npl), intent(inout)      :: array_index1_child, array_index2_child

          end subroutine symba_casehitandrun
     end interface

     interface
          subroutine symba_casemerge (t, dt, index_enc, nmergeadd, nmergesub, mergeadd_list, mergesub_list, eoffset, vbs, & 
          encounter_file, out_type, npl, ntp, symba_pla, symba_tpa, nplplenc, npltpenc, pltpenc_list, plplenc_list, &
          array_index1_child, array_index2_child, m1, m2, rad1, rad2, x1, x2, v1, v2)
          use swiftest_globals
          use swiftest_data_structures
          use helio
          use symba
          implicit none
          integer(I4B), intent(in)                         :: index_enc
          integer(I4B), intent(inout)                      :: npl, ntp, nmergeadd, nmergesub, nplplenc, npltpenc
          real(DP), intent(in)                             :: t, dt
          real(DP), intent(inout)                          :: eoffset, m1, m2, rad1, rad2
          real(DP), dimension(ndim), intent(in)            :: vbs
          real(DP), dimension(ndim), intent(inout)         :: x1, x2, v1, v2
          character(*), intent(in)                         :: encounter_file, out_type
          type(symba_plplenc), intent(inout)               :: plplenc_list
          type(symba_pltpenc), intent(inout)               :: pltpenc_list
          type(symba_merger), intent(inout)                :: mergeadd_list, mergesub_list
          type(symba_pl), intent(inout)                    :: symba_pla
          type(symba_tp), intent(inout)                    :: symba_tpa
          integer(I4B), dimension(npl), intent(inout)      :: array_index1_child, array_index2_child

          end subroutine symba_casemerge
     end interface

     interface
          subroutine symba_caseresolve (t, dt, index_enc, nmergeadd, nmergesub, mergeadd_list, mergesub_list, eoffset, vbs, & 
          encounter_file, out_type, npl, ntp, symba_pla, symba_tpa, nplplenc, npltpenc, pltpenc_list, plplenc_list, regime, &
          nplmax, ntpmax, fragmax, mres, rres, array_index1_child, array_index2_child, m1, m2, rad1, rad2, x1, x2, v1, v2)
          use swiftest_globals
          use swiftest_data_structures
          use helio
          use symba
          implicit none
          integer(I4B), intent(in)                         :: index_enc, nplmax, ntpmax
          integer(I4B), intent(inout)                      :: npl, ntp, nmergeadd, nmergesub, nplplenc, npltpenc, fragmax
          real(DP), intent(in)                             :: t, dt
          real(DP), intent(inout)                          :: eoffset, m1, m2, rad1, rad2
          real(DP), dimension(3), intent(inout)            :: mres, rres
          real(DP), dimension(ndim), intent(in)            :: vbs
          real(DP), dimension(ndim), intent(inout)         :: x1, x2, v1, v2
          character(*), intent(in)                         :: encounter_file, out_type
          type(symba_plplenc), intent(inout)               :: plplenc_list
          type(symba_pltpenc), intent(inout)               :: pltpenc_list
          type(symba_merger), intent(inout)                :: mergeadd_list, mergesub_list
          type(symba_pl), intent(inout)                    :: symba_pla
          type(symba_tp), intent(inout)                    :: symba_tpa
          integer(I4B), intent(in)                         :: regime 
          integer(I4B), dimension(npl), intent(inout)      :: array_index1_child, array_index2_child

          end subroutine symba_caseresolve
     end interface

     interface
          subroutine symba_casesupercatastrophic (t, dt, index_enc, nmergeadd, nmergesub, mergeadd_list, mergesub_list, &
               eoffset, vbs, encounter_file, out_type, npl, ntp, symba_pla, symba_tpa, nplplenc, npltpenc, pltpenc_list, &
               plplenc_list, nplmax, ntpmax, fragmax, mres, rres, array_index1_child, array_index2_child, m1, m2, rad1, &
               rad2, x1, x2, v1, v2)
          use swiftest_globals
          use swiftest_data_structures
          use helio
          use symba
          implicit none
          integer(I4B), intent(in)                         :: index_enc, nplmax, ntpmax
          integer(I4B), intent(inout)                      :: npl, ntp, nmergeadd, nmergesub, nplplenc, npltpenc, fragmax
          real(DP), intent(in)                             :: t, dt
          real(DP), intent(inout)                          :: eoffset, m1, m2, rad1, rad2
          real(DP), dimension(3), intent(inout)            :: mres, rres
          real(DP), dimension(ndim), intent(in)            :: vbs
          real(DP), dimension(ndim), intent(inout)         :: x1, x2, v1, v2
          character(*), intent(in)                         :: encounter_file, out_type
          type(symba_plplenc), intent(inout)               :: plplenc_list
          type(symba_pltpenc), intent(inout)               :: pltpenc_list
          type(symba_merger), intent(inout)                :: mergeadd_list, mergesub_list
          type(symba_pl), intent(inout)                    :: symba_pla
          type(symba_tp), intent(inout)                    :: symba_tpa
          integer(I4B), dimension(npl), intent(inout)      :: array_index1_child, array_index2_child

          end subroutine symba_casesupercatastrophic
     end interface

     interface
          subroutine symba_chk(xr, vr, rhill1, rhill2, dt, irec, lencounter, lvdotr)
               use swiftest_globals
               use swiftest_data_structures
               use helio
               use symba
               implicit none
               logical(lgt), intent(out)          :: lencounter, lvdotr
               integer(I4B), intent(in)           :: irec
               real(DP), intent(in)               :: rhill1, rhill2, dt
               real(DP), dimension(:), intent(in) :: xr, vr
          end subroutine symba_chk
     end interface

     interface 
          subroutine symba_chk_eucl(num_encounters, k_plpl, symba_pla, dt, lencounter, lvdotr, nplplenc)
               use swiftest_globals
               use swiftest_data_structures
               use helio
               use symba
               implicit none
               type(symba_pl), intent(in)                    :: symba_pla
               integer(I4B), dimension(num_encounters), intent(out) :: lencounter, lvdotr
               integer(I4B), intent(in)           :: num_encounters
               integer(I4B), dimension(2,num_encounters),intent(in)   :: k_plpl
               real(DP), intent(in)               :: dt
               integer(I4B), intent(inout)        :: nplplenc
          end subroutine symba_chk_eucl
     end interface

     interface 
          subroutine symba_chk_eucl_pltp(num_encounters, k_pltp, symba_pla, symba_tpa, dt, lencounter, lvdotr, npltpenc)
               use swiftest_globals
               use swiftest_data_structures
               use helio
               use symba
               implicit none
               type(symba_pl), intent(in)                    :: symba_pla
               type(symba_tp), intent(in)                    :: symba_tpa
               integer(I4B), dimension(num_encounters), intent(out) :: lencounter, lvdotr
               integer(I4B), intent(in)           :: num_encounters
               integer(I4B), dimension(2,num_encounters),intent(in)   :: k_pltp
               real(DP), intent(in)               :: dt
               integer(I4B), intent(inout)        :: npltpenc
          end subroutine symba_chk_eucl_pltp
     end interface

     interface
          subroutine symba_discard_merge_pl(t, npl, symba_pla, nplplenc, plplenc_list)
               use swiftest_globals
               use swiftest_data_structures
               use helio
               use symba
               implicit none
               integer(I4B), intent(in)                      :: nplplenc
               integer(I4B), intent(inout)                   :: npl
               real(DP), intent(in)                          :: t
               type(symba_pl)                                :: symba_pla
               type(symba_plplenc), intent(in)               :: plplenc_list
          end subroutine symba_discard_merge_pl
     end interface

     interface
          subroutine symba_discard_peri_pl(t, npl, symba_pla, msys, qmin, qmin_alo, qmin_ahi, qmin_coord, ldiscards)
               use swiftest_globals
               use swiftest_data_structures
               use helio
               use symba
               implicit none
               logical(lgt), intent(inout) :: ldiscards
               integer(I4B), intent(in)    :: npl
               real(DP), intent(in)        :: t, msys, qmin, qmin_alo, qmin_ahi
               character(*), intent(in)    :: qmin_coord
               type(symba_pl), intent(inout)     :: symba_pla
          end subroutine symba_discard_peri_pl
     end interface

     interface
          subroutine symba_discard_pl(t, npl, nplmax, nsp, symba_pla, rmin, rmax, rmaxu, qmin, qmin_coord,          &
               qmin_alo, qmin_ahi, j2rp2, j4rp4, eoffset)
               use swiftest_globals
               use swiftest_data_structures
               use helio
               use symba
               implicit none
               integer(I4B), intent(in)    :: nplmax
               integer(I4B), intent(inout) :: npl, nsp
               real(DP), intent(in)        :: t, rmin, rmax, rmaxu, qmin, qmin_alo, qmin_ahi, j2rp2, j4rp4
               real(DP), intent(inout)     :: eoffset
               character(*), intent(in)    :: qmin_coord
               type(symba_pl), intent(inout)     :: symba_pla
          end subroutine symba_discard_pl
     end interface

     interface
          subroutine symba_discard_sun_pl(t, npl, msys, swiftest_pla, rmin, rmax, rmaxu, ldiscards)
               use swiftest_globals
               use swiftest_data_structures
               implicit none
               logical(lgt), intent(inout) :: ldiscards
               integer(I4B), intent(in)    :: npl
               real(DP), intent(in)        :: t, msys, rmin, rmax, rmaxu
               type(swiftest_pl), intent(inout)   :: swiftest_pla
          end subroutine symba_discard_sun_pl
     end interface

     interface
          subroutine symba_discard_tp(t, npl, ntp, nsp, symba_pla, symba_tpa, dt, &
               rmin, rmax, rmaxu, qmin, qmin_coord, qmin_alo, qmin_ahi, lclose, lrhill_present)
               use swiftest_globals
               use swiftest_data_structures
               use helio
               use symba
               implicit none
               logical(lgt), intent(in)    :: lclose, lrhill_present
               integer(I4B), intent(in)    :: npl
               integer(I4B), intent(inout) :: ntp, nsp
               real(DP), intent(in)        :: t, dt, rmin, rmax, rmaxu, qmin, qmin_alo, qmin_ahi
               character(*), intent(in)    :: qmin_coord
               type(symba_pl), intent(inout)     :: symba_pla
               type(symba_tp), intent(inout)     :: symba_tpa
          end subroutine symba_discard_tp
     end interface

     interface
          subroutine symba_energy(npl, nplmax, swiftest_pla, j2rp2, j4rp4, ke, pe, te, htot)
               use swiftest_globals
               use swiftest_data_structures
               implicit none
               integer(I4B), intent(in)               :: npl, nplmax
               real(DP), intent(in)                   :: j2rp2, j4rp4
               real(DP), intent(out)                  :: ke, pe, te
               real(DP), dimension(ndim), intent(out) :: htot
               type(swiftest_pl), intent(inout)             :: swiftest_pla
          end subroutine symba_energy
     end interface

     interface
          subroutine symba_fragmentation(t, dt, index_enc, nmergeadd, nmergesub, mergeadd_list, &
               mergesub_list, eoffset, vbs, encounter_file, out_type, npl, ntp, &
               symba_pla, symba_tpa, nplplenc, npltpenc, pltpenc_list, plplenc_list, &
               nplmax, ntpmax, fragmax)
               use swiftest_globals
               use swiftest_data_structures
               use helio
               use symba
               implicit none
               integer(I4B), intent(in)                         :: index_enc, nplmax, ntpmax
               integer(I4B), intent(inout)                      :: nmergeadd, nmergesub, nplplenc, npltpenc, fragmax
               integer(I4B), intent(inout)                      :: npl, ntp
               real(DP), intent(in)                             :: t, dt
               real(DP), intent(inout)                          :: eoffset
               real(DP), dimension(ndim), intent(in)            :: vbs
               character(*), intent(in)                         :: encounter_file, out_type
               type(symba_plplenc), intent(inout)               :: plplenc_list
               type(symba_pltpenc), intent(inout)               :: pltpenc_list
               type(symba_merger), intent(inout)                :: mergeadd_list, mergesub_list
               type(symba_pl), intent(inout)                    :: symba_pla
               type(symba_tp), intent(inout)                    :: symba_tpa
          end subroutine symba_fragmentation
     end interface

     interface
          subroutine symba_getacch(lextra_force, t, npl, nplm, nplmax, symba_pla, j2rp2, j4rp4, nplplenc, &
               plplenc_list)
               use swiftest_globals
               use swiftest_data_structures
               use helio
               use symba
               implicit none
               logical(lgt), intent(in)                      :: lextra_force
               integer(I4B), intent(in)                      :: npl, nplm, nplmax, nplplenc
               real(DP), intent(in)                          :: t, j2rp2, j4rp4
               type(symba_pl), intent(inout)                 :: symba_pla
               type(symba_plplenc), intent(in)               :: plplenc_list
          end subroutine symba_getacch
     end interface

     interface
          subroutine symba_getacch_tp(lextra_force, t, npl, nplm, nplmax, ntp, ntpmax, symba_pla, symba_tpa, &
               xh, j2rp2, j4rp4,  &
               npltpenc, pltpenc_list)
               use swiftest_globals
               use swiftest_data_structures
               use helio
               use symba
               implicit none
               logical(lgt), intent(in)                      :: lextra_force
               integer(I4B), intent(in)                      :: npl, nplm, nplmax, ntp, ntpmax, npltpenc
               real(DP), intent(in)                          :: t, j2rp2, j4rp4
               real(DP), dimension(ndim, npl), intent(in)    :: xh
               type(symba_pl), intent(inout)                 :: symba_pla
               type(symba_tp), intent(inout)                 :: symba_tpa
               type(symba_pltpenc), intent(in)               :: pltpenc_list
          end subroutine symba_getacch_tp
     end interface


     interface
          subroutine symba_getacch_eucl(lextra_force, t, npl, nplm, nplmax, symba_pla, j2rp2, j4rp4, nplplenc, &
               plplenc_list, num_plpl_comparisons, k_plpl)
               use swiftest_globals
               use symba
               implicit none
               logical(lgt), intent(in)                      :: lextra_force
               integer(I4B), intent(in)                      :: npl, nplm, nplmax, nplplenc, num_plpl_comparisons
               real(DP), intent(in)                          :: t, j2rp2, j4rp4
               type(symba_pl), intent(inout)                 :: symba_pla
               type(symba_plplenc), intent(in)               :: plplenc_list
               integer(I4B), dimension(num_plpl_comparisons,2),intent(in) :: k_plpl
          end subroutine symba_getacch_eucl
     end interface

     interface
          subroutine symba_getacch_tp_eucl(lextra_force, t, npl, nplm, nplmax, ntp, ntpmax, symba_pla, symba_tpa, &
               xh, j2rp2, j4rp4, npltpenc, pltpenc_list, num_pltp_comparisons, k_pltp)
               use swiftest_globals
               use symba
               implicit none
               logical(lgt), intent(in)                      :: lextra_force
               integer(I4B), intent(in)                      :: npl, nplm, nplmax, ntp, ntpmax, npltpenc, num_pltp_comparisons
               real(DP), intent(in)                          :: t, j2rp2, j4rp4
               real(DP), dimension(ndim, npl), intent(in)    :: xh
               type(symba_pl), intent(inout)                 :: symba_pla
               type(symba_tp), intent(inout)                 :: symba_tpa
               type(symba_pltpenc), intent(in)               :: pltpenc_list
               integer(I4B), dimension(num_pltp_comparisons,2), intent(in) :: k_pltp
          end subroutine symba_getacch_tp_eucl
     end interface

     interface
          subroutine symba_helio_drift(irec, npl, symba_pla, dt)
               use swiftest_globals
               use symba
               implicit none
               integer(I4B), intent(in)      :: irec, npl
               real(DP), intent(in)          :: dt
               type(symba_pl), intent(inout) :: symba_pla
          end subroutine symba_helio_drift
     end interface

     interface
          subroutine symba_helio_drift_tp(irec, ntp, symba_tpa, mu, dt)
               use swiftest_globals
               use symba
               implicit none
               integer(I4B), intent(in)      :: irec, ntp
               real(DP), intent(in)          :: mu, dt
               type(symba_tp), intent(inout) :: symba_tpa
          end subroutine symba_helio_drift_tp
     end interface

     interface
          subroutine symba_helio_getacch(lflag, lextra_force, t, npl, nplm, nplmax, helio_pla, j2rp2, j4rp4)
               use swiftest_globals
               use helio
               implicit none
               logical(lgt), intent(in)       :: lflag, lextra_force
               integer(I4B), intent(in)       :: npl, nplm, nplmax
               real(DP), intent(in)           :: t, j2rp2, j4rp4
               type(helio_pl), intent(inout)  :: helio_pla
          end subroutine symba_helio_getacch
     end interface

     interface
          subroutine symba_helio_getacch_int(npl, nplm, helio_pla)
               use swiftest_globals
               use helio
               implicit none
               integer(I4B), intent(in)      :: npl, nplm
               type(helio_pl), intent(inout) :: helio_pla
          end subroutine symba_helio_getacch_int
     end interface

     interface
          subroutine symba_kick(irec, nplplenc, npltpenc, plplenc_list, pltpenc_list, dt, sgn, symba_pla, &
               symba_tpa)
               use swiftest_globals
               use symba
               implicit none
               integer(I4B), intent(in)        :: irec, nplplenc, npltpenc
               real(DP), intent(in)            :: dt, sgn
               type(symba_plplenc), intent(in) :: plplenc_list
               type(symba_pltpenc), intent(in) :: pltpenc_list
               type(symba_pl), intent(inout)   :: symba_pla
               type(symba_tp), intent(inout)   :: symba_tpa
          end subroutine symba_kick
     end interface

     interface
          subroutine symba_merge_pl(t, dt, index_enc, nplplenc, plplenc_list, nmergeadd, nmergesub, &
               mergeadd_list, mergesub_list, eoffset, vbs, encounter_file, out_type, npl, symba_pla, &
               symba_tpa)
               use swiftest_globals
               use symba
               implicit none
               integer(I4B), intent(in)                         :: index_enc, nplplenc
               integer(I4B), intent(inout)                      :: nmergeadd, nmergesub, npl
               real(DP), intent(in)                             :: t, dt
               real(DP), intent(inout)                          :: eoffset
               real(DP), dimension(ndim), intent(in)            :: vbs
               character(*), intent(in)                         :: encounter_file, out_type
               type(symba_plplenc), intent(inout)               :: plplenc_list
               type(symba_merger),  intent(inout)               :: mergeadd_list, mergesub_list
               type(symba_pl), intent(inout)                    :: symba_pla
               type(symba_tp), intent(inout)                    :: symba_tpa
          end subroutine symba_merge_pl
     end interface

     interface
          subroutine symba_merge_tp(t, dt, index_enc, npltpenc, pltpenc_list, vbs, encounter_file, out_type, symba_pla, symba_tpa)
               use swiftest_globals
               use module_swifter
               use helio
               use symba
               implicit none
               integer(I4B), intent(in)               :: index_enc, npltpenc
               real(DP), intent(in)                   :: t, dt
               real(DP), dimension(ndim), intent(in)  :: vbs
               character(*), intent(in)               :: encounter_file, out_type
               type(symba_pltpenc), intent(inout)     :: pltpenc_list
               type(symba_pl), intent(inout)                    :: symba_pla
               type(symba_tp), intent(inout)                    :: symba_tpa
          end subroutine symba_merge_tp
     end interface

     interface
          subroutine symba_peri(lfirst, npl, symba_pla, msys, qmin_coord)
               use swiftest_globals
               use symba
               implicit none
               logical(lgt), intent(in)       :: lfirst
               integer(I4B), intent(in)       :: npl
               real(DP), intent(in)           :: msys
               character(*), intent(in)       :: qmin_coord
               type(symba_pl), intent(inout)  :: symba_pla
          end subroutine symba_peri
     end interface

     interface
          subroutine symba_rearray(t, npl, ntp, nsppl, nsptp, symba_pla, symba_tpa, nmergeadd, &
               mergeadd_list, discard_pla, discard_tpa, nplmax, j2rp2, j4rp4)
               use swiftest_globals
               use swiftest_data_structures
               use helio
               use symba
               implicit none
               integer(I4B), intent(inout)                      :: npl, ntp, nsppl, nsptp, nmergeadd, nplmax !change to fragadd
               real(DP), intent(in)                             :: t, j2rp2, j4rp4
               type(symba_pl), intent(inout)                    :: symba_pla
               type(symba_tp), intent(inout)                    :: symba_tpa
               type(swiftest_tp), intent(inout)                 :: discard_tpa
               type(swiftest_pl), intent(inout)                 :: discard_pla
               type(symba_merger), intent(inout)                :: mergeadd_list !change to fragadd_list

          end subroutine symba_rearray
     end interface  

     interface
          subroutine symba_reorder_pl(npl, symba_pla)
               use swiftest_globals
               use swiftest_data_structures
               use helio
               use symba
               implicit none
               integer(I4B), intent(in) :: npl
               type(symba_pl), intent(inout)  :: symba_pla
               integer(I4B)                              :: i
               integer(I4B), dimension(:), allocatable   :: index
               real(DP), dimension(:), allocatable       :: mass
               real(DP), dimension(:,:), allocatable     :: symba_plwkspa
               integer(I4B), dimension(:,:), allocatable :: symba_plwkspa_id_status
          end subroutine symba_reorder_pl
     end interface

     interface
          subroutine symba_setup(npl, ntp, symba_pla, symba_tpa, symba_pl1p, symba_tp1p, swiftest_pl1p, &
               swiftest_tp1p)
               use swiftest_globals
               use swiftest_data_structures
               use helio
               use symba
               implicit none
               integer(I4B), intent(in)                            :: npl, ntp
               type(swiftest_pl), pointer                           :: swiftest_pl1p
               type(swiftest_tp), pointer                           :: swiftest_tp1p
               type(symba_pl), intent(inout) :: symba_pla
               type(symba_tp), intent(inout) :: symba_tpa
               type(symba_pl), pointer                             :: symba_pl1p
               type(symba_tp), pointer                             :: symba_tp1p
          end subroutine symba_setup
     end interface

     interface
          subroutine symba_step_eucl(lfirst, lextra_force, lclose, t, npl, nplmax, ntp, ntpmax, symba_pla, symba_tpa, j2rp2, j4rp4,&
               dt,nplplenc, npltpenc, plplenc_list, pltpenc_list, nmergeadd, nmergesub, mergeadd_list, mergesub_list, eoffset,&
               mtiny,encounter_file, out_type, num_plpl_comparisons, k_plpl, num_pltp_comparisons, k_pltp)
               use swiftest_globals
               use swiftest_data_structures
               use helio
               use symba
               implicit none
               logical(lgt), intent(in)           :: lextra_force, lclose
               logical(lgt), intent(inout)        :: lfirst
               integer(I4B), intent(in)           :: npl, nplmax, ntp, ntpmax
               integer(I4B), intent(inout)        :: nplplenc, npltpenc, nmergeadd, nmergesub
               real(DP), intent(in)               :: t, j2rp2, j4rp4, dt, mtiny
               real(DP), intent(inout)            :: eoffset
               character(*), intent(in)           :: encounter_file, out_type
               type(symba_pl), intent(inout)      :: symba_pla
               type(symba_tp), intent(inout)      :: symba_tpa
               type(symba_plplenc), intent(inout) :: plplenc_list
               type(symba_pltpenc), intent(inout) :: pltpenc_list
               type(symba_merger), intent(inout)  :: mergeadd_list, mergesub_list
               integer(I4B), intent(in)           :: num_plpl_comparisons, num_pltp_comparisons
               integer(I4B), dimension(2,num_plpl_comparisons),intent(in) :: k_plpl
               integer(I4B), dimension(2,num_pltp_comparisons),intent(in) :: k_pltp
          end subroutine symba_step_eucl
     end interface

     interface
          subroutine symba_step(lfirst, lextra_force, lclose, t, npl, nplmax, ntp, ntpmax, symba_pla, &
               symba_tpa, j2rp2, j4rp4, dt, nplplenc, npltpenc, plplenc_list, pltpenc_list, nmergeadd, &
               nmergesub, mergeadd_list, mergesub_list, eoffset, mtiny, encounter_file, out_type, &
               fragmax)
               use swiftest_globals
               use swiftest_data_structures
               use helio
               use symba
               implicit none
               logical(lgt), intent(in)           :: lextra_force, lclose
               logical(lgt), intent(inout)        :: lfirst
               integer(I4B), intent(in)           :: npl, nplmax, ntp, ntpmax
               integer(I4B), intent(inout)        :: nplplenc, npltpenc, nmergeadd, nmergesub, fragmax
               real(DP), intent(in)               :: t, j2rp2, j4rp4, dt, mtiny
               real(DP), intent(inout)            :: eoffset
               character(*), intent(in)           :: encounter_file, out_type
               type(symba_pl), intent(inout)      :: symba_pla
               type(symba_tp), intent(inout)      :: symba_tpa
               type(symba_plplenc), intent(inout) :: plplenc_list
               type(symba_pltpenc), intent(inout) :: pltpenc_list
               type(symba_merger), intent(inout)  :: mergeadd_list, mergesub_list
          end subroutine symba_step
     end interface

! for testing purposes only _ use with symba_step_test
     interface
          subroutine symba_step_helio(lfirst, lextra_force, t, npl, nplm, nplmax, ntp, ntpmax, helio_pla, helio_tpa, j2rp2,     &
               j4rp4, dt)
               use swiftest_globals
               use swiftest_data_structures
               use helio
               use symba
               implicit none
               logical(lgt), intent(in)      :: lextra_force
               logical(lgt), intent(inout)   :: lfirst
               integer(I4B), intent(in)      :: npl, nplm, nplmax, ntp, ntpmax
               real(DP), intent(in)          :: t, j2rp2, j4rp4, dt
               type(helio_pl), intent(inout) :: helio_pla
               type(helio_tp), intent(inout) :: helio_tpa
          end subroutine symba_step_helio
     end interface

     interface
          subroutine symba_step_helio_pl(lfirst, lextra_force, t, npl, nplm, nplmax, helio_pla, j2rp2, j4rp4, dt, xbeg, xend,    &
               ptb, pte)
               use swiftest_globals
               use swiftest_data_structures
               use helio
               use symba
               implicit none
               logical(lgt), intent(in)                     :: lextra_force
               logical(lgt), intent(inout)                  :: lfirst
               integer(I4B), intent(in)                     :: npl, nplm, nplmax
               real(DP), intent(in)                         :: t, j2rp2, j4rp4, dt
               real(DP), dimension(ndim, nplm), intent(out) :: xbeg, xend
               real(DP), dimension(ndim), intent(out)       :: ptb, pte
               type(helio_pl), intent(inout)                :: helio_pla
          end subroutine symba_step_helio_pl
     end interface

     interface
          subroutine symba_step_interp_eucl(lextra_force, lclose, t, npl, nplm, nplmax, ntp, ntpmax, symba_pla, symba_tpa, j2rp2,&
               j4rp4, dt, eoffset, mtiny, nplplenc, npltpenc, plplenc_list, pltpenc_list, nmergeadd, nmergesub, mergeadd_list,&
               mergesub_list, encounter_file, out_type, num_plpl_comparisons, k_plpl, num_pltp_comparisons, k_pltp)
               use swiftest_globals
               use swiftest_data_structures
               use helio
               use symba
               implicit none
               logical(lgt), intent(in)           :: lextra_force, lclose
               integer(I4B), intent(in)           :: npl, nplm, nplmax, ntp, ntpmax, nplplenc, npltpenc, num_pltp_comparisons
               integer(I4B), intent(inout)        :: nmergeadd, nmergesub
               real(DP), intent(in)               :: t, j2rp2, j4rp4, dt, mtiny
               real(DP), intent(inout)            :: eoffset
               character(*), intent(in)           :: encounter_file, out_type
               type(symba_pl), intent(inout)      :: symba_pla
               type(symba_tp), intent(inout)      :: symba_tpa
               type(symba_plplenc), intent(inout) :: plplenc_list
               type(symba_pltpenc), intent(inout) :: pltpenc_list
               type(symba_merger), intent(inout)  :: mergeadd_list, mergesub_list
               integer(I4B), intent(in)                         :: num_plpl_comparisons
               integer(I4B), dimension(num_plpl_comparisons,2),intent(in) :: k_plpl
               integer(I4B), dimension(2,num_pltp_comparisons),intent(in) :: k_pltp
          end subroutine symba_step_interp_eucl
     end interface

     interface
          subroutine symba_step_interp(lextra_force, lclose, t, npl, nplm, nplmax, ntp, ntpmax, symba_pla, symba_tpa, j2rp2,    &
               j4rp4, dt, eoffset, mtiny, nplplenc, npltpenc, plplenc_list, pltpenc_list, nmergeadd, nmergesub, mergeadd_list,    &
               mergesub_list, encounter_file, out_type, fragmax)
               use swiftest_globals
               use swiftest_data_structures
               use helio
               use symba
               implicit none
               logical(lgt), intent(in)           :: lextra_force, lclose
               integer(I4B), intent(in)           :: npl, nplm, nplmax, ntp, ntpmax, nplplenc, npltpenc
               integer(I4B), intent(inout)        :: nmergeadd, nmergesub, fragmax
               real(DP), intent(in)               :: t, j2rp2, j4rp4, dt, mtiny
               real(DP), intent(inout)            :: eoffset
               character(*), intent(in)           :: encounter_file, out_type
               type(symba_pl), intent(inout)      :: symba_pla
               type(symba_tp), intent(inout)      :: symba_tpa
               type(symba_plplenc), intent(inout) :: plplenc_list
               type(symba_pltpenc), intent(inout) :: pltpenc_list
               type(symba_merger), intent(inout)  :: mergeadd_list, mergesub_list
          end subroutine symba_step_interp
     end interface


     interface
          recursive subroutine symba_step_recur(lclose, t, ireci, npl, nplm, ntp, symba_pla, symba_tpa, dt0, eoffset, nplplenc, &
               npltpenc, plplenc_list, pltpenc_list, nmergeadd, nmergesub, mergeadd_list, mergesub_list, encounter_file, & 
               out_type, nplmax, ntpmax, fragmax)
               use swiftest_globals
               use swiftest_data_structures
               use helio
               use symba
               implicit none
               logical(lgt), intent(in)           :: lclose
               integer(I4B), intent(in)           :: ireci, npl, nplm, ntp, nplplenc, npltpenc, nplmax, ntpmax, fragmax
               integer(I4B), intent(inout)        :: nmergeadd, nmergesub
               real(DP), intent(in)               :: t, dt0
               real(DP), intent(inout)            :: eoffset
               character(*), intent(in)           :: encounter_file, out_type
               type(symba_pl), intent(inout)      :: symba_pla
               type(symba_tp), intent(inout)      :: symba_tpa
               type(symba_plplenc), intent(inout) :: plplenc_list
               type(symba_pltpenc), intent(inout) :: pltpenc_list
               type(symba_merger), intent(inout)  :: mergeadd_list, mergesub_list
          end subroutine symba_step_recur
     end interface

     interface
          subroutine symba_user_getacch(t, npl, symba_pla)
               use swiftest_globals
               use symba
               implicit none
               integer(I4B), intent(in)     :: npl
               real(DP), intent(in)         :: t
               type(symba_pl), intent(inout):: symba_pla
          end subroutine symba_user_getacch
     end interface

     interface
          subroutine symba_user_getacch_tp(t, ntp, symba_tpa)
               use swiftest_globals
               use symba
               implicit none
               integer(I4B), intent(in)       :: ntp
               real(DP), intent(in)           :: t
               type(symba_tp), intent(inout)  :: symba_tpa
          end subroutine symba_user_getacch_tp
     end interface

     interface
          subroutine util_exit(code)
               use swiftest_globals
               implicit none
               integer(I4B), intent(in) :: code
          end subroutine util_exit
     end interface

     interface
          subroutine util_dist_index_plpl(npl, nplm, num_comparisons, k_plpl)
               use swiftest_globals
               use swiftest_data_structures
               implicit none
               integer(I4B), intent(in)  :: npl, nplm
               integer(I4B), dimension(:,:),allocatable,intent(out) :: k_plpl
               integer(I4B), intent(out) :: num_comparisons
          end subroutine
     end interface

     interface
          subroutine util_dist_index_pltp(nplm, ntp, num_comparisons, k_pltp)
               use swiftest_globals
               use swiftest_data_structures
               implicit none
               integer(I4B), intent(in)  :: nplm, ntp
               integer(I4B), dimension(:,:),allocatable,intent(out) :: k_pltp
               integer(I4B), intent(out) :: num_comparisons
          end subroutine util_dist_index_pltp
     end interface

     interface
          subroutine util_dist_eucl_plpl(npl, invar, num_comparisons, k_plpl, outvar)
               use swiftest_globals
               use swiftest_data_structures
               use symba
               implicit none
               integer(I4B), intent(in)  :: npl
               integer(I4B), dimension(2,num_comparisons),intent(in) :: k_plpl
               integer(I4B), intent(in) :: num_comparisons
               real(DP),dimension(ndim,npl),intent(in) :: invar
               real(DP), dimension(ndim,num_comparisons),intent(inout) :: outvar
          end subroutine util_dist_eucl_plpl
     end interface

     interface
          subroutine util_dist_eucl_pltp(npl, ntp, planets, test_particles, num_pltp_comparisons, k_pltp, outvar)
               use swiftest_globals
               use swiftest_data_structures
               use symba
               implicit none
               integer(I4B), intent(in) :: npl, ntp
               integer(I4B), dimension(num_pltp_comparisons,2),intent(in) :: k_pltp
               integer(I4B), intent(in) :: num_pltp_comparisons
               real(DP),dimension(ndim,npl),intent(in) :: planets
               real(DP),dimension(ndim,ntp),intent(in) :: test_particles
               real(DP), dimension(ndim,num_pltp_comparisons),intent(inout) :: outvar
          end subroutine
     end interface

     interface
          subroutine util_hills(npl, swiftest_pla)
               use swiftest_globals
               use swiftest_data_structures
               implicit none
               integer(I4B), intent(in)         :: npl
               type(swiftest_pl), intent(inout) :: swiftest_pla
          end subroutine util_hills
     end interface

     interface
          subroutine util_index(arr, index)
               use swiftest_globals
               use module_nrutil
               implicit none
               integer(I4B), dimension(:), intent(out) :: index
               real(DP), dimension(:), intent(in)      :: arr
          end subroutine util_index
     end interface

     interface
          subroutine util_peri(lfirst, ntp, swiftest_tpa, mu, msys, qmin_coord)
               use swiftest_globals
               use swiftest_data_structures
               implicit none
               logical(lgt), intent(in)         :: lfirst
               integer(I4B), intent(in)         :: ntp
               real(DP), intent(in)             :: mu, msys
               character(*), intent(in)         :: qmin_coord
               type(swiftest_tp), intent(inout) :: swiftest_tpa
          end subroutine util_peri
     end interface

     interface 
          subroutine util_resize_pl(symba_pla, npl_new, npl_old)
               use swiftest_globals
               use symba
               use swiftest_data_structures
               use helio
               use module_nrutil
               implicit none
               type(symba_pl), intent(inout) :: symba_pla
               integer(I4B), intent(in)      :: npl_old, npl_new
          end subroutine util_resize_pl
     end interface

     interface util_sort
          subroutine util_sort_I4B(arr)
               use swiftest_globals
               implicit none
               integer(I4B), dimension(:), intent(inout) :: arr
          end subroutine util_sort_I4B
          subroutine util_sort_sp(arr)
               use swiftest_globals
               implicit none
               real(sp), dimension(:), intent(inout) :: arr
          end subroutine util_sort_sp
          subroutine util_sort_DP(arr)
               use swiftest_globals
               implicit none
               real(DP), dimension(:), intent(inout) :: arr
          end subroutine util_sort_DP
     end interface

     interface
          subroutine util_toupper(string)
               use swiftest_globals
               implicit none
               character(*), intent(inout) :: string
          end subroutine util_toupper
     end interface

     interface
          subroutine util_valid(npl, ntp, swiftest_pla, swiftest_tpa)
               use swiftest_globals
               use swiftest_data_structures
               implicit none
               integer(I4B), intent(in)         :: npl, ntp
               type(swiftest_pl), intent(inout) :: swiftest_pla
               type(swiftest_tp), intent(inout) :: swiftest_tpa
          end subroutine util_valid
     end interface

     interface
          subroutine util_version
               use swiftest_globals
               implicit none
          end subroutine util_version
     end interface

     interface
         function util_kahan_sum(xsum_current, xi, xerror) 
            use swiftest_globals
            implicit none
            real(DP)                :: util_kahan_sum
            real(DP), intent(in)    :: xsum_current, xi
            real(DP), intent(inout) :: xerror
         end function
     end interface

end module module_interfaces