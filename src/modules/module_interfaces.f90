module module_interfaces
  implicit none
  
   interface 
      subroutine coord_h2b(npl, swiftest_plA, msys)
         use swiftest_globals
         use swiftest_data_structures
         implicit none
         integer(I4B), intent(in)    :: npl
         real(DP), intent(out)     :: msys
         type(swiftest_pl), intent(inout) :: swiftest_plA
      end subroutine coord_h2b

      subroutine coord_h2b_tp(ntp, swiftest_tpA, swiftest_plA)
         use swiftest_globals
         use swiftest_data_structures
         implicit none
         integer(I4B), intent(in)    :: ntp
         type(swiftest_tp), intent(inout) :: swiftest_tpA
         type(swiftest_pl), intent(inout) :: swiftest_plA
      end subroutine coord_h2b_tp

      subroutine coord_vb2vh(npl, swiftest_plA)
         use swiftest_globals
         use swiftest_data_structures
         implicit none
         integer(I4B), intent(in)    :: npl
         type(swiftest_pl), intent(inout) :: swiftest_plA
      end subroutine coord_vb2vh

      subroutine coord_vb2vh_tp(ntp, swiftest_tpA, vs)
         use swiftest_globals
         use swiftest_data_structures
         implicit none
         integer(I4B), intent(in)     :: ntp
         real(DP), dimension(ndim), intent(in) :: vs
         type(swiftest_tp), intent(inout)   :: swiftest_tpA
      end subroutine coord_vb2vh_tp

      subroutine coord_vh2vb(npl, swiftest_plA, msys)
         use swiftest_globals
         use swiftest_data_structures
         implicit none
         integer(I4B), intent(in)   :: npl
         real(DP), intent(out)    :: msys
         type(swiftest_pl),intent(inout) :: swiftest_plA
      end subroutine coord_vh2vb

      subroutine coord_vh2vb_tp(ntp, swiftest_tpA, vs)
         use swiftest_globals
         use swiftest_data_structures
         implicit none
         integer(I4B), intent(in)     :: ntp
         real(DP), dimension(ndim), intent(in) :: vs
         type(swiftest_tp), intent(inout)   :: swiftest_tpA
      end subroutine coord_vh2vb_tp

      subroutine discard(t, dt, npl, ntp, swiftest_plA, swiftest_tpA, rmin, rmax, rmaxu, qmin,  &
         qmin_alo, qmin_ahi, qmin_coord, lclose, lrhill_present)
         use swiftest_globals
         use swiftest_data_structures
         implicit none
         logical(lgt), intent(in)  :: lclose, lrhill_present
         integer(I4B), intent(in)  :: npl, ntp
         real(DP), intent(in)   :: t, dt, rmin, rmax, rmaxu, qmin, qmin_alo, qmin_ahi
         character(*), intent(in)  :: qmin_coord
         type(swiftest_pl), intent(inout) :: swiftest_plA
         type(swiftest_tp), intent(inout) :: swiftest_tpA
      end subroutine discard

      subroutine discard_peri(t, npl, ntp, swiftest_plA, swiftest_tpA, msys, qmin, qmin_alo, & 
         qmin_ahi, qmin_coord, lrhill_present)
         use swiftest_globals
         use swiftest_data_structures
         implicit none
         logical(lgt), intent(in)  :: lrhill_present
         integer(I4B), intent(in)  :: npl, ntp
         real(DP), intent(in)   :: t, msys, qmin, qmin_alo, qmin_ahi
         character(*), intent(in)  :: qmin_coord
         type(swiftest_pl), intent(inout) :: swiftest_plA
         type(swiftest_tp), intent(inout) :: swiftest_tpA
      end subroutine discard_peri

      subroutine discard_pl_close(dx, dv, dt, r2crit, iflag, r2min)
         use swiftest_globals
         implicit none
         integer(I4B), intent(out)      :: iflag
         real(DP), intent(in)      :: dt, r2crit
         real(DP), dimension(ndim), intent(in) :: dx, dv
         real(DP), intent(out)      :: r2min
      end subroutine discard_pl_close

      subroutine discard_pl(t, dt, npl, ntp, swiftest_plA, swiftest_tpA)
         use swiftest_globals
         use swiftest_data_structures
         implicit none
         integer(I4B), intent(in)  :: npl, ntp
         real(DP), intent(in)   :: t, dt
         type(swiftest_pl), intent(inout) :: swiftest_plA
         type(swiftest_tp), intent(inout) :: swiftest_tpA
      end subroutine discard_pl

      subroutine discard_sun(t, ntp, msys, swifter_tpA, rmin, rmax, rmaxu)
         use swiftest_globals
         use swiftest_data_structures
         implicit none
         integer(I4B), intent(in)  :: ntp
         real(DP), intent(in)   :: t, msys, rmin, rmax, rmaxu
         type(swiftest_tp), intent(inout) :: swifter_tpA
      end subroutine discard_sun

      subroutine io_discard_write_symba(t, mtiny, npl, ntp, nsppl, nsptp, nmergeadd, nmergesub, symba_plA, & 
         discard_plA, discard_tpA, mergeadd_list, mergesub_list, fname, lbig_discard)
         use swiftest_globals
         use swiftest_data_structures
         use symba
         implicit none
         logical(lgt), intent(in)       :: lbig_discard
         integer(I4B), intent(in)       :: npl, ntp, nsppl, nsptp, nmergeadd, nmergesub
         real(DP), intent(in)         :: t, mtiny
         character(*), intent(in)       :: fname
         type(symba_pl), intent(inout)      :: symba_plA
         type(swiftest_tp), intent(inout)      :: discard_tpA
         type(swiftest_pl), intent(inout)      :: discard_plA
         type(symba_merger), intent(inout)     :: mergeadd_list, mergesub_list
      end subroutine io_discard_write_symba

      subroutine io_dump_param(nplmax, ntpmax, ntp, t, tstop, dt, in_type, istep_out, outfile, out_type, out_form,     &
         istep_dump, j2rp2, j4rp4, lclose, rmin, rmax, rmaxu, qmin, qmin_coord, qmin_alo, qmin_ahi, encounter_file,    &
         lextra_force, lbig_discard, lrhill_present, mtiny, lpython)
         use swiftest_globals
         implicit none
         logical(lgt), intent(in) :: lclose, lextra_force, lbig_discard, lrhill_present, lpython
         integer(I4B), intent(in) :: nplmax, ntpmax, ntp, istep_out, istep_dump
         real(DP), intent(in)   :: t, tstop, dt, j2rp2, j4rp4, rmin, rmax, rmaxu, qmin, qmin_alo, qmin_ahi, mtiny
         character(*), intent(in) :: qmin_coord, encounter_file, in_type, outfile, out_type, out_form
      end subroutine io_dump_param

      subroutine io_dump_pl(npl, swiftest_plA, lclose, lrhill_present)
         use swiftest_globals
         use swiftest_data_structures
         implicit none
         logical(lgt), intent(in)   :: lclose, lrhill_present
         integer(I4B), intent(in)   :: npl
         type(swiftest_pl), intent(inout):: swiftest_plA
      end subroutine io_dump_pl

      subroutine io_dump_tp(ntp, swiftest_tpA)
         use swiftest_globals
         use swiftest_data_structures
         implicit none
         integer(I4B), intent(in)   :: ntp
         type(swiftest_tp), intent(inout):: swiftest_tpA
      end subroutine io_dump_tp

      subroutine io_getn(inplfile, intpfile, in_type, npl, nplmax, ntp, ntpmax)
         use swiftest_globals
         implicit none
         integer(I4B), intent(inout) :: nplmax, ntpmax
         integer(I4B), intent(out)   :: npl, ntp
         character(*), intent(in)   :: inplfile, intpfile, in_type
      end subroutine io_getn

      subroutine io_get_token(buffer, ilength, ifirst, ilast, ierr)
         use swiftest_globals
         implicit none
         integer(I4B), intent(in)   :: ilength
         integer(I4B), intent(inout) :: ifirst
         integer(I4B), intent(out)   :: ilast, ierr
         character(*), intent(in)   :: buffer
      end subroutine io_get_token

      subroutine io_init_param(inparfile, nplmax, ntpmax, t0, tstop, dt, inplfile, intpfile, in_type, istep_out, outfile,   &
         out_type, out_form, out_stat, istep_dump, j2rp2, j4rp4, lclose, rmin, rmax, rmaxu, qmin, qmin_coord, qmin_alo,   &
         qmin_ahi, encounter_file, lextra_force, lbig_discard, lrhill_present, mtiny, lpython, lenergy)
         use swiftest_globals
         implicit none
         logical(lgt), intent(out) :: lclose, lextra_force, lbig_discard, lrhill_present, lpython, lenergy
         integer(I4B), intent(out) :: nplmax, ntpmax, istep_out, istep_dump
         real(DP), intent(out)   :: t0, tstop, dt, j2rp2, j4rp4, rmin, rmax, rmaxu, qmin, qmin_alo, qmin_ahi
         real(DP), intent(out), optional :: mtiny
         character(*), intent(in)  :: inparfile
         character(*), intent(out) :: qmin_coord, encounter_file, inplfile, intpfile, in_type, outfile, out_type, out_form, &
               out_stat
      end subroutine io_init_param

      subroutine io_init_pl(inplfile, in_type, lclose, lrhill_present, npl, symba_plA)
         use swiftest_globals
         use swiftest_data_structures
         use symba
         use helio
         implicit none
         logical(lgt), intent(in)    :: lclose, lrhill_present
         integer(I4B), intent(in)    :: npl
         character(*), intent(in)    :: inplfile, in_type
         type(symba_pl), intent(inout)   :: symba_plA
      end subroutine io_init_pl

      subroutine io_init_tp(intpfile, in_type, ntp, symba_tpA)
         use swiftest_globals
         use swiftest_data_structures
         use symba
         use helio
         implicit none
         integer(I4B), intent(in)    :: ntp
         character(*), intent(in)    :: intpfile, in_type
         type(symba_tp), intent(inout)   :: symba_tpA
      end subroutine io_init_tp

      subroutine io_open(iu, fname, fopenstat, fmt, ierr)
         use swiftest_globals
         implicit none
         integer(I4B), intent(in)  :: iu
         integer(I4B), intent(out) :: ierr
         character(*), intent(in)  :: fname, fopenstat, fmt
      end subroutine io_open

      function io_read_encounter(t, name1, name2, mass1, mass2, xh1, xh2, vh1, vh2, encounter_file, out_type)
         use swiftest_globals
         implicit none
         integer(I4B)         :: io_read_encounter
         integer(I4B), intent(out)     :: name1, name2
         real(DP), intent(out)      :: t, mass1, mass2
         real(DP), dimension(ndim), intent(out) :: xh1, xh2, vh1, vh2
         character(*), intent(in)      :: encounter_file,out_type
      end function io_read_encounter

      function io_read_hdr(iu, t, npl, ntp, iout_form, out_type)
         use swiftest_globals
         implicit none
         integer(I4B)      :: io_read_hdr
         integer(I4B), intent(in)   :: iu
         integer(I4B), intent(out)  :: npl, ntp, iout_form
         real(DP), intent(out)   :: t
         character(*), intent(in)   :: out_type
      end function io_read_hdr

      function io_read_line(iu, name, d1, d2, d3, d4, d5, d6, out_type, mass, radius)
         use swiftest_globals
         implicit none
         integer(I4B)      :: io_read_line
         integer(I4B), intent(in)   :: iu
         integer(I4B), intent(out)   :: name
         real(DP), intent(out)    :: d1, d2, d3, d4, d5, d6
         real(DP), optional, intent(out) :: mass, radius
         character(*), intent(in)   :: out_type
      end function io_read_line

      subroutine io_write_encounter(t, name1, name2, mass1, mass2, radius1, radius2, &
         xh1, xh2, vh1, vh2, encounter_file, out_type)
         use swiftest_globals
         implicit none
         integer(I4B), intent(in)     :: name1, name2
         real(DP), intent(in)      :: t, mass1, mass2, radius1, radius2
         real(DP), dimension(ndim), intent(in) :: xh1, xh2, vh1, vh2
         character(*), intent(in)     :: encounter_file, out_type
      end subroutine io_write_encounter

      subroutine io_write_frame(t, npl, ntp, swiftest_plA, swiftest_tpA, outfile, &
         out_type, out_form, out_stat)
         use swiftest_globals
         use swiftest_data_structures
         implicit none
         integer(I4B), intent(in)  :: npl, ntp
         real(DP), intent(in)   :: t
         character(*), intent(in)  :: outfile, out_type, out_form, out_stat
         type(swiftest_pl), intent(inout) :: swiftest_plA
         type(swiftest_tp), intent(inout) :: swiftest_tpA
      end subroutine io_write_frame

      subroutine io_write_hdr(iu, t, npl, ntp, iout_form, out_type)
         use swiftest_globals
         implicit none
         integer(I4B), intent(in) :: iu, npl, ntp, iout_form
         real(DP), intent(in)   :: t
         character(*), intent(in) :: out_type
      end subroutine io_write_hdr

      subroutine io_write_line(iu, name, d1, d2, d3, d4, d5, d6, out_type, mass, radius)
         use swiftest_globals
         implicit none
         integer(I4B), intent(in)   :: iu, name
         real(DP), intent(in)    :: d1, d2, d3, d4, d5, d6
         real(DP), optional, intent(in) :: mass, radius
         character(*), intent(in)   :: out_type
      end subroutine io_write_line

      subroutine obl_acc(npl, swiftest_plA, j2rp2, j4rp4, xh, irh, aobl)
         use swiftest_globals
         use swiftest_data_structures
         implicit none
         integer(I4B), intent(in)      :: npl
         real(DP), intent(in)        :: j2rp2, j4rp4
         real(DP), dimension(npl), intent(in)   :: irh
         real(DP), dimension(ndim, npl), intent(in)  :: xh
         real(DP), dimension(ndim, npl), intent(out) :: aobl
         type(swiftest_pl), intent(inout)     :: swiftest_plA
      end subroutine obl_acc

      subroutine obl_acc_tp(ntp, xht, j2rp2, j4rp4, irht, aoblt, msun)
         use swiftest_globals
         implicit none
         integer(I4B), intent(in)      :: ntp
         real(DP), intent(in)        :: j2rp2, j4rp4, msun
         real(DP), dimension(ntp), intent(in)   :: irht
         real(DP), dimension(ndim, ntp), intent(in)  :: xht
         real(DP), dimension(ndim, ntp), intent(out) :: aoblt
      end subroutine obl_acc_tp

      subroutine obl_pot(npl, swiftest_plA, j2rp2, j4rp4, xh, irh, oblpot)
         use swiftest_globals
         use swiftest_data_structures
         implicit none
         integer(I4B), intent(in)       :: npl
         real(DP), intent(in)       :: j2rp2, j4rp4
         real(DP), intent(out)        :: oblpot
         real(DP), dimension(npl), intent(in)   :: irh
         real(DP), dimension(ndim, npl), intent(in) :: xh
         type(swiftest_pl), intent(inout)    :: swiftest_plA
      end subroutine obl_pot

      pure subroutine orbel_scget(angle, sx, cx)
         use swiftest_globals
         implicit none
         real(DP), intent(in)  :: angle
         real(DP), intent(out) :: sx, cx
      end subroutine orbel_scget

      subroutine orbel_xv2aeq(x, v, mu, a, e, q)
         use swiftest_globals
         implicit none
         real(DP), intent(in)      :: mu
         real(DP), dimension(ndim), intent(in) :: x, v
         real(DP), intent(out)      :: a, e, q
      end subroutine orbel_xv2aeq

      subroutine orbel_xv2aqt(x, v, mu, a, q, capm, tperi)
         use swiftest_globals
         implicit none
         real(DP), intent(in)      :: mu
         real(DP), dimension(ndim), intent(in) :: x, v
         real(DP), intent(out)      :: a, q, capm, tperi
      end subroutine orbel_xv2aqt

      subroutine orbel_xv2el(x, v, mu, a, e, inc, capom, omega, capm)
         use swiftest_globals
         implicit none
         real(DP), intent(in)      :: mu
         real(DP), dimension(ndim), intent(in) :: x, v
         real(DP), intent(out)      :: a, e, inc, capom, omega, capm
      end subroutine orbel_xv2el

      subroutine python_io_write_frame_pl(t, symba_plA, npl, out_stat)
         use swiftest_globals
         use swiftest_data_structures
         use helio
         use symba
         implicit none
         real(DP), intent(in)   :: t
         type(symba_pl),intent(in) :: symba_plA
         integer, intent(in)   :: npl
         character(*), intent(in)  :: out_stat
      end subroutine python_io_write_frame_pl

      subroutine python_io_write_frame_tp(t, symba_tpA, ntp, out_stat)
         use swiftest_globals
         use swiftest_data_structures
         use helio
         use symba
         implicit none
         real(DP), intent(in)   :: t
         type(symba_tp),intent(in) :: symba_tpA
         integer, intent(in)   :: ntp
         character(*), intent(in)  :: out_stat
      end subroutine python_io_write_frame_tp

      subroutine rmvs_chk_ind(xr, vr, dt, r2crit, iflag)
         use swiftest_globals
         implicit none
         real(DP), intent(in)      :: dt, r2crit
         real(DP), dimension(ndim), intent(in) :: xr, vr
         integer(I4B), intent(out)      :: iflag
      end subroutine rmvs_chk_ind

      subroutine symba_casedisruption (t, dt, index_enc, nmergeadd, nmergesub, mergeadd_list, mergesub_list, eoffset, vbs, & 
         symba_plA, nplplenc, plplenc_list, nplmax, ntpmax, fragmax, mres, rres, m1, m2, rad1, rad2, x1, x2, v1, v2)
         use swiftest_globals
         use symba
         implicit none
         integer(I4B), intent(in)       :: index_enc, nplmax, ntpmax
         integer(I4B), intent(inout)        :: nmergeadd, nmergesub, nplplenc, fragmax
         real(DP), intent(in)          :: t, dt
         real(DP), intent(inout)        :: eoffset, m1, m2, rad1, rad2
         real(DP), dimension(:), intent(inout)     :: mres, rres
         real(DP), dimension(:), intent(in)     :: vbs
         real(DP), dimension(:), intent(inout)    :: x1, x2, v1, v2
         type(symba_plplenc), intent(inout)      :: plplenc_list
         type(symba_merger), intent(inout)     :: mergeadd_list, mergesub_list
         type(symba_pl), intent(inout)      :: symba_plA

      end subroutine symba_casedisruption

      subroutine symba_casehitandrun (t, dt, index_enc, nmergeadd, nmergesub, mergeadd_list, mergesub_list, eoffset, vbs, & 
         symba_plA, nplplenc, plplenc_list, nplmax, ntpmax, fragmax, mres, rres, m1, m2, rad1, rad2, x1, x2, v1, v2)
         use swiftest_globals
         use swiftest_data_structures
         use helio
         use symba
         implicit none
         integer(I4B), intent(in)       :: index_enc, nplmax, ntpmax
         integer(I4B), intent(inout)        :: nmergeadd, nmergesub, nplplenc, fragmax
         real(DP), intent(in)          :: t, dt
         real(DP), intent(inout)        :: eoffset, m1, m2, rad1, rad2
         real(DP), dimension(:), intent(inout)     :: mres, rres
         real(DP), dimension(:), intent(in)     :: vbs
         real(DP), dimension(:), intent(inout)    :: x1, x2, v1, v2
         type(symba_plplenc), intent(inout)      :: plplenc_list
         type(symba_merger), intent(inout)     :: mergeadd_list, mergesub_list
         type(symba_pl), intent(inout)      :: symba_plA

      end subroutine symba_casehitandrun

      subroutine symba_casemerge (t, index_enc, nmergeadd, nmergesub, mergeadd_list, mergesub_list, eoffset, vbs, npl, &
         symba_plA, nplplenc, plplenc_list, array_index1_child, array_index2_child, m1, m2, rad1, rad2, x1, x2, v1, v2)
         use swiftest_globals
         use swiftest_data_structures
         use helio
         use symba
         implicit none
         integer(I4B), intent(in)       :: index_enc
         integer(I4B), intent(inout)        :: npl, nmergeadd, nmergesub, nplplenc
         real(DP), intent(in)          :: t
         real(DP), intent(inout)        :: eoffset, m1, m2, rad1, rad2
         real(DP), dimension(:), intent(in)     :: vbs
         real(DP), dimension(:), intent(inout)    :: x1, x2, v1, v2
         type(symba_plplenc), intent(inout)      :: plplenc_list
         type(symba_merger), intent(inout)     :: mergeadd_list, mergesub_list
         type(symba_pl), intent(inout)      :: symba_plA
         integer(I4B), dimension(npl), intent(inout)   :: array_index1_child, array_index2_child
      end subroutine symba_casemerge

      subroutine symba_caseresolve (t, dt, index_enc, nmergeadd, nmergesub, mergeadd_list, mergesub_list, eoffset, vbs, & 
         npl, symba_plA, nplplenc, plplenc_list, regime, nplmax, ntpmax, fragmax, mres, rres, array_index1_child, &
         array_index2_child, m1, m2, rad1, rad2, x1, x2, v1, v2)
         use swiftest_globals
         use swiftest_data_structures
         use helio
         use symba
         implicit none
         integer(I4B), intent(in)       :: index_enc, nplmax, ntpmax
         integer(I4B), intent(inout)        :: npl, nmergeadd, nmergesub, nplplenc, fragmax
         real(DP), intent(in)          :: t, dt
         real(DP), intent(inout)        :: eoffset, m1, m2, rad1, rad2
         real(DP), dimension(:), intent(inout)     :: mres, rres
         real(DP), dimension(:), intent(in)     :: vbs
         real(DP), dimension(:), intent(inout)    :: x1, x2, v1, v2
         type(symba_plplenc), intent(inout)      :: plplenc_list
         type(symba_merger), intent(inout)     :: mergeadd_list, mergesub_list
         type(symba_pl), intent(inout)      :: symba_plA
         integer(I4B), intent(in)       :: regime
         integer(I4B), dimension(npl), intent(inout)   :: array_index1_child, array_index2_child
      end subroutine symba_caseresolve

      subroutine symba_casesupercatastrophic (t, dt, index_enc, nmergeadd, nmergesub, mergeadd_list, mergesub_list, eoffset, vbs, & 
         symba_plA, nplplenc, plplenc_list, nplmax, ntpmax, fragmax, mres, rres, m1, m2, rad1, rad2, x1, x2, v1, v2)
         use swiftest_globals
         use swiftest_globals
         use swiftest_data_structures
         use helio
         use symba
         implicit none
         integer(I4B), intent(in)       :: index_enc, nplmax, ntpmax
         integer(I4B), intent(inout)        :: nmergeadd, nmergesub, nplplenc, fragmax
         real(DP), intent(in)          :: t, dt
         real(DP), intent(inout)        :: eoffset, m1, m2, rad1, rad2
         real(DP), dimension(:), intent(inout)     :: mres, rres
         real(DP), dimension(:), intent(in)     :: vbs
         real(DP), dimension(:), intent(inout)    :: x1, x2, v1, v2
         type(symba_plplenc), intent(inout)      :: plplenc_list
         type(symba_merger), intent(inout)     :: mergeadd_list, mergesub_list
         TYPE(symba_pl), INTENT(INOUT)      :: symba_plA
      end subroutine symba_casesupercatastrophic

      subroutine symba_chk(xr, vr, rhill1, rhill2, dt, irec, lencounter, lvdotr)
         use swiftest_globals
         use swiftest_data_structures
         use helio
         use symba
         implicit none
         logical(lgt), intent(out)   :: lencounter, lvdotr
         integer(I4B), intent(in)    :: irec
         real(DP), intent(in)      :: rhill1, rhill2, dt
         real(DP), dimension(:), intent(in) :: xr, vr
      end subroutine symba_chk

      subroutine symba_chk_eucl(num_encounters, k_plpl, symba_plA, dt, lencounter, lvdotr, nplplenc)
         use swiftest_globals
         use swiftest_data_structures
         use helio
         use symba
         implicit none
         type(symba_pl), intent(in)      :: symba_plA
         integer(I4B), dimension(num_encounters), intent(out) :: lencounter, lvdotr
         integer(I4B), intent(in)    :: num_encounters
         integer(I4B), dimension(2,num_encounters),intent(in)   :: k_plpl
         real(DP), intent(in)      :: dt
         integer(I4B), intent(inout)   :: nplplenc
      end subroutine symba_chk_eucl

      subroutine symba_chk_eucl_pltp(num_encounters, k_pltp, symba_plA, symba_tpA, dt, lencounter, lvdotr, npltpenc)
         use swiftest_globals
         use swiftest_data_structures
         use helio
         use symba
         implicit none
         type(symba_pl), intent(in)      :: symba_plA
         type(symba_tp), intent(in)      :: symba_tpA
         integer(I4B), dimension(num_encounters), intent(out) :: lencounter, lvdotr
         integer(I4B), intent(in)    :: num_encounters
         integer(I4B), dimension(2,num_encounters),intent(in)   :: k_pltp
         real(DP), intent(in)      :: dt
         integer(I4B), intent(inout)   :: npltpenc
      end subroutine symba_chk_eucl_pltp

      subroutine symba_discard_merge_pl(t, npl, symba_plA, nplplenc, plplenc_list)
         use swiftest_globals
         use swiftest_data_structures
         use helio
         use symba
         implicit none
         integer(I4B), intent(in)        :: nplplenc
         integer(I4B), intent(inout)       :: npl
         real(DP), intent(in)        :: t
         type(symba_pl)         :: symba_plA
         type(symba_plplenc), intent(in)      :: plplenc_list
      end subroutine symba_discard_merge_pl

      subroutine symba_discard_peri_pl(t, npl, symba_plA, msys, qmin, qmin_alo, qmin_ahi, qmin_coord, ldiscards)
         use swiftest_globals
         use swiftest_data_structures
         use helio
         use symba
         implicit none
         logical(lgt), intent(inout) :: ldiscards
         integer(I4B), intent(in)   :: npl
         real(DP), intent(in)   :: t, msys, qmin, qmin_alo, qmin_ahi
         character(*), intent(in)   :: qmin_coord
         type(symba_pl), intent(inout)   :: symba_plA
      end subroutine symba_discard_peri_pl

      subroutine symba_discard_pl(t, npl, nplmax, nsp, symba_plA, rmin, rmax, rmaxu, qmin, qmin_coord,   &
         qmin_alo, qmin_ahi, j2rp2, j4rp4, eoffset)
         use swiftest_globals
         use swiftest_data_structures
         use helio
         use symba
         implicit none
         integer(I4B), intent(in)   :: nplmax
         integer(I4B), intent(inout) :: npl, nsp
         real(DP), intent(in)   :: t, rmin, rmax, rmaxu, qmin, qmin_alo, qmin_ahi, j2rp2, j4rp4
         real(DP), intent(inout)   :: eoffset
         character(*), intent(in)   :: qmin_coord
         type(symba_pl), intent(inout)   :: symba_plA
      end subroutine symba_discard_pl

      subroutine symba_discard_sun_pl(t, npl, msys, swiftest_plA, rmin, rmax, rmaxu, ldiscards)
         use swiftest_globals
         use swiftest_data_structures
         implicit none
         logical(lgt), intent(inout) :: ldiscards
         integer(I4B), intent(in)   :: npl
         real(DP), intent(in)   :: t, msys, rmin, rmax, rmaxu
         type(swiftest_pl), intent(inout)   :: swiftest_plA
      end subroutine symba_discard_sun_pl

      subroutine symba_discard_tp(t, npl, ntp, nsp, symba_plA, symba_tpA, dt, &
         rmin, rmax, rmaxu, qmin, qmin_coord, qmin_alo, qmin_ahi, lclose, lrhill_present)
         use swiftest_globals
         use swiftest_data_structures
         use helio
         use symba
         implicit none
         logical(lgt), intent(in)   :: lclose, lrhill_present
         integer(I4B), intent(in)   :: npl
         integer(I4B), intent(inout) :: ntp, nsp
         real(DP), intent(in)   :: t, dt, rmin, rmax, rmaxu, qmin, qmin_alo, qmin_ahi
         character(*), intent(in)   :: qmin_coord
         type(symba_pl), intent(inout)   :: symba_plA
         type(symba_tp), intent(inout)   :: symba_tpA
      end subroutine symba_discard_tp

      subroutine symba_energy(npl, nplmax, swiftest_plA, j2rp2, j4rp4, ke, pe, te, htot)
         use swiftest_globals
         use swiftest_data_structures
         implicit none
         integer(I4B), intent(in)      :: npl, nplmax
         real(DP), intent(in)       :: j2rp2, j4rp4
         real(DP), intent(out)      :: ke, pe, te
         real(DP), dimension(ndim), intent(out) :: htot
         type(swiftest_pl), intent(inout)      :: swiftest_plA
      end subroutine symba_energy

      subroutine symba_fragmentation(t, dt, index_enc, nmergeadd, nmergesub, mergeadd_list, &
         mergesub_list, eoffset, vbs, encounter_file, out_type, npl, ntp, &
         symba_plA, symba_tpA, nplplenc, npltpenc, pltpenc_list, plplenc_list, &
         nplmax, ntpmax, fragmax)
         use swiftest_globals
         use swiftest_data_structures
         use helio
         use symba
         implicit none
         integer(I4B), intent(in)       :: index_enc, nplmax, ntpmax
         integer(I4B), intent(inout)        :: nmergeadd, nmergesub, nplplenc, npltpenc, fragmax
         integer(I4B), intent(inout)        :: npl, ntp
         real(DP), intent(in)          :: t, dt
         real(DP), intent(inout)        :: eoffset
         real(DP), dimension(ndim), intent(in)     :: vbs
         character(*), intent(in)       :: encounter_file, out_type
         type(symba_plplenc), intent(inout)      :: plplenc_list
         type(symba_pltpenc), intent(inout)      :: pltpenc_list
         type(symba_merger), intent(inout)     :: mergeadd_list, mergesub_list
         type(symba_pl), intent(inout)      :: symba_plA
         type(symba_tp), intent(inout)      :: symba_tpA
      end subroutine symba_fragmentation

      subroutine symba_getacch(lextra_force, t, npl, nplm, nplmax, symba_plA, j2rp2, j4rp4, nplplenc, &
         plplenc_list)
         use swiftest_globals
         use swiftest_data_structures
         use helio
         use symba
         implicit none
         logical(lgt), intent(in)        :: lextra_force
         integer(I4B), intent(in)        :: npl, nplm, nplmax, nplplenc
         real(DP), intent(in)        :: t, j2rp2, j4rp4
         type(symba_pl), intent(inout)      :: symba_plA
         type(symba_plplenc), intent(in)      :: plplenc_list
      end subroutine symba_getacch

      subroutine symba_getacch_tp(lextra_force, t, npl, nplm, nplmax, ntp, ntpmax, symba_plA, symba_tpA, &
         xh, j2rp2, j4rp4,  &
         npltpenc, pltpenc_list)
         use swiftest_globals
         use swiftest_data_structures
         use helio
         use symba
         implicit none
         logical(lgt), intent(in)        :: lextra_force
         integer(I4B), intent(in)        :: npl, nplm, nplmax, ntp, ntpmax, npltpenc
         real(DP), intent(in)        :: t, j2rp2, j4rp4
         real(DP), dimension(ndim, npl), intent(in)   :: xh
         type(symba_pl), intent(inout)      :: symba_plA
         type(symba_tp), intent(inout)      :: symba_tpA
         type(symba_pltpenc), intent(in)      :: pltpenc_list
      end subroutine symba_getacch_tp


      subroutine symba_getacch_eucl(lextra_force, t, npl, nplm, nplmax, symba_plA, j2rp2, j4rp4, nplplenc, &
         plplenc_list, num_plpl_comparisons, k_plpl)
         use swiftest_globals
         use symba
         implicit none
         logical(lgt), intent(in)        :: lextra_force
         integer(I4B), intent(in)        :: npl, nplm, nplmax, nplplenc, num_plpl_comparisons
         real(DP), intent(in)        :: t, j2rp2, j4rp4
         type(symba_pl), intent(inout)      :: symba_plA
         type(symba_plplenc), intent(in)      :: plplenc_list
         integer(I4B), dimension(num_plpl_comparisons,2),intent(in) :: k_plpl
      end subroutine symba_getacch_eucl

      subroutine symba_getacch_tp_eucl(lextra_force, t, npl, nplm, nplmax, ntp, ntpmax, symba_plA, symba_tpA, &
         xh, j2rp2, j4rp4, npltpenc, pltpenc_list, num_pltp_comparisons, k_pltp)
         use swiftest_globals
         use symba
         implicit none
         logical(lgt), intent(in)        :: lextra_force
         integer(I4B), intent(in)        :: npl, nplm, nplmax, ntp, ntpmax, npltpenc, num_pltp_comparisons
         real(DP), intent(in)        :: t, j2rp2, j4rp4
         real(DP), dimension(ndim, npl), intent(in)   :: xh
         type(symba_pl), intent(inout)      :: symba_plA
         type(symba_tp), intent(inout)      :: symba_tpA
         type(symba_pltpenc), intent(in)      :: pltpenc_list
         integer(I4B), dimension(num_pltp_comparisons,2), intent(in) :: k_pltp
      end subroutine symba_getacch_tp_eucl

      subroutine symba_helio_drift(irec, npl, symba_plA, dt)
         use swiftest_globals
         use symba
         implicit none
         integer(I4B), intent(in)   :: irec, npl
         real(DP), intent(in)   :: dt
         type(symba_pl), intent(inout) :: symba_plA
      end subroutine symba_helio_drift

      subroutine symba_helio_drift_tp(irec, ntp, symba_tpA, mu, dt)
         use swiftest_globals
         use symba
         implicit none
         integer(I4B), intent(in)   :: irec, ntp
         real(DP), intent(in)   :: mu, dt
         type(symba_tp), intent(inout) :: symba_tpA
      end subroutine symba_helio_drift_tp

      subroutine symba_helio_getacch(lflag, lextra_force, t, npl, nplm, nplmax, helio_plA, j2rp2, j4rp4)
         use swiftest_globals
         use helio
         implicit none
         logical(lgt), intent(in)   :: lflag, lextra_force
         integer(I4B), intent(in)   :: npl, nplm, nplmax
         real(DP), intent(in)    :: t, j2rp2, j4rp4
         type(helio_pl), intent(inout)  :: helio_plA
      end subroutine symba_helio_getacch

      subroutine symba_helio_getacch_int(npl, nplm, helio_plA)
         use swiftest_globals
         use helio
         implicit none
         integer(I4B), intent(in)   :: npl, nplm
         type(helio_pl), intent(inout) :: helio_plA
      end subroutine symba_helio_getacch_int

      subroutine symba_kick(irec, nplplenc, npltpenc, plplenc_list, pltpenc_list, dt, sgn, symba_plA, &
         symba_tpA)
         use swiftest_globals
         use symba
         implicit none
         integer(I4B), intent(in)   :: irec, nplplenc, npltpenc
         real(DP), intent(in)     :: dt, sgn
         type(symba_plplenc), intent(in) :: plplenc_list
         type(symba_pltpenc), intent(in) :: pltpenc_list
         type(symba_pl), intent(inout)   :: symba_plA
         type(symba_tp), intent(inout)   :: symba_tpA
      end subroutine symba_kick

      subroutine symba_merge_pl(t, dt, index_enc, nplplenc, plplenc_list, nmergeadd, nmergesub, &
         mergeadd_list, mergesub_list, eoffset, vbs, encounter_file, out_type, npl, symba_plA, &
         symba_tpA)
         use swiftest_globals
         use symba
         implicit none
         integer(I4B), intent(in)       :: index_enc, nplplenc
         integer(I4B), intent(inout)        :: nmergeadd, nmergesub, npl
         real(DP), intent(in)          :: t, dt
         real(DP), intent(inout)        :: eoffset
         real(DP), dimension(ndim), intent(in)     :: vbs
         character(*), intent(in)       :: encounter_file, out_type
         type(symba_plplenc), intent(inout)      :: plplenc_list
         type(symba_merger),  intent(inout)      :: mergeadd_list, mergesub_list
         type(symba_pl), intent(inout)      :: symba_plA
         type(symba_tp), intent(inout)      :: symba_tpA
      end subroutine symba_merge_pl

      subroutine symba_merge_tp(t, dt, index_enc, npltpenc, pltpenc_list, vbs, encounter_file, out_type, symba_plA, symba_tpA)
         use swiftest_globals
         use module_swifter
         use helio
         use symba
         implicit none
         integer(I4B), intent(in)      :: index_enc, npltpenc
         real(DP), intent(in)       :: t, dt
         real(DP), dimension(ndim), intent(in)  :: vbs
         character(*), intent(in)      :: encounter_file, out_type
         type(symba_pltpenc), intent(inout)   :: pltpenc_list
         type(symba_pl), intent(inout)      :: symba_plA
         type(symba_tp), intent(inout)      :: symba_tpA
      end subroutine symba_merge_tp

      subroutine symba_peri(lfirst, npl, symba_plA, msys, qmin_coord)
         use swiftest_globals
         use symba
         implicit none
         logical(lgt), intent(in)   :: lfirst
         integer(I4B), intent(in)   :: npl
         real(DP), intent(in)    :: msys
         character(*), intent(in)   :: qmin_coord
         type(symba_pl), intent(inout)  :: symba_plA
      end subroutine symba_peri

      subroutine symba_rearray(t, npl, ntp, nsppl, nsptp, symba_plA, symba_tpA, nmergeadd, &
         mergeadd_list, discard_plA, discard_tpA, nplmax, j2rp2, j4rp4)
         use swiftest_globals
         use swiftest_data_structures
         use helio
         use symba
         implicit none
         integer(I4B), intent(inout)        :: npl, ntp, nsppl, nsptp, nmergeadd, nplmax !change to fragadd
         real(DP), intent(in)          :: t, j2rp2, j4rp4
         type(symba_pl), intent(inout)      :: symba_plA
         type(symba_tp), intent(inout)      :: symba_tpA
         type(swiftest_tp), intent(inout)      :: discard_tpA
         type(swiftest_pl), intent(inout)      :: discard_plA
         type(symba_merger), intent(inout)     :: mergeadd_list !change to fragadd_list

      end subroutine symba_rearray

      subroutine symba_reorder_pl(npl, symba_plA)
         use swiftest_globals
         use swiftest_data_structures
         use helio
         use symba
         implicit none
         integer(I4B), intent(in) :: npl
         type(symba_pl), intent(inout)  :: symba_plA
         integer(I4B)         :: i
         integer(I4B), dimension(:), allocatable   :: index
         real(DP), dimension(:), allocatable   :: mass
         real(DP), dimension(:,:), allocatable   :: symba_plwkspa
         integer(I4B), dimension(:,:), allocatable :: symba_plwkspa_id_status
      end subroutine symba_reorder_pl

      subroutine symba_setup(npl, ntp, symba_plA, symba_tpA, symba_pl1p, symba_tp1p, swiftest_pl1p, &
         swiftest_tp1p)
         use swiftest_globals
         use swiftest_data_structures
         use helio
         use symba
         implicit none
         integer(I4B), intent(in)         :: npl, ntp
         type(swiftest_pl), pointer         :: swiftest_pl1p
         type(swiftest_tp), pointer         :: swiftest_tp1p
         type(symba_pl), intent(inout) :: symba_plA
         type(symba_tp), intent(inout) :: symba_tpA
         type(symba_pl), pointer          :: symba_pl1p
         type(symba_tp), pointer          :: symba_tp1p
      end subroutine symba_setup

      subroutine symba_step_eucl(lfirst, lextra_force, lclose, t, npl, nplmax, ntp, ntpmax, symba_plA, symba_tpA, j2rp2, j4rp4,&
         dt,nplplenc, npltpenc, plplenc_list, pltpenc_list, nmergeadd, nmergesub, mergeadd_list, mergesub_list, eoffset,&
         mtiny,encounter_file, out_type, num_plpl_comparisons, k_plpl, num_pltp_comparisons, k_pltp)
         use swiftest_globals
         use swiftest_data_structures
         use helio
         use symba
         implicit none
         logical(lgt), intent(in)    :: lextra_force, lclose
         logical(lgt), intent(inout)   :: lfirst
         integer(I4B), intent(in)    :: npl, nplmax, ntp, ntpmax
         integer(I4B), intent(inout)   :: nplplenc, npltpenc, nmergeadd, nmergesub
         real(DP), intent(in)      :: t, j2rp2, j4rp4, dt, mtiny
         real(DP), intent(inout)     :: eoffset
         character(*), intent(in)    :: encounter_file, out_type
         type(symba_pl), intent(inout)   :: symba_plA
         type(symba_tp), intent(inout)   :: symba_tpA
         type(symba_plplenc), intent(inout) :: plplenc_list
         type(symba_pltpenc), intent(inout) :: pltpenc_list
         type(symba_merger), intent(inout)  :: mergeadd_list, mergesub_list
         integer(I4B), intent(in)    :: num_plpl_comparisons, num_pltp_comparisons
         integer(I4B), dimension(2,num_plpl_comparisons),intent(in) :: k_plpl
         integer(I4B), dimension(2,num_pltp_comparisons),intent(in) :: k_pltp
      end subroutine symba_step_eucl

      subroutine symba_step(lfirst, lextra_force, lclose, t, npl, nplmax, ntp, ntpmax, symba_plA, &
         symba_tpA, j2rp2, j4rp4, dt, nplplenc, npltpenc, plplenc_list, pltpenc_list, nmergeadd, &
         nmergesub, mergeadd_list, mergesub_list, eoffset, mtiny, encounter_file, out_type, &
         fragmax)
         use swiftest_globals
         use swiftest_data_structures
         use helio
         use symba
         implicit none
         logical(lgt), intent(in)    :: lextra_force, lclose
         logical(lgt), intent(inout)   :: lfirst
         integer(I4B), intent(in)    :: npl, nplmax, ntp, ntpmax
         integer(I4B), intent(inout)   :: nplplenc, npltpenc, nmergeadd, nmergesub, fragmax
         real(DP), intent(in)      :: t, j2rp2, j4rp4, dt, mtiny
         real(DP), intent(inout)     :: eoffset
         character(*), intent(in)    :: encounter_file, out_type
         type(symba_pl), intent(inout)   :: symba_plA
         type(symba_tp), intent(inout)   :: symba_tpA
         type(symba_plplenc), intent(inout) :: plplenc_list
         type(symba_pltpenc), intent(inout) :: pltpenc_list
         type(symba_merger), intent(inout)  :: mergeadd_list, mergesub_list
      end subroutine symba_step

      ! for testing purposes only _ use with symba_step_test
      subroutine symba_step_helio(lfirst, lextra_force, t, npl, nplm, nplmax, ntp, ntpmax, helio_plA, helio_tpA, j2rp2,   &
         j4rp4, dt)
         use swiftest_globals
         use swiftest_data_structures
         use helio
         use symba
         implicit none
         logical(lgt), intent(in)   :: lextra_force
         logical(lgt), intent(inout)   :: lfirst
         integer(I4B), intent(in)   :: npl, nplm, nplmax, ntp, ntpmax
         real(DP), intent(in)   :: t, j2rp2, j4rp4, dt
         type(helio_pl), intent(inout) :: helio_plA
         type(helio_tp), intent(inout) :: helio_tpA
      end subroutine symba_step_helio

      subroutine symba_step_helio_pl(lfirst, lextra_force, t, npl, nplm, nplmax, helio_plA, j2rp2, j4rp4, dt, xbeg, xend,   &
         ptb, pte)
         use swiftest_globals
         use swiftest_data_structures
         use helio
         use symba
         implicit none
         logical(lgt), intent(in)       :: lextra_force
         logical(lgt), intent(inout)      :: lfirst
         integer(I4B), intent(in)       :: npl, nplm, nplmax
         real(DP), intent(in)       :: t, j2rp2, j4rp4, dt
         real(DP), dimension(ndim, nplm), intent(out) :: xbeg, xend
         real(DP), dimension(ndim), intent(out)   :: ptb, pte
         type(helio_pl), intent(inout)     :: helio_plA
      end subroutine symba_step_helio_pl

      subroutine symba_step_interp_eucl(lextra_force, lclose, t, npl, nplm, nplmax, ntp, ntpmax, symba_plA, symba_tpA, j2rp2,&
         j4rp4, dt, eoffset, mtiny, nplplenc, npltpenc, plplenc_list, pltpenc_list, nmergeadd, nmergesub, mergeadd_list,&
         mergesub_list, encounter_file, out_type, num_plpl_comparisons, k_plpl, num_pltp_comparisons, k_pltp)
         use swiftest_globals
         use swiftest_data_structures
         use helio
         use symba
         implicit none
         logical(lgt), intent(in)    :: lextra_force, lclose
         integer(I4B), intent(in)    :: npl, nplm, nplmax, ntp, ntpmax, nplplenc, npltpenc, num_pltp_comparisons
         integer(I4B), intent(inout)   :: nmergeadd, nmergesub
         real(DP), intent(in)      :: t, j2rp2, j4rp4, dt, mtiny
         real(DP), intent(inout)     :: eoffset
         character(*), intent(in)    :: encounter_file, out_type
         type(symba_pl), intent(inout)   :: symba_plA
         type(symba_tp), intent(inout)   :: symba_tpA
         type(symba_plplenc), intent(inout) :: plplenc_list
         type(symba_pltpenc), intent(inout) :: pltpenc_list
         type(symba_merger), intent(inout)  :: mergeadd_list, mergesub_list
         integer(I4B), intent(in)       :: num_plpl_comparisons
         integer(I4B), dimension(num_plpl_comparisons,2),intent(in) :: k_plpl
         integer(I4B), dimension(2,num_pltp_comparisons),intent(in) :: k_pltp
      end subroutine symba_step_interp_eucl

      subroutine symba_step_interp(lextra_force, lclose, t, npl, nplm, nplmax, ntp, ntpmax, symba_plA, symba_tpA, j2rp2,   &
         j4rp4, dt, eoffset, mtiny, nplplenc, npltpenc, plplenc_list, pltpenc_list, nmergeadd, nmergesub, mergeadd_list,   &
         mergesub_list, encounter_file, out_type, fragmax)
         use swiftest_globals
         use swiftest_data_structures
         use helio
         use symba
         implicit none
         logical(lgt), intent(in)    :: lextra_force, lclose
         integer(I4B), intent(in)    :: npl, nplm, nplmax, ntp, ntpmax, nplplenc, npltpenc
         integer(I4B), intent(inout)   :: nmergeadd, nmergesub, fragmax
         real(DP), intent(in)      :: t, j2rp2, j4rp4, dt, mtiny
         real(DP), intent(inout)     :: eoffset
         character(*), intent(in)    :: encounter_file, out_type
         type(symba_pl), intent(inout)   :: symba_plA
         type(symba_tp), intent(inout)   :: symba_tpA
         type(symba_plplenc), intent(inout) :: plplenc_list
         type(symba_pltpenc), intent(inout) :: pltpenc_list
         type(symba_merger), intent(inout)  :: mergeadd_list, mergesub_list
      end subroutine symba_step_interp


      recursive subroutine symba_step_recur(lclose, t, ireci, npl, nplm, ntp, symba_plA, symba_tpA, dt0, eoffset, nplplenc, &
         npltpenc, plplenc_list, pltpenc_list, nmergeadd, nmergesub, mergeadd_list, mergesub_list, encounter_file, & 
         out_type, nplmax, ntpmax, fragmax)
         use swiftest_globals
         use swiftest_data_structures
         use helio
         use symba
         implicit none
         logical(lgt), intent(in)    :: lclose
         integer(I4B), intent(in)    :: ireci, npl, nplm, ntp, nplplenc, npltpenc, nplmax, ntpmax, fragmax
         integer(I4B), intent(inout)   :: nmergeadd, nmergesub
         real(DP), intent(in)      :: t, dt0
         real(DP), intent(inout)     :: eoffset
         character(*), intent(in)    :: encounter_file, out_type
         type(symba_pl), intent(inout)   :: symba_plA
         type(symba_tp), intent(inout)   :: symba_tpA
         type(symba_plplenc), intent(inout) :: plplenc_list
         type(symba_pltpenc), intent(inout) :: pltpenc_list
         type(symba_merger), intent(inout)  :: mergeadd_list, mergesub_list
      end subroutine symba_step_recur

      subroutine symba_user_getacch(t, npl, symba_plA)
         use swiftest_globals
         use symba
         implicit none
         integer(I4B), intent(in)   :: npl
         real(DP), intent(in)    :: t
         type(symba_pl), intent(inout):: symba_plA
      end subroutine symba_user_getacch

      subroutine symba_user_getacch_tp(t, ntp, symba_tpA)
         use swiftest_globals
         use symba
         implicit none
         integer(I4B), intent(in)   :: ntp
         real(DP), intent(in)    :: t
         type(symba_tp), intent(inout)  :: symba_tpA
      end subroutine symba_user_getacch_tp

      subroutine util_exit(code)
         use swiftest_globals
         implicit none
         integer(I4B), intent(in) :: code
      end subroutine util_exit

      subroutine util_dist_index_plpl(npl, nplm, num_comparisons, k_plpl)
         use swiftest_globals
         use swiftest_data_structures
         implicit none
         integer(I4B), intent(in)  :: npl, nplm
         integer(I4B), dimension(:,:),allocatable,intent(out) :: k_plpl
         integer(I4B), intent(out) :: num_comparisons
      end subroutine

      subroutine util_dist_index_pltp(nplm, ntp, num_comparisons, k_pltp)
         use swiftest_globals
         use swiftest_data_structures
         implicit none
         integer(I4B), intent(in)  :: nplm, ntp
         integer(I4B), dimension(:,:),allocatable,intent(out) :: k_pltp
         integer(I4B), intent(out) :: num_comparisons
      end subroutine util_dist_index_pltp

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

      subroutine util_dist_eucl_pltp(npl, ntp, plAnets, test_particles, num_pltp_comparisons, k_pltp, outvar)
         use swiftest_globals
         use swiftest_data_structures
         use symba
         implicit none
         integer(I4B), intent(in) :: npl, ntp
         integer(I4B), dimension(num_pltp_comparisons,2),intent(in) :: k_pltp
         integer(I4B), intent(in) :: num_pltp_comparisons
         real(DP),dimension(ndim,npl),intent(in) :: plAnets
         real(DP),dimension(ndim,ntp),intent(in) :: test_particles
         real(DP), dimension(ndim,num_pltp_comparisons),intent(inout) :: outvar
      end subroutine

      subroutine util_hills(npl, swiftest_plA)
         use swiftest_globals
         use swiftest_data_structures
         implicit none
         integer(I4B), intent(in)    :: npl
         type(swiftest_pl), intent(inout) :: swiftest_plA
      end subroutine util_hills

      subroutine util_index(arr, index)
         use swiftest_globals
         use module_nrutil
         implicit none
         integer(I4B), dimension(:), intent(out) :: index
         real(DP), dimension(:), intent(in)   :: arr
      end subroutine util_index


      subroutine util_resize_pl(symba_plA, npl_new, npl_old)
         use swiftest_globals
         use symba
         use swiftest_data_structures
         use helio
         use module_nrutil
         implicit none
         type(symba_pl), intent(inout) :: symba_plA
         integer(I4B), intent(in)   :: npl_old, npl_new
      end subroutine util_resize_pl

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

      subroutine util_toupper(string)
         use swiftest_globals
         implicit none
         character(*), intent(inout) :: string
      end subroutine util_toupper

      subroutine util_valid(npl, ntp, swiftest_plA, swiftest_tpA)
         use swiftest_globals
         use swiftest_data_structures
         implicit none
         integer(I4B), intent(in)    :: npl, ntp
         type(swiftest_pl), intent(in) :: swiftest_plA
         type(swiftest_tp), intent(in) :: swiftest_tpA
      end subroutine util_valid

      subroutine util_version
         use swiftest_globals
         implicit none
      end subroutine util_version

      function util_kahan_sum(xsum_current, xi, xerror) 
         use swiftest_globals
         implicit none
         real(DP)     :: util_kahan_sum
         real(DP), intent(in)   :: xsum_current, xi
         real(DP), intent(inout) :: xerror
      end function
   end interface
end module module_interfaces