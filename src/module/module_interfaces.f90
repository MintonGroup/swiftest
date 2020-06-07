!**********************************************************************************************************************************
!
!  Unit Name   : module_interfaces
!  Unit Type   : module
!  Project     : Swiftest
!  Package     : module
!  Language    : Fortran 90/95
!
!  Description : Definition of interfaces of subroutines and functions used in swiftest package
!
!  Input
!    Arguments : N/A
!    Terminal  : N/A
!    File      : N/A
!
!  Output
!    Arguments : N/A
!    Terminal  : N/A
!    File      : N/A
!
!  Invocation  : N/A
!
!  Notes       : 
!
!**********************************************************************************************************************************
MODULE module_interfaces

     IMPLICIT NONE

     INTERFACE
          SUBROUTINE coord_h2b(npl, swiftest_plA, msys)
               USE swiftest
               USE module_swiftest
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN)         :: npl
               REAL(DP), INTENT(OUT)            :: msys
               TYPE(swiftest_pl), INTENT(INOUT) :: swiftest_plA
          END SUBROUTINE coord_h2b
     END INTERFACE

     INTERFACE
          SUBROUTINE coord_h2b_tp(ntp, swiftest_tpA, swiftest_plA)
               USE swiftest
               USE module_swiftest
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN)         :: ntp
               TYPE(swiftest_tp), INTENT(INOUT) :: swiftest_tpA
               TYPE(swiftest_pl), INTENT(INOUT) :: swiftest_plA
          END SUBROUTINE coord_h2b_tp
     END INTERFACE

     INTERFACE
          SUBROUTINE coord_vb2vh(npl, swiftest_plA)
               USE swiftest
               USE module_swiftest
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN)         :: npl
               TYPE(swiftest_pl), INTENT(INOUT) :: swiftest_plA
          END SUBROUTINE coord_vb2vh
     END INTERFACE

     INTERFACE
          SUBROUTINE coord_vb2vh_tp(ntp, swiftest_tpA, vs)
               USE swiftest
               USE module_swiftest
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN)              :: ntp
               REAL(DP), DIMENSION(NDIM), INTENT(IN) :: vs
               TYPE(swiftest_tp), INTENT(INOUT)      :: swiftest_tpA
          END SUBROUTINE coord_vb2vh_tp
     END INTERFACE

     INTERFACE
          SUBROUTINE coord_vh2vb(npl, swiftest_plA, msys)
               USE swiftest
               USE module_swiftest
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN)        :: npl
               REAL(DP), INTENT(OUT)           :: msys
               TYPE(swiftest_pl),INTENT(INOUT) :: swiftest_plA
          END SUBROUTINE coord_vh2vb
     END INTERFACE

     INTERFACE
          SUBROUTINE coord_vh2vb_tp(ntp, swiftest_tpA, vs)
               USE swiftest
               USE module_swiftest
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN)              :: ntp
               REAL(DP), DIMENSION(NDIM), INTENT(IN) :: vs
               TYPE(swiftest_tp), INTENT(INOUT)      :: swiftest_tpA
          END SUBROUTINE coord_vh2vb_tp
     END INTERFACE

     INTERFACE
          SUBROUTINE discard(t, dt, npl, ntp, swiftest_plA, swiftest_tpA, rmin, rmax, rmaxu, qmin,  &
               qmin_alo, qmin_ahi, qmin_coord, lclose, lrhill_present)
               USE swiftest
               USE module_swiftest
               IMPLICIT NONE
               LOGICAL(LGT), INTENT(IN)  :: lclose, lrhill_present
               INTEGER(I4B), INTENT(IN)  :: npl, ntp
               REAL(DP), INTENT(IN)      :: t, dt, rmin, rmax, rmaxu, qmin, qmin_alo, qmin_ahi
               CHARACTER(*), INTENT(IN)  :: qmin_coord
               TYPE(swiftest_pl), INTENT(INOUT) :: swiftest_plA
               TYPE(swiftest_tp), INTENT(INOUT) :: swiftest_tpA
          END SUBROUTINE discard
     END INTERFACE

     INTERFACE
          SUBROUTINE discard_peri(t, npl, ntp, swiftest_plA, swiftest_tpA, msys, qmin, qmin_alo, & 
               qmin_ahi, qmin_coord, lrhill_present)
               USE swiftest
               USE module_swiftest
               IMPLICIT NONE
               LOGICAL(LGT), INTENT(IN)  :: lrhill_present
               INTEGER(I4B), INTENT(IN)  :: npl, ntp
               REAL(DP), INTENT(IN)      :: t, msys, qmin, qmin_alo, qmin_ahi
               CHARACTER(*), INTENT(IN)  :: qmin_coord
               TYPE(swiftest_pl), INTENT(INOUT) :: swiftest_plA
               TYPE(swiftest_tp), INTENT(INOUT) :: swiftest_tpA
          END SUBROUTINE discard_peri
     END INTERFACE

     INTERFACE
          SUBROUTINE discard_pl_close(dx, dv, dt, r2crit, iflag, r2min)
               USE swiftest
               IMPLICIT NONE
               INTEGER(I4B), INTENT(OUT)             :: iflag
               REAL(DP), INTENT(IN)                  :: dt, r2crit
               REAL(DP), DIMENSION(NDIM), INTENT(IN) :: dx, dv
               REAL(DP), INTENT(OUT)                 :: r2min
          END SUBROUTINE discard_pl_close
     END INTERFACE

     INTERFACE
          SUBROUTINE discard_pl(t, dt, npl, ntp, swiftest_plA, swiftest_tpA)
               USE swiftest
               USE module_swiftest
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN)  :: npl, ntp
               REAL(DP), INTENT(IN)      :: t, dt
               TYPE(swiftest_pl), INTENT(INOUT) :: swiftest_plA
               TYPE(swiftest_tp), INTENT(INOUT) :: swiftest_tpA
          END SUBROUTINE discard_pl
     END INTERFACE

     INTERFACE
          SUBROUTINE discard_sun(t, ntp, msys, swifter_tpA, rmin, rmax, rmaxu)
               USE swiftest
               USE module_swiftest
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN)  :: ntp
               REAL(DP), INTENT(IN)      :: t, msys, rmin, rmax, rmaxu
               TYPE(swiftest_tp), INTENT(INOUT) :: swifter_tpA
          END SUBROUTINE discard_sun
     END INTERFACE

     INTERFACE
          SUBROUTINE drift_dan(mu, x0, v0, dt0, iflag)
               USE swiftest
               IMPLICIT NONE
               INTEGER(I4B), INTENT(OUT)                :: iflag
               REAL(DP), INTENT(IN)                     :: mu, dt0
               REAL(DP), DIMENSION(NDIM), INTENT(INOUT) :: x0, v0
          END SUBROUTINE drift_dan
     END INTERFACE

     INTERFACE
          SUBROUTINE drift_kepmd(dm, es, ec, x, s, c)
               USE swiftest
               IMPLICIT NONE
               REAL(DP), INTENT(IN)  :: dm, es, ec
               REAL(DP), INTENT(OUT) :: x, s, c
          END SUBROUTINE drift_kepmd
     END INTERFACE

     INTERFACE
          SUBROUTINE drift_kepu(dt,r0,mu,alpha,u,fp,c1,c2,c3,iflag)
               USE swiftest
               IMPLICIT NONE
               INTEGER(I4B), INTENT(OUT) :: iflag
               REAL(DP), INTENT(IN)      :: dt, r0, mu, alpha, u
               REAL(DP), INTENT(OUT)     :: fp, c1, c2, c3
          END SUBROUTINE drift_kepu
     END INTERFACE

     INTERFACE
          SUBROUTINE drift_kepu_fchk(dt, r0, mu, alpha, u, s, f)
               USE swiftest
               IMPLICIT NONE
               REAL(DP), INTENT(IN)  :: dt, r0, mu, alpha, u, s
               REAL(DP), INTENT(OUT) :: f
          END SUBROUTINE drift_kepu_fchk
     END INTERFACE

     INTERFACE
          SUBROUTINE drift_kepu_guess(dt, r0, mu, alpha, u, s)
               USE swiftest
               IMPLICIT NONE
               REAL(DP), INTENT(IN)  :: dt, r0, mu, alpha, u
               REAL(DP), INTENT(OUT) :: s
          END SUBROUTINE drift_kepu_guess
     END INTERFACE

     INTERFACE
          SUBROUTINE drift_kepu_lag(s, dt, r0, mu, alpha, u, fp, c1, c2, c3, iflag)
               USE swiftest
               IMPLICIT NONE
               INTEGER(I4B), INTENT(OUT) :: iflag
               REAL(DP), INTENT(IN)      :: dt, r0, mu, alpha, u
               REAL(DP), INTENT(INOUT)   :: s
               REAL(DP), INTENT(OUT)     :: fp, c1, c2, c3
          END SUBROUTINE drift_kepu_lag
     END INTERFACE

     INTERFACE
          SUBROUTINE drift_kepu_new(s, dt, r0, mu, alpha, u, fp, c1, c2, c3, iflag)
               USE swiftest
               IMPLICIT NONE
               INTEGER(I4B), INTENT(OUT) :: iflag
               REAL(DP), INTENT(IN)      :: dt, r0, mu, alpha, u
               REAL(DP), INTENT(INOUT)   :: s
               REAL(DP), INTENT(OUT)     :: fp, c1, c2, c3
          END SUBROUTINE drift_kepu_new
     END INTERFACE

     INTERFACE
          SUBROUTINE drift_kepu_p3solve(dt, r0, mu, alpha, u, s, iflag)
               USE swiftest
               IMPLICIT NONE
               INTEGER(I4B), INTENT(OUT) :: iflag
               REAL(DP), INTENT(IN)      :: dt, r0, mu, alpha, u
               REAL(DP), INTENT(OUT)     :: s
          END SUBROUTINE drift_kepu_p3solve
     END INTERFACE

     INTERFACE
          SUBROUTINE drift_kepu_stumpff(x, c0, c1, c2, c3)
               USE swiftest
               IMPLICIT NONE
               REAL(DP), INTENT(INOUT) :: x
               REAL(DP), INTENT(OUT)   :: c0, c1, c2, c3
          END SUBROUTINE drift_kepu_stumpff
     END INTERFACE

     INTERFACE
          SUBROUTINE drift_one(mu, x, v, dt, iflag)
               USE swiftest
               IMPLICIT NONE
               INTEGER(I4B), INTENT(OUT)                :: iflag
               REAL(DP), INTENT(IN)                     :: mu, dt
               REAL(DP), DIMENSION(NDIM), INTENT(INOUT) :: x, v
          END SUBROUTINE drift_one
     END INTERFACE

     INTERFACE
          SUBROUTINE helio_drift(npl, swiftest_plA, dt)
               USE swiftest
               USE module_swiftest
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN)         :: npl
               REAL(DP), INTENT(IN)             :: dt
               TYPE(swiftest_pl), INTENT(INOUT) :: swiftest_plA
          END SUBROUTINE helio_drift
     END INTERFACE

     INTERFACE
          SUBROUTINE helio_drift_tp(ntp, swiftest_tpA, mu, dt)
               USE swiftest
               USE module_swiftest
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN)         :: ntp
               REAL(DP), INTENT(IN)             :: mu, dt
               TYPE(swiftest_tp), INTENT(INOUT) :: swiftest_tpA
          END SUBROUTINE helio_drift_tp
     END INTERFACE

     INTERFACE
          SUBROUTINE helio_getacch(lflag, lextra_force, t, npl, nplmax, helio_plA, j2rp2, j4rp4)
               USE swiftest
               USE module_swiftest
               USE module_helio
               IMPLICIT NONE
               LOGICAL(LGT), INTENT(IN)       :: lflag, lextra_force
               INTEGER(I4B), INTENT(IN)       :: npl, nplmax
               REAL(DP), INTENT(IN)           :: t, j2rp2, j4rp4
               TYPE(helio_pl), INTENT(INOUT)  :: helio_plA
          END SUBROUTINE helio_getacch
     END INTERFACE

     INTERFACE
          SUBROUTINE helio_getacch_int(npl, helio_plA)
               USE swiftest
               USE module_helio
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN)       :: npl
               TYPE(helio_pl), INTENT(INOUT)  :: helio_plA
          END SUBROUTINE helio_getacch_int
     END INTERFACE

     INTERFACE
          SUBROUTINE helio_getacch_int_tp(npl, ntp, swiftest_plA, helio_tpA)
               USE swiftest
               USE module_swiftest
               USE module_helio
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN)                   :: npl, ntp
               TYPE(swiftest_pl), INTENT(INOUT)           :: swiftest_plA
               TYPE(helio_tp), INTENT(INOUT)              :: helio_tpA
          END SUBROUTINE helio_getacch_int_tp
     END INTERFACE

     INTERFACE
          SUBROUTINE helio_getacch_tp(lflag, lextra_force, t, npl, nplmax, ntp, ntpmax, helio_plA, helio_tpA, xh, j2rp2, j4rp4)
               USE swiftest
               USE module_swiftest
               USE module_helio
               IMPLICIT NONE
               LOGICAL(LGT), INTENT(IN)                   :: lflag, lextra_force
               INTEGER(I4B), INTENT(IN)                   :: npl, nplmax, ntp, ntpmax
               REAL(DP), INTENT(IN)                       :: t, j2rp2, j4rp4
               REAL(DP), DIMENSION(NDIM, npl), INTENT(IN) :: xh
               TYPE(helio_pl), INTENT(INOUT)              :: helio_plA
               TYPE(helio_tp), INTENT(INOUT)              :: helio_tpA
          END SUBROUTINE helio_getacch_tp
     END INTERFACE

     INTERFACE
          SUBROUTINE helio_kickvb(npl, helio_plA, dt)
               USE swiftest
               USE module_swiftest
               USE module_helio
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN)       :: npl
               REAL(DP), INTENT(IN)           :: dt
               TYPE(helio_pl), INTENT(INOUT)  :: helio_plA
          END SUBROUTINE helio_kickvb
     END INTERFACE

     INTERFACE
          SUBROUTINE helio_kickvb_tp(ntp, helio_tpA, dt)
               USE swiftest
               USE module_swiftest
               USE module_helio
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN)       :: ntp
               REAL(DP), INTENT(IN)           :: dt
               TYPE(helio_tp), INTENT(INOUT)  :: helio_tpA
          END SUBROUTINE helio_kickvb_tp
     END INTERFACE

     INTERFACE
          SUBROUTINE helio_lindrift(npl, swiftest_plA, dt, pt)
               USE swiftest
               USE module_swiftest
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN)               :: npl
               REAL(DP), INTENT(IN)                   :: dt
               REAL(DP), DIMENSION(NDIM), INTENT(OUT) :: pt
               TYPE(swiftest_pl), INTENT(INOUT)       :: swiftest_plA
          END SUBROUTINE helio_lindrift
     END INTERFACE

     INTERFACE
          SUBROUTINE helio_lindrift_tp(ntp, swiftest_tpA, dt, pt)
               USE swiftest
               USE module_swiftest
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN)              :: ntp
               REAL(DP), INTENT(IN)                  :: dt
               REAL(DP), DIMENSION(NDIM), INTENT(IN) :: pt
               TYPE(swiftest_tp), INTENT(INOUT)      :: swiftest_tpA
          END SUBROUTINE helio_lindrift_tp
     END INTERFACE

     INTERFACE
          SUBROUTINE helio_step(lfirst, lextra_force, t, npl, nplmax, ntp, ntpmax, helio_plA, helio_tpA, j2rp2, j4rp4, dt)
               USE swiftest
               USE module_helio
               IMPLICIT NONE
               LOGICAL(LGT), INTENT(IN)      :: lextra_force
               LOGICAL(LGT), INTENT(INOUT)   :: lfirst
               INTEGER(I4B), INTENT(IN)      :: npl, nplmax, ntp, ntpmax
               REAL(DP), INTENT(IN)          :: t, j2rp2, j4rp4, dt
               TYPE(helio_pl), INTENT(INOUT) :: helio_plA
               TYPE(helio_tp), INTENT(INOUT) :: helio_tpA
          END SUBROUTINE helio_step
     END INTERFACE

     INTERFACE
          SUBROUTINE helio_step_pl(lfirst, lextra_force, t, npl, nplmax, helio_plA, j2rp2, j4rp4, dt, xbeg, xend, ptb, pte)
               USE swiftest
               USE module_swiftest
               USE module_helio
               IMPLICIT NONE
               LOGICAL(LGT), INTENT(IN)                    :: lextra_force
               LOGICAL(LGT), INTENT(INOUT)                 :: lfirst
               INTEGER(I4B), INTENT(IN)                    :: npl, nplmax
               REAL(DP), INTENT(IN)                        :: t, j2rp2, j4rp4, dt
               REAL(DP), DIMENSION(NDIM), INTENT(OUT)      :: ptb, pte
               REAL(DP), DIMENSION(NDIM, npl), INTENT(OUT) :: xbeg,xend
               TYPE(helio_pl), INTENT(INOUT)               :: helio_plA
          END SUBROUTINE helio_step_pl
     END INTERFACE

     INTERFACE
          SUBROUTINE helio_step_tp(lfirsttp, lextra_force, t, npl, nplmax, ntp, ntpmax, helio_plA, helio_tpA, j2rp2, j4rp4, dt, &
               xbeg, xend, ptb, pte)
               USE swiftest
               USE module_swiftest
               USE module_helio
               IMPLICIT NONE
               LOGICAL(LGT), INTENT(IN)                   :: lextra_force
               LOGICAL(LGT), INTENT(INOUT)                :: lfirsttp
               INTEGER(I4B), INTENT(IN)                   :: npl, nplmax, ntp, ntpmax
               REAL(DP), INTENT(IN)                       :: t, j2rp2, j4rp4, dt
               REAL(DP), DIMENSION(NDIM), INTENT(IN)      :: ptb, pte
               REAL(DP), DIMENSION(NDIM, npl), INTENT(IN) :: xbeg, xend
               TYPE(helio_pl), INTENT(INOUT)              :: helio_plA
               TYPE(helio_tp), INTENT(INOUT)              :: helio_tpA
          END SUBROUTINE helio_step_tp
     END INTERFACE

     INTERFACE
          SUBROUTINE helio_user_getacch(t, npl, helio_plA)
               USE swiftest
               USE module_helio
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN)      :: npl
               REAL(DP), INTENT(IN)          :: t
               TYPE(helio_pl), INTENT(INOUT) :: helio_plA
          END SUBROUTINE helio_user_getacch
     END INTERFACE

     INTERFACE
          SUBROUTINE helio_user_getacch_tp(t, ntp, helio_tpA)
               USE swiftest
               USE module_helio
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN)      :: ntp
               REAL(DP), INTENT(IN)          :: t
               TYPE(helio_tp), INTENT(INOUT) :: helio_tpA
          END SUBROUTINE helio_user_getacch_tp
     END INTERFACE

     INTERFACE
          SUBROUTINE io_discard_write_symba(t, mtiny, npl, ntp, nsppl, nsptp, nmergeadd, symba_plA, & 
               discard_plA, discard_tpA, mergeadd_list, mergesub_list, fname, lbig_discard)
               USE swiftest
               USE module_swiftest
               USE module_symba
               IMPLICIT NONE
               LOGICAL(LGT), INTENT(IN)                       :: lbig_discard
               INTEGER(I4B), INTENT(IN)                       :: npl, ntp, nsppl, nsptp, nmergeadd
               REAL(DP), INTENT(IN)                           :: t, mtiny
               CHARACTER(*), INTENT(IN)                       :: fname
               TYPE(symba_pl), INTENT(INOUT)                  :: symba_plA
               TYPE(swiftest_tp), INTENT(INOUT)               :: discard_tpA
               TYPE(swiftest_pl), INTENT(INOUT)               :: discard_plA
               TYPE(symba_merger), INTENT(INOUT)              :: mergeadd_list, mergesub_list
          END SUBROUTINE io_discard_write_symba
     END INTERFACE

     INTERFACE
          SUBROUTINE io_dump_param(nplmax, ntpmax, ntp, t, tstop, dt, in_type, istep_out, outfile, out_type, out_form,            &
               istep_dump, j2rp2, j4rp4,rmin, rmax, rmaxu, qmin, qmin_coord, qmin_alo, qmin_ahi, encounter_file,         &
               mtiny, feature, ring_outfile)
               USE swiftest
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN) :: nplmax, ntpmax, ntp, istep_out, istep_dump
               REAL(DP), INTENT(IN)     :: t, tstop, dt, j2rp2, j4rp4, rmin, rmax, rmaxu, qmin, qmin_alo, qmin_ahi
               CHARACTER(*), INTENT(IN) :: qmin_coord, encounter_file, in_type, outfile, out_type, out_form
               REAl(DP), INTENT(IN), OPTIONAL :: mtiny 
               TYPE(feature_list), INTENT(IN) :: feature
               CHARACTER(*), INTENT(IN), OPTIONAL :: ring_outfile
          END SUBROUTINE io_dump_param
     END INTERFACE

     INTERFACE
          SUBROUTINE io_dump_pl(npl, swiftest_plA, lclose, lrhill_present)
               USE swiftest
               USE module_swiftest
               USE module_fxdr
               IMPLICIT NONE
               LOGICAL(LGT), INTENT(IN)        :: lclose, lrhill_present
               INTEGER(I4B), INTENT(IN)        :: npl
               TYPE(swiftest_pl), INTENT(INOUT):: swiftest_plA
          END SUBROUTINE io_dump_pl
     END INTERFACE

     INTERFACE
          SUBROUTINE io_dump_tp(ntp, swiftest_tpA)
               USE swiftest
               USE module_swiftest
               USE module_fxdr
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN)        :: ntp
               TYPE(swiftest_tp), INTENT(INOUT):: swiftest_tpA
          END SUBROUTINE io_dump_tp
     END INTERFACE

     INTERFACE
          SUBROUTINE io_init_pl(inplfile, in_type, lclose, lrhill_present, npl, symba_plA)
               USE swiftest
               USE module_swiftest
               USE module_symba
               USE module_helio
               USE module_fxdr
               IMPLICIT NONE
               LOGICAL(LGT), INTENT(IN)         :: lclose, lrhill_present
               INTEGER(I4B), INTENT(IN)         :: npl
               CHARACTER(*), INTENT(IN)         :: inplfile, in_type
               TYPE(symba_pl), INTENT(INOUT)    :: symba_plA
          END SUBROUTINE io_init_pl
     END INTERFACE

     INTERFACE
          SUBROUTINE io_init_tp(intpfile, in_type, ntp, symba_tpA)
               USE swiftest
               USE module_swiftest
               USE module_symba
               USE module_helio
               USE module_fxdr
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN)         :: ntp
               CHARACTER(*), INTENT(IN)         :: intpfile, in_type
               TYPE(symba_tp), INTENT(INOUT)    :: symba_tpA
          END SUBROUTINE io_init_tp
     END INTERFACE

     INTERFACE
          SUBROUTINE io_open(iu, fname, fopenstat, fmt, ierr)
               USE swiftest
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN)  :: iu
               INTEGER(I4B), INTENT(OUT) :: ierr
               CHARACTER(*), INTENT(IN)  :: fname, fopenstat, fmt
          END SUBROUTINE io_open
     END INTERFACE

     INTERFACE
          SUBROUTINE io_open_fxdr(fname, fopenstat, lflag, iu, ierr)
               USE swiftest
               USE module_fxdr
               IMPLICIT NONE
               LOGICAL(LGT), INTENT(IN)  :: lflag
               INTEGER(I4B), INTENT(OUT) :: iu, ierr
               CHARACTER(*), INTENT(IN)  :: fname
               CHARACTER(1), INTENT(IN)  :: fopenstat
          END SUBROUTINE io_open_fxdr
     END INTERFACE

     INTERFACE
          FUNCTION io_read_encounter(t, name1, name2, mass1, mass2, xh1, xh2, vh1, vh2, encounter_file, out_type)
               USE swiftest
               USE module_fxdr
               IMPLICIT NONE
               INTEGER(I4B)                           :: io_read_encounter
               INTEGER(I4B), INTENT(OUT)              :: name1, name2
               REAL(DP), INTENT(OUT)                  :: t, mass1, mass2
               REAL(DP), DIMENSION(NDIM), INTENT(OUT) :: xh1, xh2, vh1, vh2
               CHARACTER(*), INTENT(IN)               :: encounter_file,out_type
          END FUNCTION io_read_encounter
     END INTERFACE

     INTERFACE
          FUNCTION io_read_hdr(iu, t, npl, ntp, iout_form, out_type)
               USE swiftest
               USE module_fxdr
               IMPLICIT NONE
               INTEGER(I4B)               :: io_read_hdr
               INTEGER(I4B), INTENT(IN)   :: iu
               INTEGER(I4B), INTENT(OUT)  :: npl, ntp, iout_form
               REAL(DP), INTENT(OUT)      :: t
               CHARACTER(*), INTENT(IN)   :: out_type
          END FUNCTION io_read_hdr
     END INTERFACE

     INTERFACE
          FUNCTION io_read_line(iu, name, d1, d2, d3, d4, d5, d6, out_type, MASS, RADIUS)
               USE swiftest
               USE module_fxdr
               IMPLICIT NONE
               INTEGER(I4B)                    :: io_read_line
               INTEGER(I4B), INTENT(IN)        :: iu
               INTEGER(I4B), INTENT(OUT)       :: name
               REAL(DP), INTENT(OUT)           :: d1, d2, d3, d4, d5, d6
               REAL(DP), OPTIONAL, INTENT(OUT) :: MASS, RADIUS
               CHARACTER(*), INTENT(IN)        :: out_type
          END FUNCTION io_read_line
     END INTERFACE

     INTERFACE
          SUBROUTINE io_write_encounter(t, name1, name2, mass1, mass2, radius1, radius2, &
               xh1, xh2, vh1, vh2, encounter_file, out_type)
               USE swiftest
               USE module_fxdr
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN)              :: name1, name2
               REAL(DP), INTENT(IN)                  :: t, mass1, mass2, radius1, radius2
               REAL(DP), DIMENSION(NDIM), INTENT(IN) :: xh1, xh2, vh1, vh2
               CHARACTER(*), INTENT(IN)              :: encounter_file, out_type
          END SUBROUTINE io_write_encounter
     END INTERFACE

     INTERFACE
          SUBROUTINE io_write_frame(t, npl, ntp, swiftest_plA, swiftest_tpA, outfile, &
               out_type, out_form, out_stat)
               USE swiftest
               USE module_swiftest
               USE module_fxdr
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN)  :: npl, ntp
               REAL(DP), INTENT(IN)      :: t
               CHARACTER(*), INTENT(IN)  :: outfile, out_type, out_form, out_stat
               TYPE(swiftest_pl), INTENT(INOUT) :: swiftest_plA
               TYPE(swiftest_tp), INTENT(INOUT) :: swiftest_tpA
          END SUBROUTINE io_write_frame
     END INTERFACE

     INTERFACE
          SUBROUTINE io_write_hdr(iu, t, npl, ntp, iout_form, out_type)
               USE swiftest
               USE module_fxdr
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN) :: iu, npl, ntp, iout_form
               REAL(DP), INTENT(IN)     :: t
               CHARACTER(*), INTENT(IN) :: out_type
          END SUBROUTINE io_write_hdr
     END INTERFACE

     INTERFACE
          SUBROUTINE io_write_line(iu, name, d1, d2, d3, d4, d5, d6, out_type, MASS, RADIUS)
               USE swiftest
               USE module_fxdr
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN)       :: iu, name
               REAL(DP), INTENT(IN)           :: d1, d2, d3, d4, d5, d6
               REAL(DP), OPTIONAL, INTENT(IN) :: MASS, RADIUS
               CHARACTER(*), INTENT(IN)       :: out_type
          END SUBROUTINE io_write_line
     END INTERFACE

     INTERFACE
          SUBROUTINE obl_acc(npl, swiftest_plA, j2rp2, j4rp4, xh, irh, aobl)
               USE swiftest
               USE module_swiftest
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN)                    :: npl
               REAL(DP), INTENT(IN)                        :: j2rp2, j4rp4
               REAL(DP), DIMENSION(npl), INTENT(IN)        :: irh
               REAL(DP), DIMENSION(NDIM, npl), INTENT(IN)  :: xh
               REAL(DP), DIMENSION(NDIM, npl), INTENT(OUT) :: aobl
               TYPE(swiftest_pl), INTENT(INOUT)            :: swiftest_plA
          END SUBROUTINE obl_acc
     END INTERFACE

     INTERFACE
          SUBROUTINE obl_acc_tp(ntp, xht, j2rp2, j4rp4, irht, aoblt, msun)
               USE swiftest
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN)                    :: ntp
               REAL(DP), INTENT(IN)                        :: j2rp2, j4rp4, msun
               REAL(DP), DIMENSION(ntp), INTENT(IN)        :: irht
               REAL(DP), DIMENSION(NDIM, ntp), INTENT(IN)  :: xht
               REAL(DP), DIMENSION(NDIM, ntp), INTENT(OUT) :: aoblt
          END SUBROUTINE obl_acc_tp
     END INTERFACE

     INTERFACE
          SUBROUTINE obl_pot(npl, swiftest_plA, j2rp2, j4rp4, xh, irh, oblpot)
               USE swiftest
               USE module_swiftest
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN)                   :: npl
               REAL(DP), INTENT(IN)                       :: j2rp2, j4rp4
               REAL(DP), INTENT(OUT)                      :: oblpot
               REAL(DP), DIMENSION(npl), INTENT(IN)       :: irh
               REAL(DP), DIMENSION(NDIM, npl), INTENT(IN) :: xh
               TYPE(swiftest_pl), INTENT(INOUT)           :: swiftest_plA
          END SUBROUTINE obl_pot
     END INTERFACE

     INTERFACE
          SUBROUTINE orbel_scget(angle, sx, cx)
               USE swiftest
               IMPLICIT NONE
               REAL(DP), INTENT(IN)  :: angle
               REAL(DP), INTENT(OUT) :: sx, cx
          END SUBROUTINE orbel_scget
     END INTERFACE

     INTERFACE
          SUBROUTINE orbel_xv2aeq(x, v, mu, a, e, q)
               USE swiftest
               IMPLICIT NONE
               REAL(DP), INTENT(IN)                  :: mu
               REAL(DP), DIMENSION(NDIM), INTENT(IN) :: x, v
               REAL(DP), INTENT(OUT)                 :: a, e, q
          END SUBROUTINE orbel_xv2aeq
     END INTERFACE

     INTERFACE
          SUBROUTINE orbel_xv2aqt(x, v, mu, a, q, capm, tperi)
               USE swiftest
               IMPLICIT NONE
               REAL(DP), INTENT(IN)                  :: mu
               REAL(DP), DIMENSION(NDIM), INTENT(IN) :: x, v
               REAL(DP), INTENT(OUT)                 :: a, q, capm, tperi
          END SUBROUTINE orbel_xv2aqt
     END INTERFACE

     INTERFACE
          SUBROUTINE orbel_xv2el(x, v, mu, a, e, inc, capom, omega, capm)
               USE swiftest
               IMPLICIT NONE
               REAL(DP), INTENT(IN)                  :: mu
               REAL(DP), DIMENSION(NDIM), INTENT(IN) :: x, v
               REAL(DP), INTENT(OUT)                 :: a, e, inc, capom, omega, capm
          END SUBROUTINE orbel_xv2el
     END INTERFACE

     INTERFACE
          SUBROUTINE python_io_write_frame_pl(t, symba_plA, npl, out_stat)
               use swiftest
               use module_swiftest
               use module_helio
               use module_symba
               IMPLICIT NONE
               real(DP), intent(in)      :: t
               type(symba_pl),intent(in) :: symba_plA
               integer, intent(in)       :: npl
               character(*), intent(in)  :: out_stat
          END SUBROUTINE python_io_write_frame_pl
     END INTERFACE

     INTERFACE
          SUBROUTINE python_io_write_frame_tp(t, symba_tpA, ntp, out_stat)
               use swiftest
               use module_swiftest
               use module_helio
               use module_symba
               IMPLICIT NONE
               real(DP), intent(in)      :: t
               type(symba_tp),intent(in) :: symba_tpA
               integer, intent(in)       :: ntp
               character(*), intent(in)  :: out_stat
          END SUBROUTINE python_io_write_frame_tp
     END INTERFACE



     INTERFACE
          SUBROUTINE rmvs_chk_ind(xr, vr, dt, r2crit, iflag)
               USE swiftest
               IMPLICIT NONE
               REAL(DP), INTENT(IN)                  :: dt, r2crit
               REAL(DP), DIMENSION(NDIM), INTENT(IN) :: xr, vr
               INTEGER(I4B), INTENT(OUT)             :: iflag
          END SUBROUTINE rmvs_chk_ind
     END INTERFACE

     INTERFACE
          SUBROUTINE symba_casedisruption (t, dt, index_enc, nmergeadd, nmergesub, mergeadd_list, mergesub_list, eoffset, vbs, & 
          symba_plA, nplplenc, plplenc_list, nplmax, ntpmax, fragmax, mres, rres, m1, m2, rad1, rad2, x1, x2, v1, v2)
          USE swiftest
          USE module_swiftest
          USE module_helio
          USE module_symba
          IMPLICIT NONE
          INTEGER(I4B), INTENT(IN)                         :: index_enc, nplmax, ntpmax
          INTEGER(I4B), INTENT(INOUT)                      :: nmergeadd, nmergesub, nplplenc, fragmax
          REAL(DP), INTENT(IN)                             :: t, dt
          REAL(DP), INTENT(INOUT)                          :: eoffset, m1, m2, rad1, rad2
          REAL(DP), DIMENSION(3), INTENT(INOUT)            :: mres, rres
          REAL(DP), DIMENSION(NDIM), INTENT(IN)            :: vbs
          REAL(DP), DIMENSION(NDIM), INTENT(INOUT)         :: x1, x2, v1, v2
          TYPE(symba_plplenc), INTENT(INOUT)               :: plplenc_list
          TYPE(symba_merger), INTENT(INOUT)                :: mergeadd_list, mergesub_list
          TYPE(symba_pl), INTENT(INOUT)                    :: symba_plA

          END SUBROUTINE symba_casedisruption
     END INTERFACE

     INTERFACE
          SUBROUTINE symba_casehitandrun (t, index_enc, nmergeadd, nmergesub, mergeadd_list, mergesub_list, eoffset, vbs, & 
          symba_plA, nplplenc, plplenc_list, nplmax, ntpmax, fragmax, mres, rres, m1, m2, rad1, rad2, x1, x2, v1, v2)
          USE swiftest
          USE module_swiftest
          USE module_helio
          USE module_symba
          IMPLICIT NONE
          INTEGER(I4B), INTENT(IN)                         :: index_enc, nplmax, ntpmax
          INTEGER(I4B), INTENT(INOUT)                      :: nmergeadd, nmergesub, nplplenc, fragmax
          REAL(DP), INTENT(IN)                             :: t
          REAL(DP), INTENT(INOUT)                          :: eoffset, m1, m2, rad1, rad2
          REAL(DP), DIMENSION(3), INTENT(INOUT)            :: mres, rres
          REAL(DP), DIMENSION(NDIM), INTENT(IN)            :: vbs
          REAL(DP), DIMENSION(NDIM), INTENT(INOUT)         :: x1, x2, v1, v2
          TYPE(symba_plplenc), INTENT(INOUT)               :: plplenc_list
          TYPE(symba_merger), INTENT(INOUT)                :: mergeadd_list, mergesub_list
          TYPE(symba_pl), INTENT(INOUT)                    :: symba_plA

          END SUBROUTINE symba_casehitandrun
     END INTERFACE

     INTERFACE
          SUBROUTINE symba_casemerge (t, dt, index_enc, nmergeadd, nmergesub, mergeadd_list, mergesub_list, eoffset, vbs, & 
          encounter_file, out_type, npl, ntp, symba_plA, symba_tpA, nplplenc, npltpenc, pltpenc_list, plplenc_list, &
          array_index1_child, array_index2_child, m1, m2, rad1, rad2, x1, x2, v1, v2)
          USE swiftest
          USE module_swiftest
          USE module_helio
          USE module_symba
          IMPLICIT NONE
          INTEGER(I4B), INTENT(IN)                         :: index_enc
          INTEGER(I4B), INTENT(INOUT)                      :: npl, ntp, nmergeadd, nmergesub, nplplenc, npltpenc
          REAL(DP), INTENT(IN)                             :: t, dt
          REAL(DP), INTENT(INOUT)                          :: eoffset, m1, m2, rad1, rad2
          REAL(DP), DIMENSION(NDIM), INTENT(IN)            :: vbs
          REAL(DP), DIMENSION(NDIM), INTENT(INOUT)         :: x1, x2, v1, v2
          CHARACTER(*), INTENT(IN)                         :: encounter_file, out_type
          TYPE(symba_plplenc), INTENT(INOUT)               :: plplenc_list
          TYPE(symba_pltpenc), INTENT(INOUT)               :: pltpenc_list
          TYPE(symba_merger), INTENT(INOUT)                :: mergeadd_list, mergesub_list
          TYPE(symba_pl), INTENT(INOUT)                    :: symba_plA
          TYPE(symba_tp), INTENT(INOUT)                    :: symba_tpA
          INTEGER(I4B), DIMENSION(npl), INTENT(INOUT)      :: array_index1_child, array_index2_child

          END SUBROUTINE symba_casemerge
     END INTERFACE

     INTERFACE
          SUBROUTINE symba_caseresolve (t, dt, index_enc, nmergeadd, nmergesub, mergeadd_list, mergesub_list, eoffset, vbs, & 
          encounter_file, out_type, npl, ntp, symba_plA, symba_tpA, nplplenc, npltpenc, pltpenc_list, plplenc_list, regime, &
          nplmax, ntpmax, fragmax, mres, rres, array_index1_child, array_index2_child, m1, m2, rad1, rad2, x1, x2, v1, v2)
          USE swiftest
          USE module_swiftest
          USE module_helio
          USE module_symba
          IMPLICIT NONE
          INTEGER(I4B), INTENT(IN)                         :: index_enc, nplmax, ntpmax
          INTEGER(I4B), INTENT(INOUT)                      :: npl, ntp, nmergeadd, nmergesub, nplplenc, npltpenc, fragmax
          REAL(DP), INTENT(IN)                             :: t, dt
          REAL(DP), INTENT(INOUT)                          :: eoffset, m1, m2, rad1, rad2
          REAL(DP), DIMENSION(3), INTENT(INOUT)            :: mres, rres
          REAL(DP), DIMENSION(NDIM), INTENT(IN)            :: vbs
          REAL(DP), DIMENSION(NDIM), INTENT(INOUT)         :: x1, x2, v1, v2
          CHARACTER(*), INTENT(IN)                         :: encounter_file, out_type
          TYPE(symba_plplenc), INTENT(INOUT)               :: plplenc_list
          TYPE(symba_pltpenc), INTENT(INOUT)               :: pltpenc_list
          TYPE(symba_merger), INTENT(INOUT)                :: mergeadd_list, mergesub_list
          TYPE(symba_pl), INTENT(INOUT)                    :: symba_plA
          TYPE(symba_tp), INTENT(INOUT)                    :: symba_tpA
          INTEGER(I4B), INTENT(IN)                         :: regime 
          INTEGER(I4B), DIMENSION(npl), INTENT(INOUT)      :: array_index1_child, array_index2_child

          END SUBROUTINE symba_caseresolve
     END INTERFACE

     INTERFACE
          SUBROUTINE symba_casesupercatastrophic (t, dt, index_enc, nmergeadd, nmergesub, mergeadd_list, mergesub_list, &
               eoffset, vbs, symba_plA, nplplenc, plplenc_list, nplmax, ntpmax, fragmax, mres, rres, m1, m2, rad1, &
               rad2, x1, x2, v1, v2)
          USE swiftest
          USE module_swiftest
          USE module_helio
          USE module_symba
          IMPLICIT NONE
          INTEGER(I4B), INTENT(IN)                         :: index_enc, nplmax, ntpmax
          INTEGER(I4B), INTENT(INOUT)                      :: nmergeadd, nmergesub, nplplenc, fragmax
          REAL(DP), INTENT(IN)                             :: t, dt
          REAL(DP), INTENT(INOUT)                          :: eoffset, m1, m2, rad1, rad2
          REAL(DP), DIMENSION(3), INTENT(INOUT)            :: mres, rres
          REAL(DP), DIMENSION(NDIM), INTENT(IN)            :: vbs
          REAL(DP), DIMENSION(NDIM), INTENT(INOUT)         :: x1, x2, v1, v2
          TYPE(symba_plplenc), INTENT(INOUT)               :: plplenc_list
          TYPE(symba_merger), INTENT(INOUT)                :: mergeadd_list, mergesub_list
          TYPE(symba_pl), INTENT(INOUT)                    :: symba_plA

          END SUBROUTINE symba_casesupercatastrophic
     END INTERFACE

     INTERFACE
          SUBROUTINE symba_chk(xr, vr, rhill1, rhill2, dt, irec, lencounter, lvdotr)
               USE swiftest
               USE module_swiftest
               USE module_helio
               USE module_symba
               IMPLICIT NONE
               LOGICAL(LGT), INTENT(OUT)          :: lencounter, lvdotr
               INTEGER(I4B), INTENT(IN)           :: irec
               REAL(DP), INTENT(IN)               :: rhill1, rhill2, dt
               REAL(DP), DIMENSION(:), INTENT(IN) :: xr, vr
          END SUBROUTINE symba_chk
     END INTERFACE

     INTERFACE
          SUBROUTINE symba_discard_merge_pl(t, npl, symba_plA, nplplenc, plplenc_list)
               USE swiftest
               USE module_swiftest
               USE module_helio
               USE module_symba
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN)                      :: nplplenc
               INTEGER(I4B), INTENT(INOUT)                   :: npl
               REAL(DP), INTENT(IN)                          :: t
               TYPE(symba_pl)                                :: symba_plA
               TYPE(symba_plplenc), INTENT(IN)               :: plplenc_list
          END SUBROUTINE symba_discard_merge_pl
     END INTERFACE

     INTERFACE
          SUBROUTINE symba_discard_peri_pl(t, npl, symba_plA, msys, qmin, qmin_alo, qmin_ahi, qmin_coord, ldiscards)
               USE swiftest
               USE module_swiftest
               USE module_helio
               USE module_symba
               IMPLICIT NONE
               LOGICAL(LGT), INTENT(INOUT) :: ldiscards
               INTEGER(I4B), INTENT(IN)    :: npl
               REAL(DP), INTENT(IN)        :: t, msys, qmin, qmin_alo, qmin_ahi
               CHARACTER(*), INTENT(IN)    :: qmin_coord
               TYPE(symba_pl), INTENT(INOUT)     :: symba_plA
          END SUBROUTINE symba_discard_peri_pl
     END INTERFACE

     INTERFACE
          SUBROUTINE symba_discard_pl(t, npl, nplmax, nsp, symba_plA, rmin, rmax, rmaxu, qmin, qmin_coord,          &
               qmin_alo, qmin_ahi, j2rp2, j4rp4, eoffset)
               USE swiftest
               USE module_swiftest
               USE module_helio
               USE module_symba
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN)    :: nplmax
               INTEGER(I4B), INTENT(INOUT) :: npl, nsp
               REAL(DP), INTENT(IN)        :: t, rmin, rmax, rmaxu, qmin, qmin_alo, qmin_ahi, j2rp2, j4rp4
               REAL(DP), INTENT(INOUT)     :: eoffset
               CHARACTER(*), INTENT(IN)    :: qmin_coord
               TYPE(symba_pl), INTENT(INOUT)     :: symba_plA
          END SUBROUTINE symba_discard_pl
     END INTERFACE

     INTERFACE
          SUBROUTINE symba_discard_sun_pl(t, npl, msys, swiftest_plA, rmin, rmax, rmaxu, ldiscards)
               USE swiftest
               USE module_swiftest
               IMPLICIT NONE
               LOGICAL(LGT), INTENT(INOUT) :: ldiscards
               INTEGER(I4B), INTENT(IN)    :: npl
               REAL(DP), INTENT(IN)        :: t, msys, rmin, rmax, rmaxu
               TYPE(swiftest_pl), INTENT(INOUT)   :: swiftest_plA
          END SUBROUTINE symba_discard_sun_pl
     END INTERFACE

     INTERFACE
          SUBROUTINE symba_discard_tp(t, npl, ntp, nsp, symba_plA, symba_tpA, dt, &
               rmin, rmax, rmaxu, qmin, qmin_coord, qmin_alo, qmin_ahi, lclose, lrhill_present)
               USE swiftest
               USE module_swiftest
               USE module_helio
               USE module_symba
               IMPLICIT NONE
               LOGICAL(LGT), INTENT(IN)    :: lclose, lrhill_present
               INTEGER(I4B), INTENT(IN)    :: npl
               INTEGER(I4B), INTENT(INOUT) :: ntp, nsp
               REAL(DP), INTENT(IN)        :: t, dt, rmin, rmax, rmaxu, qmin, qmin_alo, qmin_ahi
               CHARACTER(*), INTENT(IN)    :: qmin_coord
               TYPE(symba_pl), INTENT(INOUT)     :: symba_plA
               TYPE(symba_tp), INTENT(INOUT)     :: symba_tpA
          END SUBROUTINE symba_discard_tp
     END INTERFACE

     INTERFACE
          SUBROUTINE symba_energy(npl, swiftest_plA, j2rp2, j4rp4, ke, pe, te, htot)
               USE swiftest
               USE module_swiftest
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN)               :: npl
               REAL(DP), INTENT(IN)                   :: j2rp2, j4rp4
               REAL(DP), INTENT(OUT)                  :: ke, pe, te
               REAL(DP), DIMENSION(NDIM), INTENT(OUT) :: htot
               TYPE(swiftest_pl), INTENT(INOUT)             :: swiftest_plA
          END SUBROUTINE symba_energy
     END INTERFACE

     INTERFACE
          SUBROUTINE symba_fragmentation(t, dt, index_enc, nmergeadd, nmergesub, mergeadd_list, &
               mergesub_list, eoffset, vbs, encounter_file, out_type, npl, ntp, &
               symba_plA, symba_tpA, nplplenc, npltpenc, pltpenc_list, plplenc_list, &
               nplmax, ntpmax, fragmax)
               USE swiftest
               USE module_swiftest
               USE module_helio
               USE module_symba
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN)                         :: index_enc, nplmax, ntpmax
               INTEGER(I4B), INTENT(INOUT)                      :: nmergeadd, nmergesub, nplplenc, npltpenc, fragmax
               INTEGER(I4B), INTENT(INOUT)                      :: npl, ntp
               REAL(DP), INTENT(IN)                             :: t, dt
               REAL(DP), INTENT(INOUT)                          :: eoffset
               REAL(DP), DIMENSION(NDIM), INTENT(IN)            :: vbs
               CHARACTER(*), INTENT(IN)                         :: encounter_file, out_type
               TYPE(symba_plplenc), INTENT(INOUT)               :: plplenc_list
               TYPE(symba_pltpenc), INTENT(INOUT)               :: pltpenc_list
               TYPE(symba_merger), INTENT(INOUT)                :: mergeadd_list, mergesub_list
               TYPE(symba_pl), INTENT(INOUT)                    :: symba_plA
               TYPE(symba_tp), INTENT(INOUT)                    :: symba_tpA
          END SUBROUTINE symba_fragmentation
     END INTERFACE

     INTERFACE
          SUBROUTINE symba_getacch(lextra_force, t, npl, nplm, symba_plA, j2rp2, j4rp4, nplplenc, &
               plplenc_list)
               USE swiftest
               USE module_swiftest
               USE module_helio
               USE module_symba
               IMPLICIT NONE
               LOGICAL(LGT), INTENT(IN)                      :: lextra_force
               INTEGER(I4B), INTENT(IN)                      :: npl, nplm, nplplenc
               REAL(DP), INTENT(IN)                          :: t, j2rp2, j4rp4
               TYPE(symba_pl), INTENT(INOUT)                 :: symba_plA
               TYPE(symba_plplenc), INTENT(IN)               :: plplenc_list
          END SUBROUTINE symba_getacch
     END INTERFACE

     INTERFACE
          SUBROUTINE symba_getacch_tp(lextra_force, t, npl, nplm, nplmax, ntp, ntpmax, symba_plA, symba_tpA, &
               xh, j2rp2, j4rp4,  &
               npltpenc, pltpenc_list)
               USE swiftest
               USE module_swiftest
               USE module_helio
               USE module_symba
               IMPLICIT NONE
               LOGICAL(LGT), INTENT(IN)                      :: lextra_force
               INTEGER(I4B), INTENT(IN)                      :: npl, nplm, nplmax, ntp, ntpmax, npltpenc
               REAL(DP), INTENT(IN)                          :: t, j2rp2, j4rp4
               REAL(DP), DIMENSION(NDIM, npl), INTENT(IN)    :: xh
               TYPE(symba_pl), INTENT(INOUT)                 :: symba_plA
               TYPE(symba_tp), INTENT(INOUT)                 :: symba_tpA
               TYPE(symba_pltpenc), INTENT(IN)               :: pltpenc_list
          END SUBROUTINE symba_getacch_tp
     END INTERFACE

     INTERFACE
          SUBROUTINE symba_helio_drift(irec, npl, symba_plA, dt)
               USE swiftest
               USE module_swiftest
               USE module_helio
               USE module_symba
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN)      :: irec, npl
               REAL(DP), INTENT(IN)          :: dt
               TYPE(symba_pl), INTENT(INOUT) :: symba_plA
          END SUBROUTINE symba_helio_drift
     END INTERFACE

     INTERFACE
          SUBROUTINE symba_helio_drift_tp(irec, ntp, symba_tpA, mu, dt)
               USE swiftest
               USE module_swiftest
               USE module_helio
               USE module_symba
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN)      :: irec, ntp
               REAL(DP), INTENT(IN)          :: mu, dt
               TYPE(symba_tp), INTENT(INOUT) :: symba_tpA
          END SUBROUTINE symba_helio_drift_tp
     END INTERFACE

     INTERFACE
          SUBROUTINE symba_helio_getacch(lflag, lextra_force, t, npl, nplm, nplmax, helio_plA, j2rp2, j4rp4)
               USE swiftest
               USE module_swiftest
               USE module_helio
               USE module_symba
               IMPLICIT NONE
               LOGICAL(LGT), INTENT(IN)       :: lflag, lextra_force
               INTEGER(I4B), INTENT(IN)       :: npl, nplm, nplmax
               REAL(DP), INTENT(IN)           :: t, j2rp2, j4rp4
               TYPE(helio_pl), INTENT(INOUT)  :: helio_plA
          END SUBROUTINE symba_helio_getacch
     END INTERFACE

     INTERFACE
          SUBROUTINE symba_helio_getacch_int(npl, nplm, helio_plA)
               USE swiftest
               USE module_helio
               USE module_symba
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN)      :: npl, nplm
               TYPE(helio_pl), INTENT(INOUT) :: helio_plA
          END SUBROUTINE symba_helio_getacch_int
     END INTERFACE

     INTERFACE
          SUBROUTINE symba_kick(irec, nplplenc, npltpenc, plplenc_list, pltpenc_list, dt, sgn, symba_plA, &
               symba_tpA)
               USE swiftest
               USE module_swiftest
               USE module_helio
               USE module_symba
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN)        :: irec, nplplenc, npltpenc
               REAL(DP), INTENT(IN)            :: dt, sgn
               TYPE(symba_plplenc), INTENT(IN) :: plplenc_list
               TYPE(symba_pltpenc), INTENT(IN) :: pltpenc_list
               TYPE(symba_pl), INTENT(INOUT)   :: symba_plA
               TYPE(symba_tp), INTENT(INOUT)   :: symba_tpA
          END SUBROUTINE symba_kick
     END INTERFACE

     INTERFACE
          SUBROUTINE symba_merge_pl(t, dt, index_enc, nplplenc, plplenc_list, nmergeadd, nmergesub, &
               mergeadd_list, mergesub_list, eoffset, vbs, encounter_file, out_type, npl, symba_plA, &
               symba_tpA)
               USE swiftest
               USE module_swiftest
               USE module_helio
               USE module_symba
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN)                         :: index_enc, nplplenc
               INTEGER(I4B), INTENT(INOUT)                      :: nmergeadd, nmergesub, npl
               REAL(DP), INTENT(IN)                             :: t, dt
               REAL(DP), INTENT(INOUT)                          :: eoffset
               REAL(DP), DIMENSION(NDIM), INTENT(IN)            :: vbs
               CHARACTER(*), INTENT(IN)                         :: encounter_file, out_type
               TYPE(symba_plplenc), INTENT(INOUT)               :: plplenc_list
               TYPE(symba_merger),  INTENT(INOUT)               :: mergeadd_list, mergesub_list
               TYPE(symba_pl), INTENT(INOUT)                    :: symba_plA
               TYPE(symba_tp), INTENT(INOUT)                    :: symba_tpA
          END SUBROUTINE symba_merge_pl
     END INTERFACE

     INTERFACE
          SUBROUTINE symba_merge_tp(t, dt, index_enc, npltpenc, pltpenc_list, vbs, encounter_file, out_type, symba_plA, symba_tpA)
               USE swiftest
               USE module_swifter
               USE module_helio
               USE module_symba
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN)               :: index_enc, npltpenc
               REAL(DP), INTENT(IN)                   :: t, dt
               REAL(DP), DIMENSION(NDIM), INTENT(IN)  :: vbs
               CHARACTER(*), INTENT(IN)               :: encounter_file, out_type
               TYPE(symba_pltpenc), INTENT(INOUT)     :: pltpenc_list
               TYPE(symba_pl), INTENT(INOUT)                    :: symba_plA
               TYPE(symba_tp), INTENT(INOUT)                    :: symba_tpA
          END SUBROUTINE symba_merge_tp
     END INTERFACE

     INTERFACE
          SUBROUTINE symba_peri(lfirst, npl, symba_plA, msys, qmin_coord)
               USE swiftest
               USE module_swiftest
               USE module_helio
               USE module_symba
               IMPLICIT NONE
               LOGICAL(LGT), INTENT(IN)       :: lfirst
               INTEGER(I4B), INTENT(IN)       :: npl
               REAL(DP), INTENT(IN)           :: msys
               CHARACTER(*), INTENT(IN)       :: qmin_coord
               TYPE(symba_pl), INTENT(INOUT)  :: symba_plA
          END SUBROUTINE symba_peri
     END INTERFACE

     INTERFACE
          SUBROUTINE symba_rearray(npl, ntp, nsppl, nsptp, symba_plA, symba_tpA, nmergeadd, mergeadd_list, discard_plA, &
    discard_tpA,feature)
               USE swiftest
               USE module_swiftestalloc 
               USE module_swiftest
               USE module_helio
               USE module_symba
               IMPLICIT NONE
               INTEGER(I4B), INTENT(INOUT)                      :: npl, ntp, nsppl, nsptp, nmergeadd !change to fragadd
               TYPE(symba_pl), INTENT(INOUT)                    :: symba_plA
               TYPE(symba_tp), INTENT(INOUT)                    :: symba_tpA
               TYPE(swiftest_tp), INTENT(INOUT)                 :: discard_tpA
               TYPE(swiftest_pl), INTENT(INOUT)                 :: discard_plA
               TYPE(symba_merger), INTENT(INOUT)                :: mergeadd_list !change to fragadd_list
               TYPE(feature_list),intent(in)                    :: feature

          END SUBROUTINE symba_rearray
     END INTERFACE  

     INTERFACE
          SUBROUTINE symba_reorder_pl(npl, symba_plA)
               USE swiftest
               USE module_swiftest
               USE module_helio
               USE module_symba
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN) :: npl
               TYPE(symba_pl), INTENT(INOUT)  :: symba_plA
               INTEGER(I4B)                              :: i
               INTEGER(I4B), DIMENSION(:), ALLOCATABLE   :: index
               REAL(DP), DIMENSION(:), ALLOCATABLE       :: mass
               REAL(DP), DIMENSION(:,:), allocatable     :: symba_plwkspA
               INTEGER(I4B), DIMENSION(:,:), allocatable :: symba_plwkspA_id_status
          END SUBROUTINE symba_reorder_pl
     END INTERFACE

     INTERFACE
          SUBROUTINE symba_setup(npl, ntp, symba_plA, symba_tpA, symba_pl1P, symba_tp1P, swiftest_pl1P, &
               swiftest_tp1P)
               USE swiftest
               USE module_swiftest
               USE module_helio
               USE module_symba
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN)                            :: npl, ntp
               TYPE(swiftest_pl), POINTER                           :: swiftest_pl1P
               TYPE(swiftest_tp), POINTER                           :: swiftest_tp1P
               TYPE(symba_pl), INTENT(INOUT) :: symba_plA
               TYPE(symba_tp), INTENT(INOUT) :: symba_tpA
               TYPE(symba_pl), POINTER                             :: symba_pl1P
               TYPE(symba_tp), POINTER                             :: symba_tp1P
          END SUBROUTINE symba_setup
     END INTERFACE

     INTERFACE
          SUBROUTINE symba_step(lfirst, lextra_force, lclose, t, npl, nplmax, ntp, ntpmax, symba_plA, &
               symba_tpA, j2rp2, j4rp4, dt, nplplenc, npltpenc, plplenc_list, pltpenc_list, nmergeadd, &
               nmergesub, mergeadd_list, mergesub_list, eoffset, mtiny, encounter_file, out_type, &
               fragmax, feature)
               USE swiftest
               USE module_swiftest
               USE module_helio
               USE module_symba
               IMPLICIT NONE
               LOGICAL(LGT), INTENT(IN)           :: lextra_force, lclose
               LOGICAL(LGT), INTENT(INOUT)        :: lfirst
               INTEGER(I4B), INTENT(IN)           :: npl, nplmax, ntp, ntpmax
               INTEGER(I4B), INTENT(INOUT)        :: nplplenc, npltpenc, nmergeadd, nmergesub, fragmax
               REAL(DP), INTENT(IN)               :: t, j2rp2, j4rp4, dt, mtiny
               REAL(DP), INTENT(INOUT)            :: eoffset
               CHARACTER(*), INTENT(IN)           :: encounter_file, out_type
               TYPE(symba_pl), INTENT(INOUT)      :: symba_plA
               TYPE(symba_tp), INTENT(INOUT)      :: symba_tpA
               TYPE(symba_plplenc), INTENT(INOUT) :: plplenc_list
               TYPE(symba_pltpenc), INTENT(INOUT) :: pltpenc_list
               TYPE(symba_merger), INTENT(INOUT)  :: mergeadd_list, mergesub_list
               TYPE(feature_list):: feature        ! Derived type containing logical flags to turn on or off various features of the code
          END SUBROUTINE symba_step
     END INTERFACE

     INTERFACE
          SUBROUTINE symba_step_helio(lfirst, lextra_force, t, npl, nplm, nplmax, ntp, ntpmax, helio_plA, helio_tpA, j2rp2,     &
               j4rp4, dt)
               USE swiftest
               USE module_swiftest
               USE module_helio
               USE module_symba
               IMPLICIT NONE
               LOGICAL(LGT), INTENT(IN)      :: lextra_force
               LOGICAL(LGT), INTENT(INOUT)   :: lfirst
               INTEGER(I4B), INTENT(IN)      :: npl, nplm, nplmax, ntp, ntpmax
               REAL(DP), INTENT(IN)          :: t, j2rp2, j4rp4, dt
               TYPE(helio_pl), INTENT(INOUT) :: helio_plA
               TYPE(helio_tp), INTENT(INOUT) :: helio_tpA
          END SUBROUTINE symba_step_helio
     END INTERFACE

     INTERFACE
          SUBROUTINE symba_step_helio_pl(lfirst, lextra_force, t, npl, nplm, nplmax, helio_plA, j2rp2, j4rp4, dt, xbeg, xend,    &
               ptb, pte)
               USE swiftest
               USE module_swiftest
               USE module_helio
               USE module_symba
               IMPLICIT NONE
               LOGICAL(LGT), INTENT(IN)                     :: lextra_force
               LOGICAL(LGT), INTENT(INOUT)                  :: lfirst
               INTEGER(I4B), INTENT(IN)                     :: npl, nplm, nplmax
               REAL(DP), INTENT(IN)                         :: t, j2rp2, j4rp4, dt
               REAL(DP), DIMENSION(NDIM, nplm), INTENT(OUT) :: xbeg, xend
               REAL(DP), DIMENSION(NDIM), INTENT(OUT)       :: ptb, pte
               TYPE(helio_pl), INTENT(INOUT)                :: helio_plA
          END SUBROUTINE symba_step_helio_pl
     END INTERFACE

     INTERFACE
          SUBROUTINE symba_step_interp(lextra_force, lclose, t, npl, nplm, nplmax, ntp, ntpmax, symba_plA, symba_tpA, j2rp2,    &
               j4rp4, dt, eoffset, nplplenc, npltpenc, plplenc_list, pltpenc_list, nmergeadd, nmergesub, mergeadd_list,    &
               mergesub_list, encounter_file, out_type, fragmax, feature)
               USE swiftest
               USE module_swiftest
               USE module_helio
               USE module_symba
               IMPLICIT NONE
               LOGICAL(LGT), INTENT(IN)           :: lextra_force, lclose
               INTEGER(I4B), INTENT(IN)           :: npl, nplm, nplmax, ntp, ntpmax, nplplenc, npltpenc
               INTEGER(I4B), INTENT(INOUT)        :: nmergeadd, nmergesub, fragmax
               REAL(DP), INTENT(IN)               :: t, j2rp2, j4rp4, dt
               REAL(DP), INTENT(INOUT)            :: eoffset
               CHARACTER(*), INTENT(IN)           :: encounter_file, out_type
               TYPE(symba_pl), INTENT(INOUT)      :: symba_plA
               TYPE(symba_tp), INTENT(INOUT)      :: symba_tpA
               TYPE(symba_plplenc), INTENT(INOUT) :: plplenc_list
               TYPE(symba_pltpenc), INTENT(INOUT) :: pltpenc_list
               TYPE(symba_merger), INTENT(INOUT)  :: mergeadd_list, mergesub_list
               type(feature_list), intent(in)     :: feature
          END SUBROUTINE symba_step_interp
     END INTERFACE

     INTERFACE
          RECURSIVE SUBROUTINE symba_step_recur(lclose, t, ireci, npl, nplm, ntp, symba_plA, symba_tpA, dt0, eoffset, nplplenc, &
               npltpenc, plplenc_list, pltpenc_list, nmergeadd, nmergesub, mergeadd_list, mergesub_list, encounter_file, & 
               out_type, nplmax, ntpmax, fragmax, feature)
               USE swiftest
               USE module_swiftest
               USE module_helio
               USE module_symba
               IMPLICIT NONE
               LOGICAL(LGT), INTENT(IN)           :: lclose
               INTEGER(I4B), INTENT(IN)           :: ireci, npl, nplm, ntp, nplplenc, npltpenc, nplmax, ntpmax, fragmax
               INTEGER(I4B), INTENT(INOUT)        :: nmergeadd, nmergesub
               REAL(DP), INTENT(IN)               :: t, dt0
               REAL(DP), INTENT(INOUT)            :: eoffset
               CHARACTER(*), INTENT(IN)           :: encounter_file, out_type
               TYPE(symba_pl), INTENT(INOUT)      :: symba_plA
               TYPE(symba_tp), INTENT(INOUT)      :: symba_tpA
               TYPE(symba_plplenc), INTENT(INOUT) :: plplenc_list
               TYPE(symba_pltpenc), INTENT(INOUT) :: pltpenc_list
               TYPE(symba_merger), INTENT(INOUT)  :: mergeadd_list, mergesub_list
               type(feature_list), intent(in)     :: feature
          END SUBROUTINE symba_step_recur
     END INTERFACE

     INTERFACE
          SUBROUTINE symba_user_getacch(t, npl, symba_plA)
               USE swiftest
               USE module_symba
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN)     :: npl
               REAL(DP), INTENT(IN)         :: t
               TYPE(symba_pl), INTENT(INOUT):: symba_plA
          END SUBROUTINE symba_user_getacch
     END INTERFACE

     INTERFACE
          SUBROUTINE symba_user_getacch_tp(t, ntp, symba_tpA)
               USE swiftest
               USE module_symba
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN)       :: ntp
               REAL(DP), INTENT(IN)           :: t
               TYPE(symba_tp), INTENT(INOUT)  :: symba_tpA
          END SUBROUTINE symba_user_getacch_tp
     END INTERFACE

     INTERFACE
          SUBROUTINE util_exit(code)
               USE swiftest
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN) :: code
          END SUBROUTINE util_exit
     END INTERFACE

     INTERFACE
          SUBROUTINE util_hills(npl, swiftest_plA)
               USE swiftest
               USE module_swiftest
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN)         :: npl
               TYPE(swiftest_pl), INTENT(INOUT) :: swiftest_plA
          END SUBROUTINE util_hills
     END INTERFACE

     INTERFACE
          SUBROUTINE util_index(arr, index)
               USE swiftest
               USE module_nrutil
               IMPLICIT NONE
               INTEGER(I4B), DIMENSION(:), INTENT(OUT) :: index
               REAL(DP), DIMENSION(:), INTENT(IN)      :: arr
          END SUBROUTINE util_index
     END INTERFACE

     INTERFACE
          SUBROUTINE util_peri(lfirst, ntp, swiftest_tpA, mu, msys, qmin_coord)
               USE swiftest
               USE module_swiftest
               IMPLICIT NONE
               LOGICAL(LGT), INTENT(IN)         :: lfirst
               INTEGER(I4B), INTENT(IN)         :: ntp
               REAL(DP), INTENT(IN)             :: mu, msys
               CHARACTER(*), INTENT(IN)         :: qmin_coord
               TYPE(swiftest_tp), INTENT(INOUT) :: swiftest_tpA
          END SUBROUTINE util_peri
     END INTERFACE

     INTERFACE 
          SUBROUTINE util_resize_pl(symba_plA, npl_new, npl_old)
               USE swiftest
               USE module_symba
               USE module_swiftest
               USE module_helio
               USE module_nrutil
               USE module_swiftestalloc
               IMPLICIT NONE
               TYPE(symba_pl), INTENT(INOUT) :: symba_plA
               INTEGER(I4B), INTENT(IN)      :: npl_old, npl_new
          END SUBROUTINE util_resize_pl
     END INTERFACE

     INTERFACE util_sort
          SUBROUTINE util_sort_i4b(arr)
               USE swiftest
               IMPLICIT NONE
               INTEGER(I4B), DIMENSION(:), INTENT(INOUT) :: arr
          END SUBROUTINE util_sort_i4b
          SUBROUTINE util_sort_sp(arr)
               USE swiftest
               IMPLICIT NONE
               REAL(SP), DIMENSION(:), INTENT(INOUT) :: arr
          END SUBROUTINE util_sort_sp
          SUBROUTINE util_sort_dp(arr)
               USE swiftest
               IMPLICIT NONE
               REAL(DP), DIMENSION(:), INTENT(INOUT) :: arr
          END SUBROUTINE util_sort_dp
     END INTERFACE

     INTERFACE
          SUBROUTINE util_toupper(string)
               USE swiftest
               IMPLICIT NONE
               CHARACTER(*), INTENT(INOUT) :: string
          END SUBROUTINE util_toupper
     END INTERFACE

     INTERFACE
          SUBROUTINE util_valid(npl, ntp, swiftest_plA, swiftest_tpA)
               USE swiftest
               USE module_swiftest
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN)         :: npl, ntp
               TYPE(swiftest_pl), INTENT(INOUT) :: swiftest_plA
               TYPE(swiftest_tp), INTENT(INOUT) :: swiftest_tpA
          END SUBROUTINE util_valid
     END INTERFACE

     INTERFACE
          SUBROUTINE util_version
               USE swiftest
               IMPLICIT NONE
          END SUBROUTINE util_version
     END INTERFACE

     ! Added by D. Minton
     INTERFACE
         FUNCTION util_kahan_sum(xsum_current, xi, xerror) 
            USE swiftest
            IMPLICIT NONE
            REAL(DP)                :: util_kahan_sum
            REAL(DP), INTENT(IN)    :: xsum_current, xi
            REAL(DP), INTENT(INOUT) :: xerror
         END FUNCTION
     END INTERFACE

     INTERFACE
         FUNCTION collresolve_resolve(model,m1,m2,r1,r2,p1,p2,v1,v2,n,mres,rres,pres,vres)
         USE swiftest
         IMPLICIT NONE
         INTEGER(I4B) :: collresolve_resolve
         INTEGER(I4B), INTENT(IN) :: model               ! collision model to apply
         REAL(DP), INTENT(IN) :: m1                      ! mass of the target
         REAL(DP), INTENT(IN) :: m2                      ! mass of the impactor
         REAL(DP), INTENT(IN) :: r1                      ! radius of the target
         REAL(DP), INTENT(IN) :: r2                      ! radius of the impactor
         REAL(DP), DIMENSION(3), INTENT(IN) :: p1        ! position of the target
         REAL(DP), DIMENSION(3), INTENT(IN) :: p2        ! position of the impactor
         REAL(DP), DIMENSION(3), INTENT(IN) :: v1        ! velocity of the target
         REAL(DP), DIMENSION(3), INTENT(IN) :: v2        ! velocity of the impactor
         INTEGER, INTENT(IN) :: n                                ! number of bodies to return
         REAL(DP), DIMENSION(n+1), INTENT(OUT) :: mres   ! mass of the resulting bodies
         REAL(DP), DIMENSION(n+1), INTENT(OUT) :: rres   ! radius of the resulting bodies
         REAL(DP), DIMENSION(3,n+1), INTENT(OUT) :: pres ! position of the resulting bodies
         REAL(DP), DIMENSION(3,n+1), INTENT(OUT) :: vres ! velocity of the resulting bodies
         END FUNCTION
      END INTERFACE

     INTERFACE
         SUBROUTINE util_regime(symba_plA, index1, index2, regime, Mlr, Mslr)
         USE swiftest
         USE module_swiftest
         USE module_symba
         IMPLICIT NONE
         TYPE(symba_pl), INTENT(INOUT) :: symba_plA
         INTEGER(I4B), INTENT(IN)      :: index1, index2
         INTEGER(I4B), INTENT(OUT)     :: regime
         REAL(DP), INTENT(OUT)         :: Mlr, Mslr
         END SUBROUTINE util_regime
     END INTERFACE

END MODULE module_interfaces
!**********************************************************************************************************************************
!
!  Author(s)   : David E. Kaufmann
!
!  Revision Control System (RCS) Information
!
!  Source File : $RCSfile$
!  Full Path   : $Source$
!  Revision    : $Revision$
!  Date        : $Date$
!  Programmer  : $Author$
!  Locked By   : $Locker$
!  State       : $State$
!
!  Modification History:
!
!  $Log$
!**********************************************************************************************************************************
