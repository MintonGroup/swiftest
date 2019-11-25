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
          SUBROUTINE coord_b2h(npl, swiftest_plA)
               USE module_parameters
               USE module_swiftest
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN)  :: npl
               TYPE(swiftest_pl), INTENT(INOUT) :: swiftest_plA
          END SUBROUTINE coord_b2h
     END INTERFACE

     INTERFACE
          SUBROUTINE coord_b2h_tp(ntp, swiftest_tpA, swiftest_plA)
               USE module_parameters
               USE module_swiftest
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN)         :: ntp
               TYPE(swiftest_tp), INTENT(INOUT) :: swiftest_tpA
               TYPE(swiftest_pl), INTENT(INOUT) :: swiftest_plA
          END SUBROUTINE coord_b2h_tp
     END INTERFACE

     INTERFACE
          SUBROUTINE coord_h2b(npl, swiftest_plA, msys)
               USE module_parameters
               USE module_swiftest
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN)         :: npl
               REAL(DP), INTENT(OUT)            :: msys
               TYPE(swiftest_pl), INTENT(INOUT) :: swiftest_plA
          END SUBROUTINE coord_h2b
     END INTERFACE

     INTERFACE
          SUBROUTINE coord_h2b_tp(ntp, swiftest_tpA, swiftest_plA)
               USE module_parameters
               USE module_swiftest
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN)         :: ntp
               TYPE(swiftest_tp), INTENT(INOUT) :: swiftest_tpA
               TYPE(swiftest_pl), INTENT(INOUT) :: swiftest_plA
          END SUBROUTINE coord_h2b_tp
     END INTERFACE

     INTERFACE
          SUBROUTINE coord_vb2vh(npl, swiftest_plA)
               USE module_parameters
               USE module_swiftest
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN)         :: npl
               TYPE(swiftest_pl), INTENT(INOUT) :: swiftest_plA
          END SUBROUTINE coord_vb2vh
     END INTERFACE

     INTERFACE
          SUBROUTINE coord_vb2vh_tp(ntp, swiftest_tpA, vs)
               USE module_parameters
               USE module_swiftest
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN)              :: ntp
               REAL(DP), DIMENSION(NDIM), INTENT(IN) :: vs
               TYPE(swiftest_tp), INTENT(INOUT)      :: swiftest_tpA
          END SUBROUTINE coord_vb2vh_tp
     END INTERFACE

     INTERFACE
          SUBROUTINE coord_vh2vb(npl, swiftest_plA, msys)
               USE module_parameters
               USE module_swiftest
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN)        :: npl
               REAL(DP), INTENT(OUT)           :: msys
               TYPE(swiftest_pl),INTENT(INOUT) :: swiftest_plA
          END SUBROUTINE coord_vh2vb
     END INTERFACE

     INTERFACE
          SUBROUTINE coord_vh2vb_tp(ntp, swiftest_tpA, vs)
               USE module_parameters
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
               USE module_parameters
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
               USE module_parameters
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
               USE module_parameters
               IMPLICIT NONE
               INTEGER(I4B), INTENT(OUT)             :: iflag
               REAL(DP), INTENT(IN)                  :: dt, r2crit
               REAL(DP), DIMENSION(NDIM), INTENT(IN) :: dx, dv
               REAL(DP), INTENT(OUT)                 :: r2min
          END SUBROUTINE discard_pl_close
     END INTERFACE

     INTERFACE
          SUBROUTINE discard_pl(t, dt, npl, ntp, swiftest_plA, swiftest_tpA)
               USE module_parameters
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
               USE module_parameters
               USE module_swiftest
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN)  :: ntp
               REAL(DP), INTENT(IN)      :: t, msys, rmin, rmax, rmaxu
               TYPE(swiftest_tp), INTENT(INOUT) :: swifter_tpA
          END SUBROUTINE discard_sun
     END INTERFACE

     INTERFACE
          SUBROUTINE drift_dan(mu, x0, v0, dt0, iflag)
               USE module_parameters
               IMPLICIT NONE
               INTEGER(I4B), INTENT(OUT)                :: iflag
               REAL(DP), INTENT(IN)                     :: mu, dt0
               REAL(DP), DIMENSION(NDIM), INTENT(INOUT) :: x0, v0
          END SUBROUTINE drift_dan
     END INTERFACE

     INTERFACE
          SUBROUTINE drift_kepmd(dm, es, ec, x, s, c)
               USE module_parameters
               IMPLICIT NONE
               REAL(DP), INTENT(IN)  :: dm, es, ec
               REAL(DP), INTENT(OUT) :: x, s, c
          END SUBROUTINE drift_kepmd
     END INTERFACE

     INTERFACE
          SUBROUTINE drift_kepu(dt,r0,mu,alpha,u,fp,c1,c2,c3,iflag)
               USE module_parameters
               IMPLICIT NONE
               INTEGER(I4B), INTENT(OUT) :: iflag
               REAL(DP), INTENT(IN)      :: dt, r0, mu, alpha, u
               REAL(DP), INTENT(OUT)     :: fp, c1, c2, c3
          END SUBROUTINE drift_kepu
     END INTERFACE

     INTERFACE
          SUBROUTINE drift_kepu_fchk(dt, r0, mu, alpha, u, s, f)
               USE module_parameters
               IMPLICIT NONE
               REAL(DP), INTENT(IN)  :: dt, r0, mu, alpha, u, s
               REAL(DP), INTENT(OUT) :: f
          END SUBROUTINE drift_kepu_fchk
     END INTERFACE

     INTERFACE
          SUBROUTINE drift_kepu_guess(dt, r0, mu, alpha, u, s)
               USE module_parameters
               IMPLICIT NONE
               REAL(DP), INTENT(IN)  :: dt, r0, mu, alpha, u
               REAL(DP), INTENT(OUT) :: s
          END SUBROUTINE drift_kepu_guess
     END INTERFACE

     INTERFACE
          SUBROUTINE drift_kepu_lag(s, dt, r0, mu, alpha, u, fp, c1, c2, c3, iflag)
               USE module_parameters
               IMPLICIT NONE
               INTEGER(I4B), INTENT(OUT) :: iflag
               REAL(DP), INTENT(IN)      :: dt, r0, mu, alpha, u
               REAL(DP), INTENT(INOUT)   :: s
               REAL(DP), INTENT(OUT)     :: fp, c1, c2, c3
          END SUBROUTINE drift_kepu_lag
     END INTERFACE

     INTERFACE
          SUBROUTINE drift_kepu_new(s, dt, r0, mu, alpha, u, fp, c1, c2, c3, iflag)
               USE module_parameters
               IMPLICIT NONE
               INTEGER(I4B), INTENT(OUT) :: iflag
               REAL(DP), INTENT(IN)      :: dt, r0, mu, alpha, u
               REAL(DP), INTENT(INOUT)   :: s
               REAL(DP), INTENT(OUT)     :: fp, c1, c2, c3
          END SUBROUTINE drift_kepu_new
     END INTERFACE

     INTERFACE
          SUBROUTINE drift_kepu_p3solve(dt, r0, mu, alpha, u, s, iflag)
               USE module_parameters
               IMPLICIT NONE
               INTEGER(I4B), INTENT(OUT) :: iflag
               REAL(DP), INTENT(IN)      :: dt, r0, mu, alpha, u
               REAL(DP), INTENT(OUT)     :: s
          END SUBROUTINE drift_kepu_p3solve
     END INTERFACE

     INTERFACE
          SUBROUTINE drift_kepu_stumpff(x, c0, c1, c2, c3)
               USE module_parameters
               IMPLICIT NONE
               REAL(DP), INTENT(INOUT) :: x
               REAL(DP), INTENT(OUT)   :: c0, c1, c2, c3
          END SUBROUTINE drift_kepu_stumpff
     END INTERFACE

     INTERFACE
          SUBROUTINE drift_one(mu, x, v, dt, iflag)
               USE module_parameters
               IMPLICIT NONE
               INTEGER(I4B), INTENT(OUT)                :: iflag
               REAL(DP), INTENT(IN)                     :: mu, dt
               REAL(DP), DIMENSION(NDIM), INTENT(INOUT) :: x, v
          END SUBROUTINE drift_one
     END INTERFACE

     INTERFACE
          SUBROUTINE helio_discard(t, npl, ntp, nsp, helio_pl1P, helio_tp1P, helio_tpd1P, dt, rmin, rmax, rmaxu, qmin,            &
               qmin_coord, qmin_alo, qmin_ahi, lclose, lrhill_present)
               USE module_parameters
               USE module_swifter
               USE module_helio
               IMPLICIT NONE
               LOGICAL(LGT), INTENT(IN)    :: lclose, lrhill_present
               INTEGER(I4B), INTENT(IN)    :: npl
               INTEGER(I4B), INTENT(INOUT) :: ntp, nsp
               REAL(DP), INTENT(IN)        :: t, dt, rmin, rmax, rmaxu, qmin, qmin_alo, qmin_ahi
               CHARACTER(*), INTENT(IN)    :: qmin_coord
               TYPE(helio_pl), POINTER     :: helio_pl1P
               TYPE(helio_tp), POINTER     :: helio_tp1P, helio_tpd1P
          END SUBROUTINE helio_discard
     END INTERFACE

     INTERFACE
          SUBROUTINE helio_discard_spill(ntp, nsp, helio_tp1P, helio_tpd1P, helio_tpspP)
               USE module_parameters
               USE module_swifter
               USE module_helio
               IMPLICIT NONE
               INTEGER(I4B), INTENT(INOUT) :: ntp, nsp
               TYPE(helio_tp), POINTER     :: helio_tp1P, helio_tpd1P, helio_tpspP
          END SUBROUTINE helio_discard_spill
     END INTERFACE

     INTERFACE
          SUBROUTINE helio_drift(npl, swiftest_plA, dt)
               USE module_parameters
               USE module_swiftest
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN)         :: npl
               REAL(DP), INTENT(IN)             :: dt
               TYPE(swiftest_pl), INTENT(INOUT) :: swiftest_plA
          END SUBROUTINE helio_drift
     END INTERFACE

     INTERFACE
          SUBROUTINE helio_drift_tp(ntp, swiftest_tpA, mu, dt)
               USE module_parameters
               USE module_swiftest
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN)         :: ntp
               REAL(DP), INTENT(IN)             :: mu, dt
               TYPE(swiftest_tp), INTENT(INOUT) :: swiftest_tpA
          END SUBROUTINE helio_drift_tp
     END INTERFACE

     INTERFACE
          SUBROUTINE helio_getacch(lflag, lextra_force, t, npl, nplmax, helio_plA, j2rp2, j4rp4)
               USE module_parameters
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
               USE module_parameters
               USE module_helio
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN)       :: npl
               TYPE(helio_pl), INTENT(INOUT)  :: helio_plA
          END SUBROUTINE helio_getacch_int
     END INTERFACE

     INTERFACE
          SUBROUTINE helio_getacch_int_tp(npl, ntp, swiftest_plA, helio_tpA, xh)
               USE module_parameters
               USE module_swiftest
               USE module_helio
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN)                   :: npl, ntp
               REAL(DP), DIMENSION(NDIM, npl), INTENT(IN) :: xh
               TYPE(swiftest_pl), INTENT(INOUT)           :: swiftest_plA
               TYPE(helio_tp), INTENT(INOUT)              :: helio_tpA
          END SUBROUTINE helio_getacch_int_tp
     END INTERFACE

     INTERFACE
          SUBROUTINE helio_getacch_tp(lflag, lextra_force, t, npl, nplmax, ntp, ntpmax, helio_plA, helio_tpA, xh, j2rp2, j4rp4)
               USE module_parameters
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
               USE module_parameters
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
               USE module_parameters
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
               USE module_parameters
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
               USE module_parameters
               USE module_swiftest
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN)              :: ntp
               REAL(DP), INTENT(IN)                  :: dt
               REAL(DP), DIMENSION(NDIM), INTENT(IN) :: pt
               TYPE(swiftest_tp), INTENT(INOUT)      :: swiftest_tpA
          END SUBROUTINE helio_lindrift_tp
     END INTERFACE

     INTERFACE
          SUBROUTINE helio_setup(npl, ntp, helio_plA, helio_tpA, helio_pl1P, helio_tp1P, swiftest_pl1P, swiftest_tp1P)
               USE module_parameters
               USE module_swiftest
               USE module_helio
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN)                            :: npl, ntp
               TYPE(swiftest_pl), POINTER                           :: swiftest_pl1P
               TYPE(swiftest_tp), POINTER                           :: swiftest_tp1P
               TYPE(helio_pl), DIMENSION(:), TARGET, INTENT(INOUT) :: helio_plA
               TYPE(helio_tp), DIMENSION(:), TARGET, INTENT(INOUT) :: helio_tpA
               TYPE(helio_pl), POINTER                             :: helio_pl1P
               TYPE(helio_tp), POINTER                             :: helio_tp1P
          END SUBROUTINE helio_setup
     END INTERFACE

     INTERFACE
          SUBROUTINE helio_step(lfirst, lextra_force, t, npl, nplmax, ntp, ntpmax, helio_plA, helio_tpA, j2rp2, j4rp4, dt)
               USE module_parameters
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
               USE module_parameters
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
               USE module_parameters
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
               USE module_parameters
               USE module_helio
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN)      :: npl
               REAL(DP), INTENT(IN)          :: t
               TYPE(helio_pl), INTENT(INOUT) :: helio_plA
          END SUBROUTINE helio_user_getacch
     END INTERFACE

     INTERFACE
          SUBROUTINE helio_user_getacch_tp(t, ntp, helio_tpA)
               USE module_parameters
               USE module_helio
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN)      :: ntp
               REAL(DP), INTENT(IN)          :: t
               TYPE(helio_tp), INTENT(INOUT) :: helio_tpA
          END SUBROUTINE helio_user_getacch_tp
     END INTERFACE

     INTERFACE
          SUBROUTINE io_discard_write(t, npl, nsp, swifter_pl1P, swifter_tpd1P, fname, lbig_discard)
               USE module_parameters
               USE module_swifter
               IMPLICIT NONE
               LOGICAL(LGT), INTENT(IN)  :: lbig_discard
               INTEGER(I4B), INTENT(IN)  :: npl, nsp
               REAL(DP), INTENT(IN)      :: t
               CHARACTER(*), INTENT(IN)  :: fname
               TYPE(swifter_pl), POINTER :: swifter_pl1P
               TYPE(swifter_tp), POINTER :: swifter_tpd1P
          END SUBROUTINE io_discard_write
     END INTERFACE

     INTERFACE
          SUBROUTINE io_discard_write_symba(t, mtiny, npl, ntp, nsppl, nsptp, nmergeadd, nmergesub, symba_plA,discard_plA, &
                    discard_tpA, mergeadd_list, mergesub_list, fname, lbig_discard, discard_plA_id_status, &
                    discard_tpA_id_status)
               USE module_parameters
               USE module_swiftest
               USE module_symba
               IMPLICIT NONE
               LOGICAL(LGT), INTENT(IN)                       :: lbig_discard
               INTEGER(I4B), INTENT(IN)                       :: npl, ntp, nsppl, nsptp, nmergeadd, nmergesub
               REAL(DP), INTENT(IN)                           :: t, mtiny
               CHARACTER(*), INTENT(IN)                       :: fname
               TYPE(symba_pl), INTENT(INOUT)                  :: symba_plA
               REAL(DP), DIMENSION(8,npl), INTENT(IN)         :: discard_plA
               REAL(DP), DIMENSION(8,ntp), INTENT(IN)         :: discard_tpA
               TYPE(symba_merger), INTENT(INOUT)              :: mergeadd_list, mergesub_list
               INTEGER(I4B), DIMENSION(2,npl), INTENT(OUT)    :: discard_plA_id_status
               INTEGER(I4B), DIMENSION(2,ntp), INTENT(OUT)    :: discard_tpA_id_status
          END SUBROUTINE io_discard_write_symba
     END INTERFACE

     INTERFACE
          SUBROUTINE io_dump_param(nplmax, ntpmax, ntp, t, tstop, dt, in_type, istep_out, outfile, out_type, out_form,            &
               istep_dump, j2rp2, j4rp4, lclose, rmin, rmax, rmaxu, qmin, qmin_coord, qmin_alo, qmin_ahi, encounter_file,         &
               lextra_force, lbig_discard, lrhill_present)
               USE module_parameters
               IMPLICIT NONE
               LOGICAL(LGT), INTENT(IN) :: lclose, lextra_force, lbig_discard, lrhill_present
               INTEGER(I4B), INTENT(IN) :: nplmax, ntpmax, ntp, istep_out, istep_dump
               REAL(DP), INTENT(IN)     :: t, tstop, dt, j2rp2, j4rp4, rmin, rmax, rmaxu, qmin, qmin_alo, qmin_ahi
               CHARACTER(*), INTENT(IN) :: qmin_coord, encounter_file, in_type, outfile, out_type, out_form
          END SUBROUTINE io_dump_param
     END INTERFACE

     INTERFACE
          SUBROUTINE io_dump_pl(npl, swiftest_plA, lclose, lrhill_present)
               USE module_parameters
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
               USE module_parameters
               USE module_swiftest
               USE module_fxdr
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN)        :: ntp
               TYPE(swiftest_tp), INTENT(INOUT):: swiftest_tpA
          END SUBROUTINE io_dump_tp
     END INTERFACE

     INTERFACE
          SUBROUTINE io_getn(inplfile, intpfile, in_type, npl, nplmax, ntp, ntpmax)
               USE module_parameters
               USE module_fxdr
               IMPLICIT NONE
               INTEGER(I4B), INTENT(INOUT) :: nplmax, ntpmax
               INTEGER(I4B), INTENT(OUT)   :: npl, ntp
               CHARACTER(*), INTENT(IN)    :: inplfile, intpfile, in_type
          END SUBROUTINE io_getn
     END INTERFACE

     INTERFACE
          SUBROUTINE io_get_token(buffer, ilength, ifirst, ilast, ierr)
               USE module_parameters
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN)    :: ilength
               INTEGER(I4B), INTENT(INOUT) :: ifirst
               INTEGER(I4B), INTENT(OUT)   :: ilast, ierr
               CHARACTER(*), INTENT(IN)    :: buffer
          END SUBROUTINE io_get_token
     END INTERFACE

     INTERFACE
          SUBROUTINE io_init_param(inparfile, nplmax, ntpmax, t0, tstop, dt, inplfile, intpfile, in_type, istep_out, outfile,     &
               out_type, out_form, out_stat, istep_dump, j2rp2, j4rp4, lclose, rmin, rmax, rmaxu, qmin, qmin_coord, qmin_alo,     &
               qmin_ahi, encounter_file, lextra_force, lbig_discard, lrhill_present)
               USE module_parameters
               IMPLICIT NONE
               LOGICAL(LGT), INTENT(OUT) :: lclose, lextra_force, lbig_discard, lrhill_present
               INTEGER(I4B), INTENT(OUT) :: nplmax, ntpmax, istep_out, istep_dump
               REAL(DP), INTENT(OUT)     :: t0, tstop, dt, j2rp2, j4rp4, rmin, rmax, rmaxu, qmin, qmin_alo, qmin_ahi
               CHARACTER(*), INTENT(IN)  :: inparfile
               CHARACTER(*), INTENT(OUT) :: qmin_coord, encounter_file, inplfile, intpfile, in_type, outfile, out_type, out_form, &
                                            out_stat
          END SUBROUTINE io_init_param
     END INTERFACE

     INTERFACE
          SUBROUTINE io_init_pl(inplfile, in_type, lclose, lrhill_present, npl, symba_plA)
               USE module_parameters
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
               USE module_parameters
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
               USE module_parameters
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN)  :: iu
               INTEGER(I4B), INTENT(OUT) :: ierr
               CHARACTER(*), INTENT(IN)  :: fname, fopenstat, fmt
          END SUBROUTINE io_open
     END INTERFACE

     INTERFACE
          SUBROUTINE io_open_fxdr(fname, fopenstat, lflag, iu, ierr)
               USE module_parameters
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
               USE module_parameters
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
               USE module_parameters
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
               USE module_parameters
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
               USE module_parameters
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
               USE module_parameters
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
               USE module_parameters
               USE module_fxdr
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN) :: iu, npl, ntp, iout_form
               REAL(DP), INTENT(IN)     :: t
               CHARACTER(*), INTENT(IN) :: out_type
          END SUBROUTINE io_write_hdr
     END INTERFACE

     INTERFACE
          SUBROUTINE io_write_line(iu, name, d1, d2, d3, d4, d5, d6, out_type, MASS, RADIUS)
               USE module_parameters
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
               USE module_parameters
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
               USE module_parameters
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
               USE module_parameters
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
               USE module_parameters
               IMPLICIT NONE
               REAL(DP), INTENT(IN)  :: angle
               REAL(DP), INTENT(OUT) :: sx, cx
          END SUBROUTINE orbel_scget
     END INTERFACE

     INTERFACE
          SUBROUTINE orbel_xv2aeq(x, v, mu, a, e, q)
               USE module_parameters
               IMPLICIT NONE
               REAL(DP), INTENT(IN)                  :: mu
               REAL(DP), DIMENSION(NDIM), INTENT(IN) :: x, v
               REAL(DP), INTENT(OUT)                 :: a, e, q
          END SUBROUTINE orbel_xv2aeq
     END INTERFACE

     INTERFACE
          SUBROUTINE orbel_xv2aqt(x, v, mu, a, q, capm, tperi)
               USE module_parameters
               IMPLICIT NONE
               REAL(DP), INTENT(IN)                  :: mu
               REAL(DP), DIMENSION(NDIM), INTENT(IN) :: x, v
               REAL(DP), INTENT(OUT)                 :: a, q, capm, tperi
          END SUBROUTINE orbel_xv2aqt
     END INTERFACE

     INTERFACE
          SUBROUTINE orbel_xv2el(x, v, mu, a, e, inc, capom, omega, capm)
               USE module_parameters
               IMPLICIT NONE
               REAL(DP), INTENT(IN)                  :: mu
               REAL(DP), DIMENSION(NDIM), INTENT(IN) :: x, v
               REAL(DP), INTENT(OUT)                 :: a, e, inc, capom, omega, capm
          END SUBROUTINE orbel_xv2el
     END INTERFACE

     INTERFACE
          SUBROUTINE rmvs_chk_ind(xr, vr, dt, r2crit, iflag)
               USE module_parameters
               IMPLICIT NONE
               REAL(DP), INTENT(IN)                  :: dt, r2crit
               REAL(DP), DIMENSION(NDIM), INTENT(IN) :: xr, vr
               INTEGER(I4B), INTENT(OUT)             :: iflag
          END SUBROUTINE rmvs_chk_ind
     END INTERFACE

     INTERFACE
          SUBROUTINE symba_add(npl, mergeadd_list, nmergeadd, symba_pl1P, swifter_pl1P, mtiny)
               USE module_parameters
               USE module_swifter
               USE module_helio
               USE module_symba
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN)                             :: npl, nmergeadd
               REAL(DP), INTENT(IN)                                 :: mtiny
               TYPE(swifter_pl), POINTER                            :: swifter_pl1P
               TYPE(symba_merger), DIMENSION(:), TARGET, INTENT(IN) :: mergeadd_list
               TYPE(symba_pl), POINTER                              :: symba_pl1P
          END SUBROUTINE
     END INTERFACE

     INTERFACE
          SUBROUTINE symba_chk(xr, vr, rhill1, rhill2, dt, irec, lencounter, lvdotr)
               USE module_parameters
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
               USE module_parameters
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
               USE module_parameters
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
               USE module_parameters
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
          SUBROUTINE symba_discard_spill_pl(npl, nsp, symba_pld1P, symba_plspP)
               USE module_parameters
               USE module_swifter
               USE module_helio
               USE module_symba
               IMPLICIT NONE
               INTEGER(I4B), INTENT(INOUT) :: npl, nsp
               TYPE(symba_pl), POINTER     :: symba_pld1P, symba_plspP
          END SUBROUTINE symba_discard_spill_pl
     END INTERFACE

     INTERFACE
          SUBROUTINE symba_discard_spill_tp(ntp, nsp, symba_tp1P, symba_tpd1P, symba_tpspP)
               USE module_parameters
               USE module_swifter
               USE module_helio
               USE module_symba
               IMPLICIT NONE
               INTEGER(I4B), INTENT(INOUT) :: ntp, nsp
               TYPE(symba_tp), POINTER     :: symba_tp1P, symba_tpd1P, symba_tpspP
          END SUBROUTINE symba_discard_spill_tp
     END INTERFACE

     INTERFACE
          SUBROUTINE symba_discard_sun_pl(t, npl, msys, swiftest_plA, rmin, rmax, rmaxu, ldiscards)
               USE module_parameters
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
               USE module_parameters
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
          SUBROUTINE symba_energy(npl, nplmax, swiftest_plA, j2rp2, j4rp4, ke, pe, te, htot)
               USE module_parameters
               USE module_swiftest
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN)               :: npl, nplmax
               REAL(DP), INTENT(IN)                   :: j2rp2, j4rp4
               REAL(DP), INTENT(OUT)                  :: ke, pe, te
               REAL(DP), DIMENSION(NDIM), INTENT(OUT) :: htot
               TYPE(swiftest_pl), INTENT(INOUT)             :: swiftest_plA
          END SUBROUTINE symba_energy
     END INTERFACE

     INTERFACE
          SUBROUTINE symba_fragmentation(t, npl, nplmax, ntp, ntpmax, symba_plA, nplplenc, plplenc_list)
               USE module_parameters
               USE module_swiftest
               USE module_helio
               USE module_symba
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN)                      :: nplplenc
               INTEGER(I4B), INTENT(INOUT)                   :: npl, nplmax, ntp, ntpmax
               REAL(DP), INTENT(IN)                          :: t
               TYPE(symba_pl), INTENT(INOUT)                 :: symba_plA
               TYPE(symba_plplenc), INTENT(IN)               :: plplenc_list
          END SUBROUTINE symba_fragmentation
     END INTERFACE

     INTERFACE
          SUBROUTINE symba_fragmentation_pl(t, dt, index, nplplenc, plplenc_list, nmergeadd, nmergesub, & 
               mergeadd_list, mergesub_list, eoffset, vbs, encounter_file, out_type)
               USE module_parameters
               USE module_swiftest
               USE module_helio
               USE module_symba
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN)                         :: index, nplplenc
               INTEGER(I4B), INTENT(INOUT)                      :: nmergeadd, nmergesub
               REAL(DP), INTENT(IN)                             :: t, dt
               REAL(DP), INTENT(INOUT)                          :: eoffset
               REAL(DP), DIMENSION(NDIM), INTENT(IN)            :: vbs
               CHARACTER(*), INTENT(IN)                         :: encounter_file, out_type
               TYPE(symba_plplenc), INTENT(INOUT)               :: plplenc_list
               TYPE(symba_merger), INTENT(INOUT)                :: mergeadd_list, mergesub_list
          END SUBROUTINE symba_fragmentation_pl
     END INTERFACE

     INTERFACE
          SUBROUTINE symba_getacch(lextra_force, t, npl, nplm, nplmax, symba_plA, j2rp2, j4rp4, nplplenc, &
               plplenc_list)
               USE module_parameters
               USE module_swiftest
               USE module_helio
               USE module_symba
               IMPLICIT NONE
               LOGICAL(LGT), INTENT(IN)                      :: lextra_force
               INTEGER(I4B), INTENT(IN)                      :: npl, nplm, nplmax, nplplenc
               REAL(DP), INTENT(IN)                          :: t, j2rp2, j4rp4
               TYPE(symba_pl), INTENT(INOUT)                 :: symba_plA
               TYPE(symba_plplenc), INTENT(IN)               :: plplenc_list
          END SUBROUTINE symba_getacch
     END INTERFACE

     INTERFACE
          SUBROUTINE symba_getacch_tp(lextra_force, t, npl, nplm, nplmax, ntp, ntpmax, symba_plA, symba_tpA, &
               xh, j2rp2, j4rp4,  &
               npltpenc, pltpenc_list)
               USE module_parameters
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
               USE module_parameters
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
               USE module_parameters
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
               USE module_parameters
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
               USE module_parameters
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
               USE module_parameters
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
               mergeadd_list, mergesub_list, eoffset, vbs, encounter_file, out_type, npl)
               USE module_parameters
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
          END SUBROUTINE symba_merge_pl
     END INTERFACE

     INTERFACE
          SUBROUTINE symba_merge_tp(t, dt, index_enc, npltpenc, pltpenc_list, vbs, encounter_file, out_type)
               USE module_parameters
               USE module_swifter
               USE module_helio
               USE module_symba
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN)               :: index_enc, npltpenc
               REAL(DP), INTENT(IN)                   :: t, dt
               REAL(DP), DIMENSION(NDIM), INTENT(IN)  :: vbs
               CHARACTER(*), INTENT(IN)               :: encounter_file, out_type
               TYPE(symba_pltpenc), INTENT(INOUT)     :: pltpenc_list
          END SUBROUTINE symba_merge_tp
     END INTERFACE

     INTERFACE
          SUBROUTINE symba_peri(lfirst, npl, symba_plA, msys, qmin_coord)
               USE module_parameters
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
          SUBROUTINE symba_rearray(t, npl, ntp, nsppl, nsptp, symba_plA, symba_tpA, nmergeadd, &
               mergeadd_list, discard_plA, discard_tpA, discard_plA_id_status, & 
               discard_tpA_id_status, NPLMAX, j2rp2, j4rp4)
               USE module_parameters
               USE module_swiftestalloc 
               USE module_swiftest
               USE module_helio
               USE module_symba
               IMPLICIT NONE
               INTEGER(I4B), INTENT(INOUT)                       :: npl, ntp, nsppl, nsptp, nmergeadd
               INTEGER(I4B), INTENT(INOUT)                       :: NPLMAX
               REAL(DP), INTENT(IN)                              :: t, j2rp2, j4rp4
               TYPE(symba_pl), INTENT(INOUT)                     :: symba_plA
               TYPE(symba_tp), INTENT(INOUT)                     :: symba_tpA
               TYPE(symba_merger), INTENT(INOUT)                 :: mergeadd_list
               REAL(DP), DIMENSION(8,npl), INTENT(OUT)           :: discard_plA
               REAL(DP), DIMENSION(8,ntp), INTENT(OUT)           :: discard_tpA
               INTEGER(I4B), DIMENSION(2,npl), INTENT(OUT)       :: discard_plA_id_status
               INTEGER(I4B), DIMENSION(2,ntp), INTENT(OUT)       :: discard_tpA_id_status
          END SUBROUTINE symba_rearray
     END INTERFACE  

     INTERFACE
          SUBROUTINE symba_reorder_pl(npl, symba_pl1P)
               USE module_parameters
               USE module_swifter
               USE module_helio
               USE module_symba
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN) :: npl
               TYPE(symba_pl), POINTER  :: symba_pl1P
          END SUBROUTINE symba_reorder_pl
     END INTERFACE

     INTERFACE
          SUBROUTINE symba_setup(npl, ntp, symba_plA, symba_tpA, symba_pl1P, symba_tp1P, swiftest_pl1P, &
               swiftest_tp1P)
               USE module_parameters
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
               nmergesub, mergeadd_list, mergesub_list, eoffset, mtiny, encounter_file, out_type)
               USE module_parameters
               USE module_swiftest
               USE module_helio
               USE module_symba
               IMPLICIT NONE
               LOGICAL(LGT), INTENT(IN)           :: lextra_force, lclose
               LOGICAL(LGT), INTENT(INOUT)        :: lfirst
               INTEGER(I4B), INTENT(IN)           :: npl, nplmax, ntp, ntpmax
               INTEGER(I4B), INTENT(INOUT)        :: nplplenc, npltpenc, nmergeadd, nmergesub
               REAL(DP), INTENT(IN)               :: t, j2rp2, j4rp4, dt, mtiny
               REAL(DP), INTENT(INOUT)            :: eoffset
               CHARACTER(*), INTENT(IN)           :: encounter_file, out_type
               TYPE(symba_pl), INTENT(INOUT)      :: symba_plA
               TYPE(symba_tp), INTENT(INOUT)      :: symba_tpA
               TYPE(symba_plplenc), INTENT(INOUT) :: plplenc_list
               TYPE(symba_pltpenc), INTENT(INOUT) :: pltpenc_list
               TYPE(symba_merger), INTENT(INOUT)  :: mergeadd_list, mergesub_list
          END SUBROUTINE symba_step
     END INTERFACE


! FOR TESTING PURPOSES ONLY _ USE WITH SYMBA_STEP_TEST
     INTERFACE
          SUBROUTINE symba_step_test(lfirst, lextra_force, lclose, t, npl, nplmax, ntp, ntpmax, symba_plA, symba_tpA,  & 
               j2rp2, j4rp4, dt, nplplenc, npltpenc, plplenc_list, pltpenc_list, nmergeadd, nmergesub, mergeadd_list, & 
               mergesub_list, eoffset, mtiny, encounter_file, out_type)
               USE module_parameters
               USE module_swiftest
               USE module_helio
               USE module_symba
               IMPLICIT NONE
               LOGICAL(LGT), INTENT(IN)           :: lextra_force, lclose
               LOGICAL(LGT), INTENT(INOUT)        :: lfirst
               INTEGER(I4B), INTENT(IN)           :: npl, nplmax, ntp, ntpmax
               INTEGER(I4B), INTENT(INOUT)        :: nplplenc, npltpenc, nmergeadd, nmergesub
               REAL(DP), INTENT(IN)               :: t, j2rp2, j4rp4, dt, mtiny
               REAL(DP), INTENT(INOUT)            :: eoffset
               CHARACTER(*), INTENT(IN)           :: encounter_file, out_type
               TYPE(symba_pl), INTENT(INOUT)      :: symba_plA
               TYPE(symba_tp), INTENT(INOUT)      :: symba_tpA
               TYPE(symba_plplenc), INTENT(INOUT) :: plplenc_list
               TYPE(symba_pltpenc), INTENT(INOUT) :: pltpenc_list
               TYPE(symba_merger), INTENT(INOUT)  :: mergeadd_list, mergesub_list
          END SUBROUTINE symba_step_test
     END INTERFACE




     INTERFACE
          SUBROUTINE symba_step_helio(lfirst, lextra_force, t, npl, nplm, nplmax, ntp, ntpmax, helio_plA, helio_tpA, j2rp2,     &
               j4rp4, dt)
               USE module_parameters
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
               USE module_parameters
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
               j4rp4, dt, eoffset, mtiny, nplplenc, npltpenc, plplenc_list, pltpenc_list, nmergeadd, nmergesub, mergeadd_list,    &
               mergesub_list, encounter_file, out_type)
               USE module_parameters
               USE module_swiftest
               USE module_helio
               USE module_symba
               IMPLICIT NONE
               LOGICAL(LGT), INTENT(IN)           :: lextra_force, lclose
               INTEGER(I4B), INTENT(IN)           :: npl, nplm, nplmax, ntp, ntpmax, nplplenc, npltpenc
               INTEGER(I4B), INTENT(INOUT)        :: nmergeadd, nmergesub
               REAL(DP), INTENT(IN)               :: t, j2rp2, j4rp4, dt, mtiny
               REAL(DP), INTENT(INOUT)            :: eoffset
               CHARACTER(*), INTENT(IN)           :: encounter_file, out_type
               TYPE(symba_pl), INTENT(INOUT)      :: symba_plA
               TYPE(symba_tp), INTENT(INOUT)      :: symba_tpA
               TYPE(symba_plplenc), INTENT(INOUT) :: plplenc_list
               TYPE(symba_pltpenc), INTENT(INOUT) :: pltpenc_list
               TYPE(symba_merger), INTENT(INOUT)  :: mergeadd_list, mergesub_list
          END SUBROUTINE symba_step_interp
     END INTERFACE

     INTERFACE
          RECURSIVE SUBROUTINE symba_step_recur(lclose, t, ireci, npl, nplm, ntp, symba_plA, symba_tpA, dt0, eoffset, nplplenc, &
               npltpenc, plplenc_list, pltpenc_list, nmergeadd, nmergesub, mergeadd_list, mergesub_list, encounter_file, out_type)
               USE module_parameters
               USE module_swiftest
               USE module_helio
               USE module_symba
               IMPLICIT NONE
               LOGICAL(LGT), INTENT(IN)           :: lclose
               INTEGER(I4B), INTENT(IN)           :: ireci, npl, nplm, ntp, nplplenc, npltpenc
               INTEGER(I4B), INTENT(INOUT)        :: nmergeadd, nmergesub
               REAL(DP), INTENT(IN)               :: t, dt0
               REAL(DP), INTENT(INOUT)            :: eoffset
               CHARACTER(*), INTENT(IN)           :: encounter_file, out_type
               TYPE(symba_pl), INTENT(INOUT)      :: symba_plA
               TYPE(symba_tp), INTENT(INOUT)      :: symba_tpA
               TYPE(symba_plplenc), INTENT(INOUT) :: plplenc_list
               TYPE(symba_pltpenc), INTENT(INOUT) :: pltpenc_list
               TYPE(symba_merger), INTENT(INOUT)  :: mergeadd_list, mergesub_list
          END SUBROUTINE symba_step_recur
     END INTERFACE

     INTERFACE
          SUBROUTINE symba_user_getacch(t, npl, symba_plA)
               USE module_parameters
               USE module_symba
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN)     :: npl
               REAL(DP), INTENT(IN)         :: t
               TYPE(symba_pl), INTENT(INOUT):: symba_plA
          END SUBROUTINE symba_user_getacch
     END INTERFACE

     INTERFACE
          SUBROUTINE symba_user_getacch_tp(t, ntp, symba_tpA)
               USE module_parameters
               USE module_symba
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN)       :: ntp
               REAL(DP), INTENT(IN)           :: t
               TYPE(symba_tp), INTENT(INOUT)  :: symba_tpA
          END SUBROUTINE symba_user_getacch_tp
     END INTERFACE

     INTERFACE
          SUBROUTINE util_exit(code)
               USE module_parameters
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN) :: code
          END SUBROUTINE util_exit
     END INTERFACE

     INTERFACE
          SUBROUTINE util_hills(npl, swiftest_plA)
               USE module_parameters
               USE module_swiftest
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN)         :: npl
               TYPE(swiftest_pl), INTENT(INOUT) :: swiftest_plA
          END SUBROUTINE util_hills
     END INTERFACE

     INTERFACE
          SUBROUTINE util_index(arr, index)
               USE module_parameters
               USE module_nrutil
               IMPLICIT NONE
               INTEGER(I4B), DIMENSION(:), INTENT(OUT) :: index
               REAL(DP), DIMENSION(:), INTENT(IN)      :: arr
          END SUBROUTINE util_index
     END INTERFACE

     INTERFACE
          SUBROUTINE util_peri(lfirst, ntp, swiftest_tpA, mu, msys, qmin_coord)
               USE module_parameters
               USE module_swiftest
               IMPLICIT NONE
               LOGICAL(LGT), INTENT(IN)         :: lfirst
               INTEGER(I4B), INTENT(IN)         :: ntp
               REAL(DP), INTENT(IN)             :: mu, msys
               CHARACTER(*), INTENT(IN)         :: qmin_coord
               TYPE(swiftest_tp), INTENT(INOUT) :: swiftest_tpA
          END SUBROUTINE util_peri
     END INTERFACE

     INTERFACE util_sort
          SUBROUTINE util_sort_i4b(arr)
               USE module_parameters
               IMPLICIT NONE
               INTEGER(I4B), DIMENSION(:), INTENT(INOUT) :: arr
          END SUBROUTINE util_sort_i4b
          SUBROUTINE util_sort_sp(arr)
               USE module_parameters
               IMPLICIT NONE
               REAL(SP), DIMENSION(:), INTENT(INOUT) :: arr
          END SUBROUTINE util_sort_sp
          SUBROUTINE util_sort_dp(arr)
               USE module_parameters
               IMPLICIT NONE
               REAL(DP), DIMENSION(:), INTENT(INOUT) :: arr
          END SUBROUTINE util_sort_dp
     END INTERFACE

     INTERFACE
          SUBROUTINE util_toupper(string)
               USE module_parameters
               IMPLICIT NONE
               CHARACTER(*), INTENT(INOUT) :: string
          END SUBROUTINE util_toupper
     END INTERFACE

     INTERFACE
          SUBROUTINE util_valid(npl, ntp, swiftest_plA, swiftest_tpA)
               USE module_parameters
               USE module_swiftest
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN)         :: npl, ntp
               TYPE(swiftest_pl), INTENT(INOUT) :: swiftest_plA
               TYPE(swiftest_tp), INTENT(INOUT) :: swiftest_tpA
          END SUBROUTINE util_valid
     END INTERFACE

     INTERFACE
          SUBROUTINE util_version
               USE module_parameters
               IMPLICIT NONE
          END SUBROUTINE util_version
     END INTERFACE

     ! Added by D. Minton
     INTERFACE
         FUNCTION util_kahan_sum(xsum_current, xi, xerror) 
            USE module_parameters
            IMPLICIT NONE
            REAL(DP)                :: util_kahan_sum
            REAL(DP), INTENT(IN)    :: xsum_current, xi
            REAL(DP), INTENT(INOUT) :: xerror
         END FUNCTION
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
