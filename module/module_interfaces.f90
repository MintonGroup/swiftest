!**********************************************************************************************************************************
!
!  Unit Name   : module_interfaces
!  Unit Type   : module
!  Project     : Swifter
!  Package     : module
!  Language    : Fortran 90/95
!
!  Description : Definition of interfaces of subroutines and functions used in Swifter package
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
          SUBROUTINE coord_b2h(npl, swifter_pl1P)
               USE module_parameters
               USE module_swifter
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN)  :: npl
               TYPE(swifter_pl), POINTER :: swifter_pl1P
          END SUBROUTINE coord_b2h
     END INTERFACE

     INTERFACE
          SUBROUTINE coord_b2h_tp(ntp, swifter_tp1P, swifter_pl1P)
               USE module_parameters
               USE module_swifter
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN)  :: ntp
               TYPE(swifter_tp), POINTER :: swifter_tp1P
               TYPE(swifter_pl), POINTER :: swifter_pl1P
          END SUBROUTINE coord_b2h_tp
     END INTERFACE

     INTERFACE
          SUBROUTINE coord_h2b(npl, swifter_pl1P, msys)
               USE module_parameters
               USE module_swifter
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN)  :: npl
               REAL(DP), INTENT(OUT)     :: msys
               TYPE(swifter_pl), POINTER :: swifter_pl1P
          END SUBROUTINE coord_h2b
     END INTERFACE

     INTERFACE
          SUBROUTINE coord_h2b_tp(ntp, swifter_tp1P, swifter_pl1P)
               USE module_parameters
               USE module_swifter
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN)  :: ntp
               TYPE(swifter_tp), POINTER :: swifter_tp1P
               TYPE(swifter_pl), POINTER :: swifter_pl1P
          END SUBROUTINE coord_h2b_tp
     END INTERFACE

     INTERFACE
          SUBROUTINE coord_h2j(npl, whm_pl1P)
               USE module_parameters
               USE module_swifter
               USE module_whm
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN) :: npl
               TYPE(whm_pl), POINTER    :: whm_pl1P
          END SUBROUTINE coord_h2j
     END INTERFACE

     INTERFACE
          SUBROUTINE coord_j2h(npl, whm_pl1P)
               USE module_parameters
               USE module_swifter
               USE module_whm
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN) :: npl
               TYPE(whm_pl), POINTER    :: whm_pl1P
          END SUBROUTINE coord_j2h
     END INTERFACE

     INTERFACE
          SUBROUTINE coord_vb2vh(npl, swifter_pl1P)
               USE module_parameters
               USE module_swifter
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN)  :: npl
               TYPE(swifter_pl), POINTER :: swifter_pl1P
          END SUBROUTINE coord_vb2vh
     END INTERFACE

     INTERFACE
          SUBROUTINE coord_vb2vh_tp(ntp, swifter_tp1P, vs)
               USE module_parameters
               USE module_swifter
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN)              :: ntp
               REAL(DP), DIMENSION(NDIM), INTENT(IN) :: vs
               TYPE(swifter_tp), POINTER             :: swifter_tp1P
          END SUBROUTINE coord_vb2vh_tp
     END INTERFACE

     INTERFACE
          SUBROUTINE coord_vh2vb(npl, swifter_pl1P, msys)
               USE module_parameters
               USE module_swifter
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN)  :: npl
               REAL(DP), INTENT(OUT)     :: msys
               TYPE(swifter_pl), POINTER :: swifter_pl1P
          END SUBROUTINE coord_vh2vb
     END INTERFACE

     INTERFACE
          SUBROUTINE coord_vh2vb_tp(ntp, swifter_tp1P, vs)
               USE module_parameters
               USE module_swifter
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN)              :: ntp
               REAL(DP), DIMENSION(NDIM), INTENT(IN) :: vs
               TYPE(swifter_tp), POINTER             :: swifter_tp1P
          END SUBROUTINE coord_vh2vb_tp
     END INTERFACE

     INTERFACE
          SUBROUTINE coord_vh2vj(npl, whm_pl1P)
               USE module_parameters
               USE module_swifter
               USE module_whm
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN) :: npl
               TYPE(whm_pl), POINTER    :: whm_pl1P
          END SUBROUTINE coord_vh2vj
     END INTERFACE

     INTERFACE
          SUBROUTINE discard(t, dt, npl, ntp, swifter_pl1P, swifter_tp1P, rmin, rmax, rmaxu, qmin, qmin_alo, qmin_ahi,            &
               qmin_coord, lclose, lrhill_present)
               USE module_parameters
               USE module_swifter
               IMPLICIT NONE
               LOGICAL(LGT), INTENT(IN)  :: lclose, lrhill_present
               INTEGER(I4B), INTENT(IN)  :: npl, ntp
               REAL(DP), INTENT(IN)      :: t, dt, rmin, rmax, rmaxu, qmin, qmin_alo, qmin_ahi
               CHARACTER(*), INTENT(IN)  :: qmin_coord
               TYPE(swifter_pl), POINTER :: swifter_pl1P
               TYPE(swifter_tp), POINTER :: swifter_tp1P
          END SUBROUTINE discard
     END INTERFACE

     INTERFACE
          SUBROUTINE discard_peri(t, npl, ntp, swifter_pl1P, swifter_tp1P, msys, qmin, qmin_alo, qmin_ahi, qmin_coord,            &
               lrhill_present)
               USE module_parameters
               USE module_swifter
               IMPLICIT NONE
               LOGICAL(LGT), INTENT(IN)  :: lrhill_present
               INTEGER(I4B), INTENT(IN)  :: npl, ntp
               REAL(DP), INTENT(IN)      :: t, msys, qmin, qmin_alo, qmin_ahi
               CHARACTER(*), INTENT(IN)  :: qmin_coord
               TYPE(swifter_pl), POINTER :: swifter_pl1P
               TYPE(swifter_tp), POINTER :: swifter_tp1P
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
          SUBROUTINE discard_pl(t, dt, npl, ntp, swifter_pl1P, swifter_tp1P)
               USE module_parameters
               USE module_swifter
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN)  :: npl, ntp
               REAL(DP), INTENT(IN)      :: t, dt
               TYPE(swifter_pl), POINTER :: swifter_pl1P
               TYPE(swifter_tp), POINTER :: swifter_tp1P
          END SUBROUTINE discard_pl
     END INTERFACE

     INTERFACE
          SUBROUTINE discard_sun(t, ntp, msys, swifter_tp1P, rmin, rmax, rmaxu)
               USE module_parameters
               USE module_swifter
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN)  :: ntp
               REAL(DP), INTENT(IN)      :: t, msys, rmin, rmax, rmaxu
               TYPE(swifter_tp), POINTER :: swifter_tp1P
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
          SUBROUTINE helio_drift(npl, swifter_pl1P, dt)
               USE module_parameters
               USE module_swifter
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN)  :: npl
               REAL(DP), INTENT(IN)      :: dt
               TYPE(swifter_pl), POINTER :: swifter_pl1P
          END SUBROUTINE helio_drift
     END INTERFACE

     INTERFACE
          SUBROUTINE helio_drift_tp(ntp, swifter_tp1P, mu, dt)
               USE module_parameters
               USE module_swifter
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN)  :: ntp
               REAL(DP), INTENT(IN)      :: mu, dt
               TYPE(swifter_tp), POINTER :: swifter_tp1P
          END SUBROUTINE helio_drift_tp
     END INTERFACE

     INTERFACE
          SUBROUTINE helio_getacch(lflag, lextra_force, t, npl, nplmax, helio_pl1P, j2rp2, j4rp4)
               USE module_parameters
               USE module_swifter
               USE module_helio
               IMPLICIT NONE
               LOGICAL(LGT), INTENT(IN) :: lflag, lextra_force
               INTEGER(I4B), INTENT(IN) :: npl, nplmax
               REAL(DP), INTENT(IN)     :: t, j2rp2, j4rp4
               TYPE(helio_pl), POINTER  :: helio_pl1P
          END SUBROUTINE helio_getacch
     END INTERFACE

     INTERFACE
          SUBROUTINE helio_getacch_int(npl, helio_pl1P)
               USE module_parameters
               USE module_helio
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN) :: npl
               TYPE(helio_pl), POINTER  :: helio_pl1P
          END SUBROUTINE helio_getacch_int
     END INTERFACE

     INTERFACE
          SUBROUTINE helio_getacch_int_tp(npl, ntp, swifter_pl1P, helio_tp1P, xh)
               USE module_parameters
               USE module_swifter
               USE module_helio
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN)                   :: npl, ntp
               REAL(DP), DIMENSION(NDIM, npl), INTENT(IN) :: xh
               TYPE(swifter_pl), POINTER                  :: swifter_pl1P
               TYPE(helio_tp), POINTER                    :: helio_tp1P
          END SUBROUTINE helio_getacch_int_tp
     END INTERFACE

     INTERFACE
          SUBROUTINE helio_getacch_tp(lflag, lextra_force, t, npl, nplmax, ntp, ntpmax, helio_pl1P, helio_tp1P, xh, j2rp2, j4rp4)
               USE module_parameters
               USE module_swifter
               USE module_helio
               IMPLICIT NONE
               LOGICAL(LGT), INTENT(IN)                   :: lflag, lextra_force
               INTEGER(I4B), INTENT(IN)                   :: npl, nplmax, ntp, ntpmax
               REAL(DP), INTENT(IN)                       :: t, j2rp2, j4rp4
               REAL(DP), DIMENSION(NDIM, npl), INTENT(IN) :: xh
               TYPE(helio_pl), POINTER                    :: helio_pl1P
               TYPE(helio_tp), POINTER                    :: helio_tp1P
          END SUBROUTINE helio_getacch_tp
     END INTERFACE

     INTERFACE
          SUBROUTINE helio_kickvb(npl, helio_pl1P, dt)
               USE module_parameters
               USE module_swifter
               USE module_helio
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN) :: npl
               REAL(DP), INTENT(IN)     :: dt
               TYPE(helio_pl), POINTER  :: helio_pl1P
          END SUBROUTINE helio_kickvb
     END INTERFACE

     INTERFACE
          SUBROUTINE helio_kickvb_tp(ntp, helio_tp1P, dt)
               USE module_parameters
               USE module_swifter
               USE module_helio
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN) :: ntp
               REAL(DP), INTENT(IN)     :: dt
               TYPE(helio_tp), POINTER  :: helio_tp1P
          END SUBROUTINE helio_kickvb_tp
     END INTERFACE

     INTERFACE
          SUBROUTINE helio_lindrift(npl, swifter_pl1P, dt, pt)
               USE module_parameters
               USE module_swifter
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN)               :: npl
               REAL(DP), INTENT(IN)                   :: dt
               REAL(DP), DIMENSION(NDIM), INTENT(OUT) :: pt
               TYPE(swifter_pl), POINTER              :: swifter_pl1P
          END SUBROUTINE helio_lindrift
     END INTERFACE

     INTERFACE
          SUBROUTINE helio_lindrift_tp(ntp, swifter_tp1P, dt, pt)
               USE module_parameters
               USE module_swifter
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN)              :: ntp
               REAL(DP), INTENT(IN)                  :: dt
               REAL(DP), DIMENSION(NDIM), INTENT(IN) :: pt
               TYPE(swifter_tp), POINTER             :: swifter_tp1P
          END SUBROUTINE helio_lindrift_tp
     END INTERFACE

     INTERFACE
          SUBROUTINE helio_setup(npl, ntp, helio_plA, helio_tpA, helio_pl1P, helio_tp1P, swifter_pl1P, swifter_tp1P)
               USE module_parameters
               USE module_swifter
               USE module_helio
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN)                            :: npl, ntp
               TYPE(swifter_pl), POINTER                           :: swifter_pl1P
               TYPE(swifter_tp), POINTER                           :: swifter_tp1P
               TYPE(helio_pl), DIMENSION(:), TARGET, INTENT(INOUT) :: helio_plA
               TYPE(helio_tp), DIMENSION(:), TARGET, INTENT(INOUT) :: helio_tpA
               TYPE(helio_pl), POINTER                             :: helio_pl1P
               TYPE(helio_tp), POINTER                             :: helio_tp1P
          END SUBROUTINE helio_setup
     END INTERFACE

     INTERFACE
          SUBROUTINE helio_step(lfirst, lextra_force, t, npl, nplmax, ntp, ntpmax, helio_pl1P, helio_tp1P, j2rp2, j4rp4, dt)
               USE module_parameters
               USE module_helio
               IMPLICIT NONE
               LOGICAL(LGT), INTENT(IN)    :: lextra_force
               LOGICAL(LGT), INTENT(INOUT) :: lfirst
               INTEGER(I4B), INTENT(IN)    :: npl, nplmax, ntp, ntpmax
               REAL(DP), INTENT(IN)        :: t, j2rp2, j4rp4, dt
               TYPE(helio_pl), POINTER     :: helio_pl1P
               TYPE(helio_tp), POINTER     :: helio_tp1P
          END SUBROUTINE helio_step
     END INTERFACE

     INTERFACE
          SUBROUTINE helio_step_pl(lfirst, lextra_force, t, npl, nplmax, helio_pl1P, j2rp2, j4rp4, dt, xbeg, xend, ptb, pte)
               USE module_parameters
               USE module_swifter
               USE module_helio
               IMPLICIT NONE
               LOGICAL(LGT), INTENT(IN)                    :: lextra_force
               LOGICAL(LGT), INTENT(INOUT)                 :: lfirst
               INTEGER(I4B), INTENT(IN)                    :: npl, nplmax
               REAL(DP), INTENT(IN)                        :: t, j2rp2, j4rp4, dt
               REAL(DP), DIMENSION(NDIM), INTENT(OUT)      :: ptb, pte
               REAL(DP), DIMENSION(NDIM, npl), INTENT(OUT) :: xbeg,xend
               TYPE(helio_pl), POINTER                     :: helio_pl1P
          END SUBROUTINE helio_step_pl
     END INTERFACE

     INTERFACE
          SUBROUTINE helio_step_tp(lfirsttp, lextra_force, t, npl, nplmax, ntp, ntpmax, helio_pl1P, helio_tp1P, j2rp2, j4rp4, dt, &
               xbeg, xend, ptb, pte)
               USE module_parameters
               USE module_swifter
               USE module_helio
               IMPLICIT NONE
               LOGICAL(LGT), INTENT(IN)                   :: lextra_force
               LOGICAL(LGT), INTENT(INOUT)                :: lfirsttp
               INTEGER(I4B), INTENT(IN)                   :: npl, nplmax, ntp, ntpmax
               REAL(DP), INTENT(IN)                       :: t, j2rp2, j4rp4, dt
               REAL(DP), DIMENSION(NDIM), INTENT(IN)      :: ptb, pte
               REAL(DP), DIMENSION(NDIM, npl), INTENT(IN) :: xbeg, xend
               TYPE(helio_pl), POINTER                    :: helio_pl1P
               TYPE(helio_tp), POINTER                    :: helio_tp1P
          END SUBROUTINE helio_step_tp
     END INTERFACE

     INTERFACE
          SUBROUTINE helio_user_getacch(t, npl, helio_pl1P)
               USE module_parameters
               USE module_helio
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN) :: npl
               REAL(DP), INTENT(IN)     :: t
               TYPE(helio_pl), POINTER  :: helio_pl1P
          END SUBROUTINE helio_user_getacch
     END INTERFACE

     INTERFACE
          SUBROUTINE helio_user_getacch_tp(t, ntp, helio_tp1P)
               USE module_parameters
               USE module_helio
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN) :: ntp
               REAL(DP), INTENT(IN)     :: t
               TYPE(helio_tp), POINTER  :: helio_tp1P
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
          SUBROUTINE io_discard_write_symba(t, mtiny, npl, nsppl, nsptp, nmergeadd, nmergesub, symba_pl1P, symba_pld1P,           &
               symba_tpd1P, mergeadd_list, mergesub_list, fname, lbig_discard)
               USE module_parameters
               USE module_swifter
               USE module_symba
               IMPLICIT NONE
               LOGICAL(LGT), INTENT(IN)                     :: lbig_discard
               INTEGER(I4B), INTENT(IN)                     :: npl, nsppl, nsptp, nmergeadd, nmergesub
               REAL(DP), INTENT(IN)                         :: t, mtiny
               CHARACTER(*), INTENT(IN)                     :: fname
               TYPE(symba_pl), POINTER                      :: symba_pl1P, symba_pld1P
               TYPE(symba_tp), POINTER                      :: symba_tpd1P
               TYPE(symba_merger), DIMENSION(:), INTENT(IN) :: mergeadd_list, mergesub_list
          END SUBROUTINE io_discard_write_symba
     END INTERFACE

     INTERFACE
          SUBROUTINE io_getn(inplfile, intpfile, in_type, npl, nplmax, ntp, ntpmax)
               USE module_parameters
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
               out_type, out_form, out_stat, j2rp2, j4rp4, lclose, rmin, rmax, rmaxu, qmin, qmin_coord, qmin_alo,     &
               qmin_ahi, encounter_file, lextra_force, lbig_discard, lrhill_present, mtiny, ring_outfile)
               USE module_parameters
               IMPLICIT NONE
               LOGICAL(LGT), INTENT(OUT) :: lclose, lextra_force, lbig_discard, lrhill_present
               INTEGER(I4B), INTENT(OUT) :: nplmax, ntpmax, istep_out
               REAL(DP), INTENT(OUT)     :: t0, tstop, dt, j2rp2, j4rp4, rmin, rmax, rmaxu, qmin, qmin_alo, qmin_ahi
               CHARACTER(*), INTENT(IN)  :: inparfile
               CHARACTER(*), INTENT(OUT) :: qmin_coord, encounter_file, inplfile, intpfile, in_type, outfile, out_type, out_form, &
                                            out_stat
               REAL(DP), INTENT(OUT), OPTIONAL :: mtiny
               CHARACTER(*), INTENT(OUT), OPTIONAL  :: ring_outfile
          END SUBROUTINE io_init_param
     END INTERFACE

     INTERFACE
          SUBROUTINE io_init_pl(inplfile, in_type, lclose, lrhill_present, npl, swifter_pl1P)
               USE module_parameters
               USE module_swifter
               IMPLICIT NONE
               LOGICAL(LGT), INTENT(IN)  :: lclose, lrhill_present
               INTEGER(I4B), INTENT(IN)  :: npl
               CHARACTER(*), INTENT(IN)  :: inplfile, in_type
               TYPE(swifter_pl), POINTER :: swifter_pl1P
          END SUBROUTINE io_init_pl
     END INTERFACE

     INTERFACE
          SUBROUTINE io_init_tp(intpfile, in_type, ntp, swifter_tp1P)
               USE module_parameters
               USE module_swifter
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN)  :: ntp
               CHARACTER(*), INTENT(IN)  :: intpfile, in_type
               TYPE(swifter_tp), POINTER :: swifter_tp1P
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
          FUNCTION io_read_encounter(t, id1, id2, mass1, mass2, xh1, xh2, vh1, vh2, encounter_file, out_type)
               USE module_parameters
               IMPLICIT NONE
               INTEGER(I4B)                           :: io_read_encounter
               INTEGER(I4B), INTENT(OUT)              :: id1, id2
               REAL(DP), INTENT(OUT)                  :: t, mass1, mass2
               REAL(DP), DIMENSION(NDIM), INTENT(OUT) :: xh1, xh2, vh1, vh2
               CHARACTER(*), INTENT(IN)               :: encounter_file,out_type
          END FUNCTION io_read_encounter
     END INTERFACE

     INTERFACE
          FUNCTION io_read_hdr(iu, t, npl, ntp, iout_form, out_type)
               USE module_parameters
               IMPLICIT NONE
               INTEGER(I4B)               :: io_read_hdr
               INTEGER(I4B), INTENT(IN)   :: iu
               INTEGER(I4B), INTENT(OUT)  :: npl, ntp, iout_form
               REAL(DP), INTENT(OUT)      :: t
               CHARACTER(*), INTENT(IN)   :: out_type
          END FUNCTION io_read_hdr
     END INTERFACE

     INTERFACE
          FUNCTION io_read_line(iu, id, d1, d2, d3, d4, d5, d6, out_type, MASS, RADIUS)
               USE module_parameters
               IMPLICIT NONE
               INTEGER(I4B)                    :: io_read_line
               INTEGER(I4B), INTENT(IN)        :: iu
               INTEGER(I4B), INTENT(OUT)       :: id
               REAL(DP), INTENT(OUT)           :: d1, d2, d3, d4, d5, d6
               REAL(DP), OPTIONAL, INTENT(OUT) :: MASS, RADIUS
               CHARACTER(*), INTENT(IN)        :: out_type
          END FUNCTION io_read_line
     END INTERFACE

     INTERFACE
          SUBROUTINE io_write_encounter(t, id1, id2, mass1, mass2, radius1, radius2, xh1, xh2, vh1, vh2, encounter_file, out_type)
               USE module_parameters
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN)              :: id1, id2
               REAL(DP), INTENT(IN)                  :: t, mass1, mass2, radius1, radius2
               REAL(DP), DIMENSION(NDIM), INTENT(IN) :: xh1, xh2, vh1, vh2
               CHARACTER(*), INTENT(IN)              :: encounter_file, out_type
          END SUBROUTINE io_write_encounter
     END INTERFACE

     INTERFACE
          SUBROUTINE io_write_frame(t, npl, ntp, swifter_pl1P, swifter_tp1P, outfile, out_type, out_form, out_stat)
               USE module_parameters
               USE module_swifter
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN)  :: npl, ntp
               REAL(DP), INTENT(IN)      :: t
               CHARACTER(*), INTENT(IN)  :: outfile, out_type, out_form, out_stat
               TYPE(swifter_pl), POINTER :: swifter_pl1P
               TYPE(swifter_tp), POINTER :: swifter_tp1P
          END SUBROUTINE io_write_frame
     END INTERFACE

     INTERFACE
          SUBROUTINE io_write_hdr(iu, t, npl, ntp, iout_form, out_type)
               USE module_parameters
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN) :: iu, npl, ntp, iout_form
               REAL(DP), INTENT(IN)     :: t
               CHARACTER(*), INTENT(IN) :: out_type
          END SUBROUTINE io_write_hdr
     END INTERFACE

     INTERFACE
          SUBROUTINE io_write_line(iu, id, d1, d2, d3, d4, d5, d6, out_type, MASS, RADIUS)
               USE module_parameters
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN)       :: iu, id
               REAL(DP), INTENT(IN)           :: d1, d2, d3, d4, d5, d6
               REAL(DP), OPTIONAL, INTENT(IN) :: MASS, RADIUS
               CHARACTER(*), INTENT(IN)       :: out_type
          END SUBROUTINE io_write_line
     END INTERFACE

     INTERFACE
          SUBROUTINE obl_acc(npl, swifter_pl1P, j2rp2, j4rp4, xh, irh, aobl)
               USE module_parameters
               USE module_swifter
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN)                    :: npl
               REAL(DP), INTENT(IN)                        :: j2rp2, j4rp4
               REAL(DP), DIMENSION(npl), INTENT(IN)        :: irh
               REAL(DP), DIMENSION(NDIM, npl), INTENT(IN)  :: xh
               REAL(DP), DIMENSION(NDIM, npl), INTENT(OUT) :: aobl
               TYPE(swifter_pl), POINTER                   :: swifter_pl1P
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
          SUBROUTINE obl_pot(npl, swifter_pl1P, j2rp2, j4rp4, xh, irh, oblpot)
               USE module_parameters
               USE module_swifter
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN)                   :: npl
               REAL(DP), INTENT(IN)                       :: j2rp2, j4rp4
               REAL(DP), INTENT(OUT)                      :: oblpot
               REAL(DP), DIMENSION(npl), INTENT(IN)       :: irh
               REAL(DP), DIMENSION(NDIM, npl), INTENT(IN) :: xh
               TYPE(swifter_pl), POINTER                  :: swifter_pl1P
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
          SUBROUTINE rmvs_chk(npl, ntp, rmvs_pl1P, rmvs_tp1P, xh, vh, dt, rts, lencounter)
               USE module_parameters
               USE module_rmvs
               IMPLICIT NONE
               LOGICAL(LGT), INTENT(OUT)                  :: lencounter
               INTEGER(I4B), INTENT(IN)                   :: npl, ntp
               REAL(DP), INTENT(IN)                       :: dt, rts
               REAL(DP), DIMENSION(NDIM, npl), INTENT(IN) :: xh,vh
               TYPE(rmvs_pl), POINTER                     :: rmvs_pl1P
               TYPE(rmvs_tp), POINTER                     :: rmvs_tp1P
          END SUBROUTINE rmvs_chk
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
          SUBROUTINE rmvs_discard(t, npl, ntp, nsp, rmvs_pl1P, rmvs_tp1P, rmvs_tpd1P, dt, rmin, rmax, rmaxu, qmin, qmin_coord,    &
               qmin_alo, qmin_ahi, lclose, lrhill_present)
               USE module_parameters
               USE module_swifter
               USE module_rmvs
               IMPLICIT NONE
               LOGICAL(LGT), INTENT(IN)    :: lclose, lrhill_present
               INTEGER(I4B), INTENT(IN)    :: npl
               INTEGER(I4B), INTENT(INOUT) :: ntp, nsp
               REAL(DP), INTENT(IN)        :: t, dt, rmin, rmax, rmaxu, qmin, qmin_alo, qmin_ahi
               CHARACTER(*), INTENT(IN)    :: qmin_coord
               TYPE(rmvs_pl), POINTER      :: rmvs_pl1P
               TYPE(rmvs_tp), POINTER      :: rmvs_tp1P, rmvs_tpd1P
          END SUBROUTINE rmvs_discard
     END INTERFACE

     INTERFACE
          SUBROUTINE rmvs_discard_pl(t, ntp, rmvs_tp1P)
               USE module_parameters
               USE module_rmvs
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN) :: ntp
               REAL(DP), INTENT(IN)     :: t
               TYPE(rmvs_tp), POINTER   :: rmvs_tp1P
          END SUBROUTINE rmvs_discard_pl
     END INTERFACE

     INTERFACE
          SUBROUTINE rmvs_discard_spill(ntp, nsp, rmvs_tp1P, rmvs_tpd1P, rmvs_tpspP)
               USE module_parameters
               USE module_swifter
               USE module_whm
               USE module_rmvs
               IMPLICIT NONE
               INTEGER(I4B), INTENT(INOUT) :: ntp, nsp
               TYPE(rmvs_tp), POINTER      :: rmvs_tp1P, rmvs_tpd1P, rmvs_tpspP
          END SUBROUTINE rmvs_discard_spill
     END INTERFACE

     INTERFACE
          SUBROUTINE rmvs_drift_tp(nenc, rmvs_tpenc1P, mu, dt)
               USE module_parameters
               USE module_swifter
               USE module_rmvs
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN) :: nenc
               REAL(DP), INTENT(IN)     :: mu, dt
               TYPE(rmvs_tp), POINTER   :: rmvs_tpenc1P
          END SUBROUTINE rmvs_drift_tp
     END INTERFACE

     INTERFACE
          SUBROUTINE rmvs_getaccp_ah3_tp(index, npl, nenc, rmvs_pl1P, rmvs_pleP)
               USE module_parameters
               USE module_rmvs
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN) :: index, npl, nenc
               TYPE(rmvs_pl), POINTER   :: rmvs_pl1P, rmvs_pleP
          END SUBROUTINE rmvs_getaccp_ah3_tp
     END INTERFACE

     INTERFACE
          SUBROUTINE rmvs_getaccp_tp(index, lextra_force, t, npl, nplmax, nenc, ntpmax, rmvs_pl1P, rmvs_pleP, j2rp2, j4rp4)
               USE module_parameters
               USE module_rmvs
               IMPLICIT NONE
               LOGICAL(LGT), INTENT(IN) :: lextra_force
               INTEGER(I4B), INTENT(IN) :: index, npl, nplmax, nenc, ntpmax
               REAL(DP), INTENT(IN)     :: t, j2rp2, j4rp4
               TYPE(rmvs_pl), POINTER   :: rmvs_pl1P, rmvs_pleP
          END SUBROUTINE rmvs_getaccp_tp
     END INTERFACE

     INTERFACE
          SUBROUTINE rmvs_interp_in(npl, rmvs_pl1P, dt)
               USE module_parameters
               USE module_rmvs
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN) :: npl
               REAL(DP), INTENT(IN)     :: dt
               TYPE(rmvs_pl), POINTER   :: rmvs_pl1P
          END SUBROUTINE rmvs_interp_in
     END INTERFACE

     INTERFACE
          SUBROUTINE rmvs_interp_out(npl, rmvs_pl1P, dt)
               USE module_parameters
               USE module_rmvs
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN) :: npl
               REAL(DP), INTENT(IN)     :: dt
               TYPE(rmvs_pl), POINTER   :: rmvs_pl1P
          END SUBROUTINE rmvs_interp_out
     END INTERFACE

     INTERFACE
          SUBROUTINE rmvs_kickvp_tp(nenc, rmvs_tpenc1P, dt)
               USE module_parameters
               USE module_swifter
               USE module_rmvs
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN) :: nenc
               REAL(DP), INTENT(IN)     :: dt
               TYPE(rmvs_tp), POINTER   :: rmvs_tpenc1P
          END SUBROUTINE rmvs_kickvp_tp
     END INTERFACE

     INTERFACE
          SUBROUTINE rmvs_peri(lfirst, index, nenc, rmvs_pleP, rmvs_tpenc1P, mu, rhill, t, dt, encounter_file, out_type)
               USE module_parameters
               USE module_rmvs
               IMPLICIT NONE
               LOGICAL(LGT), INTENT(IN) :: lfirst
               INTEGER(I4B), INTENT(IN) :: index, nenc
               REAL(DP), INTENT(IN)     :: mu, rhill, t, dt
               CHARACTER(*), INTENT(IN) :: encounter_file, out_type
               TYPE(rmvs_pl), POINTER   :: rmvs_pleP
               TYPE(rmvs_tp), POINTER   :: rmvs_tpenc1P
          END SUBROUTINE rmvs_peri
     END INTERFACE

     INTERFACE
          SUBROUTINE rmvs_setup(npl, ntp, rmvs_plA, rmvs_tpA, rmvs_pl1P, rmvs_tp1P, swifter_pl1P, swifter_tp1P)
               USE module_parameters
               USE module_swifter
               USE module_whm
               USE module_rmvs
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN)                           :: npl, ntp
               TYPE(swifter_pl), POINTER                          :: swifter_pl1P
               TYPE(swifter_tp), POINTER                          :: swifter_tp1P
               TYPE(rmvs_pl), DIMENSION(:), TARGET, INTENT(INOUT) :: rmvs_plA
               TYPE(rmvs_tp), DIMENSION(:), TARGET, INTENT(INOUT) :: rmvs_tpA
               TYPE(rmvs_pl), POINTER                             :: rmvs_pl1P
               TYPE(rmvs_tp), POINTER                             :: rmvs_tp1P
          END SUBROUTINE rmvs_setup
     END INTERFACE

     INTERFACE
          SUBROUTINE rmvs_step(lfirst, lextra_force, t, npl, nplmax, ntp, ntpmax, rmvs_pl1P, rmvs_tp1P, j2rp2, j4rp4, dt,         &
               encounter_file, out_type)
               USE module_parameters
               USE module_swifter
               USE module_whm
               USE module_rmvs
               IMPLICIT NONE
               LOGICAL(LGT), INTENT(IN)    :: lextra_force
               LOGICAL(LGT), INTENT(INOUT) :: lfirst
               INTEGER(I4B), INTENT(IN)    :: npl, nplmax, ntp, ntpmax
               REAL(DP), INTENT(IN)        :: t, j2rp2, j4rp4, dt
               CHARACTER(*), INTENT(IN)    :: encounter_file, out_type
               TYPE(rmvs_pl), POINTER      :: rmvs_pl1P
               TYPE(rmvs_tp), POINTER      :: rmvs_tp1P
          END SUBROUTINE rmvs_step
     END INTERFACE

     INTERFACE
          SUBROUTINE rmvs_step_in(lextra_force, t, npl, nplmax, ntp, ntpmax, rmvs_pl1P, rmvs_tp1P, j2rp2, j4rp4, dt,              &
               encounter_file, out_type)
               USE module_parameters
               USE module_swifter
               USE module_rmvs
               IMPLICIT NONE
               LOGICAL(LGT), INTENT(IN) :: lextra_force
               INTEGER(I4B), INTENT(IN) :: npl, nplmax, ntp, ntpmax
               REAL(DP), INTENT(IN)     :: t, j2rp2, j4rp4, dt
               CHARACTER(*), INTENT(IN) :: encounter_file, out_type
               TYPE(rmvs_pl), POINTER   :: rmvs_pl1P
               TYPE(rmvs_tp), POINTER   :: rmvs_tp1P
          END SUBROUTINE rmvs_step_in
     END INTERFACE

     INTERFACE
          SUBROUTINE rmvs_step_in_tp(index, lfirst, lextra_force, t, npl, nplmax, nenc, ntpmax, rmvs_pl1P, rmvs_pleP, j2rp2,      &
               j4rp4, dt)
               USE module_parameters
               USE module_rmvs
               IMPLICIT NONE
               LOGICAL(LGT), INTENT(IN)    :: lextra_force
               LOGICAL(LGT), INTENT(INOUT) :: lfirst
               INTEGER(I4B), INTENT(IN)    :: index, npl, nplmax, nenc, ntpmax
               REAL(DP), INTENT(IN)        :: t, j2rp2, j4rp4, dt
               TYPE(rmvs_pl), POINTER      :: rmvs_pl1P, rmvs_pleP
          END SUBROUTINE rmvs_step_in_tp
     END INTERFACE

     INTERFACE
          SUBROUTINE rmvs_step_out2(index, lfirst, lextra_force, t, npl, nplmax, ntp, ntpmax, rmvs_pl1P, rmvs_tp1P, j2rp2, j4rp4, &
               dt, encounter_file, out_type)
               USE module_parameters
               USE module_whm
               USE module_rmvs
               IMPLICIT NONE
               LOGICAL(LGT), INTENT(IN) :: lfirst, lextra_force
               INTEGER(I4B), INTENT(IN) :: index, npl, nplmax, ntp, ntpmax
               REAL(DP), INTENT(IN)     :: t, j2rp2, j4rp4, dt
               CHARACTER(*), INTENT(IN) :: encounter_file, out_type
               TYPE(rmvs_pl), POINTER   :: rmvs_pl1P
               TYPE(rmvs_tp), POINTER   :: rmvs_tp1P
          END SUBROUTINE rmvs_step_out2
     END INTERFACE

     INTERFACE
          SUBROUTINE rmvs_step_out(lfirst, lextra_force, t, npl, nplmax, ntp, ntpmax, rmvs_pl1P, rmvs_tp1P, j2rp2, j4rp4, dt,     &
               encounter_file, out_type)
               USE module_parameters
               USE module_rmvs
               IMPLICIT NONE
               LOGICAL(LGT), INTENT(IN)    :: lfirst, lextra_force
               INTEGER(I4B), INTENT(IN)    :: npl, nplmax, ntp, ntpmax
               REAL(DP), INTENT(IN)        :: t, j2rp2, j4rp4, dt
               CHARACTER(*), INTENT(IN)    :: encounter_file, out_type
               TYPE(rmvs_pl), POINTER      :: rmvs_pl1P
               TYPE(rmvs_tp), POINTER      :: rmvs_tp1P
          END SUBROUTINE rmvs_step_out
     END INTERFACE

     INTERFACE
          SUBROUTINE rmvs_user_getacch_tp(t, nenc, rmvs_tpenc1P)
               USE module_parameters
               USE module_rmvs
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN) :: nenc
               REAL(DP), INTENT(IN)     :: t
               TYPE(rmvs_tp), POINTER   :: rmvs_tpenc1P
          END SUBROUTINE rmvs_user_getacch_tp
     END INTERFACE

     INTERFACE
          SUBROUTINE symba_chk(xr, vr, rhill1, rhill2, dt, irec, lencounter, lvdotr)
               USE module_parameters
               USE module_swifter
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
          SUBROUTINE symba_discard_merge_pl(t, npl, nsppl, symba_pl1P, symba_pld1P, nplplenc, plplenc_list)
               USE module_parameters
               USE module_swifter
               USE module_helio
               USE module_symba
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN)                      :: nplplenc
               INTEGER(I4B), INTENT(INOUT)                   :: npl, nsppl
               REAL(DP), INTENT(IN)                          :: t
               TYPE(symba_pl), POINTER                       :: symba_pl1P, symba_pld1P
               TYPE(symba_plplenc), DIMENSION(:), INTENT(IN) :: plplenc_list
          END SUBROUTINE symba_discard_merge_pl
     END INTERFACE

     INTERFACE
          SUBROUTINE symba_discard_peri_pl(t, npl, symba_pl1P, msys, qmin, qmin_alo, qmin_ahi, qmin_coord, ldiscards)
               USE module_parameters
               USE module_swifter
               USE module_helio
               USE module_symba
               IMPLICIT NONE
               LOGICAL(LGT), INTENT(INOUT) :: ldiscards
               INTEGER(I4B), INTENT(IN)    :: npl
               REAL(DP), INTENT(IN)        :: t, msys, qmin, qmin_alo, qmin_ahi
               CHARACTER(*), INTENT(IN)    :: qmin_coord
               TYPE(symba_pl), POINTER     :: symba_pl1P
          END SUBROUTINE symba_discard_peri_pl
     END INTERFACE

     INTERFACE
          SUBROUTINE symba_discard_pl(t, npl, nplmax, nsp, symba_pl1P, symba_pld1P, rmin, rmax, rmaxu, qmin, qmin_coord,          &
               qmin_alo, qmin_ahi, j2rp2, j4rp4, eoffset)
               USE module_parameters
               USE module_swifter
               USE module_helio
               USE module_symba
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN)    :: nplmax
               INTEGER(I4B), INTENT(INOUT) :: npl, nsp
               REAL(DP), INTENT(IN)        :: t, rmin, rmax, rmaxu, qmin, qmin_alo, qmin_ahi, j2rp2, j4rp4
               REAL(DP), INTENT(INOUT)     :: eoffset
               CHARACTER(*), INTENT(IN)    :: qmin_coord
               TYPE(symba_pl), POINTER     :: symba_pl1P, symba_pld1P
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
          SUBROUTINE symba_discard_sun_pl(t, npl, msys, swifter_pl1P, rmin, rmax, rmaxu, ldiscards)
               USE module_parameters
               USE module_swifter
               IMPLICIT NONE
               LOGICAL(LGT), INTENT(INOUT) :: ldiscards
               INTEGER(I4B), INTENT(IN)    :: npl
               REAL(DP), INTENT(IN)        :: t, msys, rmin, rmax, rmaxu
               TYPE(swifter_pl), POINTER   :: swifter_pl1P
          END SUBROUTINE symba_discard_sun_pl
     END INTERFACE

     INTERFACE
          SUBROUTINE symba_discard_tp(t, npl, ntp, nsp, symba_pl1P, symba_tp1P, symba_tpd1P, dt, rmin, rmax, rmaxu, qmin,         &
               qmin_coord, qmin_alo, qmin_ahi, lclose, lrhill_present)
               USE module_parameters
               USE module_swifter
               USE module_helio
               USE module_symba
               IMPLICIT NONE
               LOGICAL(LGT), INTENT(IN)    :: lclose, lrhill_present
               INTEGER(I4B), INTENT(IN)    :: npl
               INTEGER(I4B), INTENT(INOUT) :: ntp, nsp
               REAL(DP), INTENT(IN)        :: t, dt, rmin, rmax, rmaxu, qmin, qmin_alo, qmin_ahi
               CHARACTER(*), INTENT(IN)    :: qmin_coord
               TYPE(symba_pl), POINTER     :: symba_pl1P
               TYPE(symba_tp), POINTER     :: symba_tp1P,symba_tpd1P
          END SUBROUTINE symba_discard_tp
     END INTERFACE

     INTERFACE
          SUBROUTINE symba_energy(npl, nplmax, swifter_pl1P, j2rp2, j4rp4, ke, pe, te, htot)
               USE module_parameters
               USE module_swifter
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN)               :: npl, nplmax
               REAL(DP), INTENT(IN)                   :: j2rp2, j4rp4
               REAL(DP), INTENT(OUT)                  :: ke, pe, te
               REAL(DP), DIMENSION(NDIM), INTENT(OUT) :: htot
               TYPE(swifter_pl), POINTER              :: swifter_pl1P
          END SUBROUTINE symba_energy
     END INTERFACE

     INTERFACE
          SUBROUTINE symba_getacch(lextra_force, t, npl, nplm, nplmax, symba_pl1P, j2rp2, j4rp4, nplplenc, plplenc_list)
               USE module_parameters
               USE module_swifter
               USE module_helio
               USE module_symba
               IMPLICIT NONE
               LOGICAL(LGT), INTENT(IN)                      :: lextra_force
               INTEGER(I4B), INTENT(IN)                      :: npl, nplm, nplmax, nplplenc
               REAL(DP), INTENT(IN)                          :: t, j2rp2, j4rp4
               TYPE(symba_pl), POINTER                       :: symba_pl1P
               TYPE(symba_plplenc), DIMENSION(:), INTENT(IN) :: plplenc_list
          END SUBROUTINE symba_getacch
     END INTERFACE

     INTERFACE
          SUBROUTINE symba_getacch_tp(lextra_force, t, npl, nplm, nplmax, ntp, ntpmax, symba_pl1P, symba_tp1P, xh, j2rp2, j4rp4,  &
               npltpenc, pltpenc_list)
               USE module_parameters
               USE module_swifter
               USE module_helio
               USE module_symba
               IMPLICIT NONE
               LOGICAL(LGT), INTENT(IN)                      :: lextra_force
               INTEGER(I4B), INTENT(IN)                      :: npl, nplm, nplmax, ntp, ntpmax, npltpenc
               REAL(DP), INTENT(IN)                          :: t, j2rp2, j4rp4
               REAL(DP), DIMENSION(NDIM, npl), INTENT(IN)    :: xh
               TYPE(symba_pl), POINTER                       :: symba_pl1P
               TYPE(symba_tp), POINTER                       :: symba_tp1P
               TYPE(symba_pltpenc), DIMENSION(:), INTENT(IN) :: pltpenc_list
          END SUBROUTINE symba_getacch_tp
     END INTERFACE

     INTERFACE
          SUBROUTINE symba_helio_drift(irec, npl, symba_pl1P, dt)
               USE module_parameters
               USE module_swifter
               USE module_helio
               USE module_symba
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN) :: irec, npl
               REAL(DP), INTENT(IN)     :: dt
               TYPE(symba_pl), POINTER  :: symba_pl1P
          END SUBROUTINE symba_helio_drift
     END INTERFACE

     INTERFACE
          SUBROUTINE symba_helio_drift_tp(irec, ntp, symba_tp1P, mu, dt)
               USE module_parameters
               USE module_swifter
               USE module_helio
               USE module_symba
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN) :: irec, ntp
               REAL(DP), INTENT(IN)     :: mu, dt
               TYPE(symba_tp), POINTER  :: symba_tp1P
          END SUBROUTINE symba_helio_drift_tp
     END INTERFACE

     INTERFACE
          SUBROUTINE symba_helio_getacch(lflag, lextra_force, t, npl, nplm, nplmax, helio_pl1P, j2rp2, j4rp4)
               USE module_parameters
               USE module_swifter
               USE module_helio
               USE module_symba
               IMPLICIT NONE
               LOGICAL(LGT), INTENT(IN) :: lflag, lextra_force
               INTEGER(I4B), INTENT(IN) :: npl, nplm, nplmax
               REAL(DP), INTENT(IN)     :: t, j2rp2, j4rp4
               TYPE(helio_pl), POINTER  :: helio_pl1P
          END SUBROUTINE symba_helio_getacch
     END INTERFACE

     INTERFACE
          SUBROUTINE symba_helio_getacch_int(npl, nplm, helio_pl1P)
               USE module_parameters
               USE module_helio
               USE module_symba
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN) :: npl, nplm
               TYPE(helio_pl), POINTER  :: helio_pl1P
          END SUBROUTINE symba_helio_getacch_int
     END INTERFACE

     INTERFACE
          SUBROUTINE symba_kick(irec, nplplenc, npltpenc, plplenc_list, pltpenc_list, dt, sgn)
               USE module_parameters
               USE module_swifter
               USE module_helio
               USE module_symba
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN)                      :: irec, nplplenc, npltpenc
               REAL(DP), INTENT(IN)                          :: dt, sgn
               TYPE(symba_plplenc), DIMENSION(:), INTENT(IN) :: plplenc_list
               TYPE(symba_pltpenc), DIMENSION(:), INTENT(IN) :: pltpenc_list
          END SUBROUTINE symba_kick
     END INTERFACE

     INTERFACE
          SUBROUTINE symba_merge_pl(t, dt, index, nplplenc, plplenc_list, nmergeadd, nmergesub, mergeadd_list, mergesub_list,     &
               eoffset, vbs, encounter_file, out_type)
               USE module_parameters
               USE module_swifter
               USE module_helio
               USE module_symba
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN)                         :: index, nplplenc
               INTEGER(I4B), INTENT(INOUT)                      :: nmergeadd, nmergesub
               REAL(DP), INTENT(IN)                             :: t, dt
               REAL(DP), INTENT(INOUT)                          :: eoffset
               REAL(DP), DIMENSION(NDIM), INTENT(IN)            :: vbs
               CHARACTER(*), INTENT(IN)                         :: encounter_file, out_type
               TYPE(symba_plplenc), DIMENSION(:), INTENT(INOUT) :: plplenc_list
               TYPE(symba_merger), DIMENSION(:), INTENT(INOUT)  :: mergeadd_list, mergesub_list
          END SUBROUTINE symba_merge_pl
     END INTERFACE

     INTERFACE
          SUBROUTINE symba_merge_tp(t, dt, index, npltpenc, pltpenc_list, vbs, encounter_file, out_type)
               USE module_parameters
               USE module_swifter
               USE module_helio
               USE module_symba
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN)                         :: index, npltpenc
               REAL(DP), INTENT(IN)                             :: t, dt
               REAL(DP), DIMENSION(NDIM), INTENT(IN)            :: vbs
               CHARACTER(*), INTENT(IN)                         :: encounter_file, out_type
               TYPE(symba_pltpenc), DIMENSION(:), INTENT(INOUT) :: pltpenc_list
          END SUBROUTINE symba_merge_tp
     END INTERFACE

     INTERFACE
          SUBROUTINE symba_peri(lfirst, npl, symba_pl1P, msys, qmin_coord)
               USE module_parameters
               USE module_swifter
               USE module_helio
               USE module_symba
               IMPLICIT NONE
               LOGICAL(LGT), INTENT(IN) :: lfirst
               INTEGER(I4B), INTENT(IN) :: npl
               REAL(DP), INTENT(IN)     :: msys
               CHARACTER(*), INTENT(IN) :: qmin_coord
               TYPE(symba_pl), POINTER  :: symba_pl1P
          END SUBROUTINE symba_peri
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
          SUBROUTINE symba_setup(npl, ntp, symba_plA, symba_tpA, symba_pl1P, symba_tp1P, swifter_pl1P, swifter_tp1P)
               USE module_parameters
               USE module_swifter
               USE module_helio
               USE module_symba
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN)                            :: npl, ntp
               TYPE(swifter_pl), POINTER                           :: swifter_pl1P
               TYPE(swifter_tp), POINTER                           :: swifter_tp1P
               TYPE(symba_pl), DIMENSION(:), TARGET, INTENT(INOUT) :: symba_plA
               TYPE(symba_tp), DIMENSION(:), TARGET, INTENT(INOUT) :: symba_tpA
               TYPE(symba_pl), POINTER                             :: symba_pl1P
               TYPE(symba_tp), POINTER                             :: symba_tp1P
          END SUBROUTINE symba_setup
     END INTERFACE

     INTERFACE
          SUBROUTINE symba_step(lfirst, lextra_force, lclose, t, npl, nplmax, ntp, ntpmax, symba_pl1P, symba_tp1P, j2rp2, j4rp4,  &
               dt, nplplenc, npltpenc, plplenc_list, pltpenc_list, nmergeadd, nmergesub, mergeadd_list, mergesub_list, eoffset,   &
               mtiny, encounter_file, out_type)
               USE module_parameters
               USE module_swifter
               USE module_helio
               USE module_symba
               IMPLICIT NONE
               LOGICAL(LGT), INTENT(IN)                         :: lextra_force, lclose
               LOGICAL(LGT), INTENT(INOUT)                      :: lfirst
               INTEGER(I4B), INTENT(IN)                         :: npl, nplmax, ntp, ntpmax
               INTEGER(I4B), INTENT(INOUT)                      :: nplplenc, npltpenc, nmergeadd, nmergesub
               REAL(DP), INTENT(IN)                             :: t, j2rp2, j4rp4, dt, mtiny
               REAL(DP), INTENT(INOUT)                          :: eoffset
               CHARACTER(*), INTENT(IN)                         :: encounter_file, out_type
               TYPE(symba_pl), POINTER                          :: symba_pl1P
               TYPE(symba_tp), POINTER                          :: symba_tp1P
               TYPE(symba_plplenc), DIMENSION(:), INTENT(INOUT) :: plplenc_list
               TYPE(symba_pltpenc), DIMENSION(:), INTENT(INOUT) :: pltpenc_list
               TYPE(symba_merger), DIMENSION(:), INTENT(INOUT)  :: mergeadd_list, mergesub_list
          END SUBROUTINE symba_step
     END INTERFACE

     INTERFACE
          SUBROUTINE symba_step_helio(lfirst, lextra_force, t, npl, nplm, nplmax, ntp, ntpmax, helio_pl1P, helio_tp1P, j2rp2,     &
               j4rp4, dt)
               USE module_parameters
               USE module_swifter
               USE module_helio
               USE module_symba
               IMPLICIT NONE
               LOGICAL(LGT), INTENT(IN)    :: lextra_force
               LOGICAL(LGT), INTENT(INOUT) :: lfirst
               INTEGER(I4B), INTENT(IN)    :: npl, nplm, nplmax, ntp, ntpmax
               REAL(DP), INTENT(IN)        :: t, j2rp2, j4rp4, dt
               TYPE(helio_pl), POINTER     :: helio_pl1P
               TYPE(helio_tp), POINTER     :: helio_tp1P
          END SUBROUTINE symba_step_helio
     END INTERFACE

     INTERFACE
          SUBROUTINE symba_step_helio_pl(lfirst, lextra_force, t, npl, nplm, nplmax, helio_pl1P, j2rp2, j4rp4, dt, xbeg, xend,    &
               ptb, pte)
               USE module_parameters
               USE module_swifter
               USE module_helio
               USE module_symba
               IMPLICIT NONE
               LOGICAL(LGT), INTENT(IN)                     :: lextra_force
               LOGICAL(LGT), INTENT(INOUT)                  :: lfirst
               INTEGER(I4B), INTENT(IN)                     :: npl, nplm, nplmax
               REAL(DP), INTENT(IN)                         :: t, j2rp2, j4rp4, dt
               REAL(DP), DIMENSION(NDIM, nplm), INTENT(OUT) :: xbeg, xend
               REAL(DP), DIMENSION(NDIM), INTENT(OUT)       :: ptb, pte
               TYPE(helio_pl), POINTER                      :: helio_pl1P
          END SUBROUTINE symba_step_helio_pl
     END INTERFACE

     INTERFACE
          SUBROUTINE symba_step_interp(lextra_force, lclose, t, npl, nplm, nplmax, ntp, ntpmax, symba_pl1P, symba_tp1P, j2rp2,    &
               j4rp4, dt, eoffset, mtiny, nplplenc, npltpenc, plplenc_list, pltpenc_list, nmergeadd, nmergesub, mergeadd_list,    &
               mergesub_list, encounter_file, out_type)
               USE module_parameters
               USE module_swifter
               USE module_helio
               USE module_symba
               IMPLICIT NONE
               LOGICAL(LGT), INTENT(IN)                         :: lextra_force, lclose
               INTEGER(I4B), INTENT(IN)                         :: npl, nplm, nplmax, ntp, ntpmax, nplplenc, npltpenc
               INTEGER(I4B), INTENT(INOUT)                      :: nmergeadd, nmergesub
               REAL(DP), INTENT(IN)                             :: t, j2rp2, j4rp4, dt, mtiny
               REAL(DP), INTENT(INOUT)                          :: eoffset
               CHARACTER(*), INTENT(IN)                         :: encounter_file, out_type
               TYPE(symba_pl), POINTER                          :: symba_pl1P
               TYPE(symba_tp), POINTER                          :: symba_tp1P
               TYPE(symba_plplenc), DIMENSION(:), INTENT(INOUT) :: plplenc_list
               TYPE(symba_pltpenc), DIMENSION(:), INTENT(INOUT) :: pltpenc_list
               TYPE(symba_merger), DIMENSION(:), INTENT(INOUT)  :: mergeadd_list, mergesub_list
          END SUBROUTINE symba_step_interp
     END INTERFACE

     INTERFACE
          RECURSIVE SUBROUTINE symba_step_recur(lclose, t, ireci, npl, nplm, ntp, symba_pl1P, symba_tp1P, dt0, eoffset, nplplenc, &
               npltpenc, plplenc_list, pltpenc_list, nmergeadd, nmergesub, mergeadd_list, mergesub_list, encounter_file, out_type)
               USE module_parameters
               USE module_swifter
               USE module_helio
               USE module_symba
               IMPLICIT NONE
               LOGICAL(LGT), INTENT(IN)                         :: lclose
               INTEGER(I4B), INTENT(IN)                         :: ireci, npl, nplm, ntp, nplplenc, npltpenc
               INTEGER(I4B), INTENT(INOUT)                      :: nmergeadd, nmergesub
               REAL(DP), INTENT(IN)                             :: t, dt0
               REAL(DP), INTENT(INOUT)                          :: eoffset
               CHARACTER(*), INTENT(IN)                         :: encounter_file, out_type
               TYPE(symba_pl), POINTER                          :: symba_pl1P
               TYPE(symba_tp), POINTER                          :: symba_tp1P
               TYPE(symba_plplenc), DIMENSION(:), INTENT(INOUT) :: plplenc_list
               TYPE(symba_pltpenc), DIMENSION(:), INTENT(INOUT) :: pltpenc_list
               TYPE(symba_merger), DIMENSION(:), INTENT(INOUT)  :: mergeadd_list, mergesub_list
          END SUBROUTINE symba_step_recur
     END INTERFACE

     INTERFACE
          SUBROUTINE symba_user_getacch(t, npl, symba_pl1P)
               USE module_parameters
               USE module_symba
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN) :: npl
               REAL(DP), INTENT(IN)     :: t
               TYPE(symba_pl), POINTER  :: symba_pl1P
          END SUBROUTINE symba_user_getacch
     END INTERFACE

     INTERFACE
          SUBROUTINE symba_user_getacch_tp(t, ntp, symba_tp1P)
               USE module_parameters
               USE module_symba
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN) :: ntp
               REAL(DP), INTENT(IN)     :: t
               TYPE(symba_tp), POINTER  :: symba_tp1P
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
          SUBROUTINE util_hills(npl, swifter_pl1P)
               USE module_parameters
               USE module_swifter
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN)  :: npl
               TYPE(swifter_pl), POINTER :: swifter_pl1P
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
          SUBROUTINE util_peri(lfirst, ntp, swifter_tp1P, mu, msys, qmin_coord)
               USE module_parameters
               USE module_swifter
               IMPLICIT NONE
               LOGICAL(LGT), INTENT(IN)  :: lfirst
               INTEGER(I4B), INTENT(IN)  :: ntp
               REAL(DP), INTENT(IN)      :: mu, msys
               CHARACTER(*), INTENT(IN)  :: qmin_coord
               TYPE(swifter_tp), POINTER :: swifter_tp1P
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
          SUBROUTINE util_valid(npl, ntp, swifter_pl1P, swifter_tp1P)
               USE module_parameters
               USE module_swifter
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN)  :: npl, ntp
               TYPE(swifter_pl), POINTER :: swifter_pl1P
               TYPE(swifter_tp), POINTER :: swifter_tp1P
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
            REAL(DP) :: util_kahan_sum
            REAL(DP), INTENT(IN) :: xsum_current, xi
            REAL(DP), INTENT(INOUT) :: xerror
         END FUNCTION
     END INTERFACE

     INTERFACE
          SUBROUTINE whm_discard(t, npl, ntp, nsp, whm_pl1P, whm_tp1P, whm_tpd1P, dt, rmin, rmax, rmaxu, qmin, qmin_coord,        &
               qmin_alo, qmin_ahi, lclose, lrhill_present)
               USE module_parameters
               USE module_swifter
               USE module_whm
               IMPLICIT NONE
               LOGICAL(LGT), INTENT(IN)    :: lclose, lrhill_present
               INTEGER(I4B), INTENT(IN)    :: npl
               INTEGER(I4B), INTENT(INOUT) :: ntp, nsp
               REAL(DP), INTENT(IN)        :: t, dt, rmin, rmax, rmaxu, qmin, qmin_alo, qmin_ahi
               CHARACTER(*), INTENT(IN)    :: qmin_coord
               TYPE(whm_pl), POINTER       :: whm_pl1P
               TYPE(whm_tp), POINTER       :: whm_tp1P, whm_tpd1P
          END SUBROUTINE whm_discard
     END INTERFACE

     INTERFACE
          SUBROUTINE whm_discard_spill(ntp, nsp, whm_tp1P, whm_tpd1P, whm_tpspP)
               USE module_parameters
               USE module_swifter
               USE module_whm
               IMPLICIT NONE
               INTEGER(I4B), INTENT(INOUT) :: ntp, nsp
               TYPE(whm_tp), POINTER       :: whm_tp1P, whm_tpd1P, whm_tpspP
          END SUBROUTINE whm_discard_spill
     END INTERFACE

     INTERFACE
          SUBROUTINE whm_drift(npl, whm_pl1P, dt)
               USE module_parameters
               USE module_swifter
               USE module_whm
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN) :: npl
               REAL(DP), INTENT(IN)     :: dt
               TYPE(whm_pl), POINTER    :: whm_pl1P
          END SUBROUTINE whm_drift
     END INTERFACE

     INTERFACE
          SUBROUTINE whm_drift_tp(ntp, whm_tp1P, mu, dt)
               USE module_parameters
               USE module_swifter
               USE module_whm
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN) :: ntp
               REAL(DP), INTENT(IN)     :: mu, dt
               TYPE(whm_tp), POINTER    :: whm_tp1P
          END SUBROUTINE whm_drift_tp
     END INTERFACE

     INTERFACE
          SUBROUTINE whm_getacch_ah1(npl, whm_pl1P, ir3h, ir3j)
               USE module_parameters
               USE module_whm
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN)           :: npl
               REAL(DP), DIMENSION(:), INTENT(IN) :: ir3h, ir3j
               TYPE(whm_pl), POINTER              :: whm_pl1P
          END SUBROUTINE whm_getacch_ah1
     END INTERFACE

     INTERFACE
          SUBROUTINE whm_getacch_ah2(npl, whm_pl1P, ir3j)
               USE module_parameters
               USE module_whm
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN)           :: npl
               REAL(DP), DIMENSION(:), INTENT(IN) :: ir3j
               TYPE(whm_pl), POINTER              :: whm_pl1P
          END SUBROUTINE whm_getacch_ah2
     END INTERFACE

     INTERFACE
          SUBROUTINE whm_getacch_ah3(npl, whm_pl1P)
               USE module_parameters
               USE module_whm
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN) :: npl
               TYPE(whm_pl), POINTER    :: whm_pl1P
          END SUBROUTINE whm_getacch_ah3
     END INTERFACE

     INTERFACE
          SUBROUTINE whm_getacch_ah3_tp(npl, ntp, whm_pl1P, whm_tp1P, xh)
               USE module_parameters
               USE module_whm
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN)                   :: npl, ntp
               REAL(DP), DIMENSION(NDIM, npl), INTENT(IN) :: xh
               TYPE(whm_pl), POINTER                      :: whm_pl1P
               TYPE(whm_tp), POINTER                      :: whm_tp1P
          END SUBROUTINE whm_getacch_ah3_tp
     END INTERFACE

     INTERFACE
          SUBROUTINE whm_getacch(lextra_force, t, npl, nplmax, whm_pl1P, j2rp2, j4rp4)
               USE module_parameters
               USE module_swifter
               USE module_whm
               IMPLICIT NONE
               LOGICAL(LGT), INTENT(IN) :: lextra_force
               INTEGER(I4B), INTENT(IN) :: npl, nplmax
               REAL(DP), INTENT(IN)     :: t, j2rp2, j4rp4
               TYPE(whm_pl), POINTER    :: whm_pl1P
          END SUBROUTINE whm_getacch
     END INTERFACE

     INTERFACE
          SUBROUTINE whm_getacch_tp(lextra_force, t, npl, nplmax, ntp, ntpmax, whm_pl1P, whm_tp1P, xh, j2rp2, j4rp4)
               USE module_parameters
               USE module_swifter
               USE module_whm
               IMPLICIT NONE
               LOGICAL(LGT), INTENT(IN)                   :: lextra_force
               INTEGER(I4B), INTENT(IN)                   :: npl, nplmax, ntp, ntpmax
               REAL(DP), INTENT(IN)                       :: t, j2rp2, j4rp4
               REAL(DP), DIMENSION(NDIM, npl), INTENT(IN) :: xh
               TYPE(whm_pl), POINTER                      :: whm_pl1P
               TYPE(whm_tp), POINTER                      :: whm_tp1P
          END SUBROUTINE whm_getacch_tp
     END INTERFACE

     INTERFACE
          SUBROUTINE whm_kickvh(npl, whm_pl1P, dt)
               USE module_parameters
               USE module_swifter
               USE module_whm
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN) :: npl
               REAL(DP), INTENT(IN)     :: dt
               TYPE(whm_pl), POINTER    :: whm_pl1P
          END SUBROUTINE whm_kickvh
     END INTERFACE

     INTERFACE
          SUBROUTINE whm_kickvh_tp(ntp, whm_tp1P, dt)
               USE module_parameters
               USE module_swifter
               USE module_whm
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN) :: ntp
               REAL(DP), INTENT(IN)     :: dt
               TYPE(whm_tp), POINTER    :: whm_tp1P
          END SUBROUTINE whm_kickvh_tp
     END INTERFACE

     INTERFACE
          SUBROUTINE whm_setup(npl, ntp, whm_plA, whm_tpA, whm_pl1P, whm_tp1P, swifter_pl1P, swifter_tp1P)
               USE module_parameters
               USE module_swifter
               USE module_whm
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN)                          :: npl, ntp
               TYPE(swifter_pl), POINTER                         :: swifter_pl1P
               TYPE(swifter_tp), POINTER                         :: swifter_tp1P
               TYPE(whm_pl), DIMENSION(:), TARGET, INTENT(INOUT) :: whm_plA
               TYPE(whm_tp), DIMENSION(:), TARGET, INTENT(INOUT) :: whm_tpA
               TYPE(whm_pl), POINTER                             :: whm_pl1P
               TYPE(whm_tp), POINTER                             :: whm_tp1P
          END SUBROUTINE whm_setup
     END INTERFACE

     INTERFACE
          SUBROUTINE whm_step(lfirst, lextra_force, t, npl, nplmax, ntp, ntpmax, whm_pl1P, whm_tp1P, j2rp2, j4rp4, dt)
               USE module_parameters
               USE module_whm
               IMPLICIT NONE
               LOGICAL(LGT), INTENT(IN)    :: lextra_force
               LOGICAL(LGT), INTENT(INOUT) :: lfirst
               INTEGER(I4B), INTENT(IN)    :: npl, nplmax, ntp, ntpmax
               REAL(DP), INTENT(IN)        :: t, j2rp2, j4rp4, dt
               TYPE(whm_pl), POINTER       :: whm_pl1P
               TYPE(whm_tp), POINTER       :: whm_tp1P
          END SUBROUTINE whm_step
     END INTERFACE

     INTERFACE
          SUBROUTINE whm_step_pl(lfirst, lextra_force, t, npl, nplmax, whm_pl1P, j2rp2, j4rp4, dt)
               USE module_parameters
               USE module_whm
               IMPLICIT NONE
               LOGICAL(LGT), INTENT(IN)    :: lextra_force
               LOGICAL(LGT), INTENT(INOUT) :: lfirst
               INTEGER(I4B), INTENT(IN)    :: npl, nplmax
               REAL(DP), INTENT(IN)        :: t, j2rp2, j4rp4, dt
               TYPE(whm_pl), POINTER       :: whm_pl1P
          END SUBROUTINE whm_step_pl
     END INTERFACE

     INTERFACE
          SUBROUTINE whm_step_tp(lfirsttp, lextra_force, t, npl, nplmax, ntp, ntpmax, whm_pl1P, whm_tp1P, xbeg, xend, j2rp2,      &
               j4rp4, dt)
               USE module_parameters
               USE module_whm
               IMPLICIT NONE
               LOGICAL(LGT), INTENT(IN)                   :: lfirsttp, lextra_force
               INTEGER(I4B), INTENT(IN)                   :: npl, nplmax, ntp, ntpmax
               REAL(DP), INTENT(IN)                       :: t, j2rp2, j4rp4, dt
               REAL(DP), DIMENSION(NDIM, npl), INTENT(IN) :: xbeg, xend
               TYPE(whm_pl), POINTER                      :: whm_pl1P
               TYPE(whm_tp), POINTER                      :: whm_tp1P
          END SUBROUTINE whm_step_tp
     END INTERFACE

     INTERFACE
          SUBROUTINE whm_user_getacch(t, npl, whm_pl1P)
               USE module_parameters
               USE module_whm
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN) :: npl
               REAL(DP), INTENT(IN)     :: t
               TYPE(whm_pl), POINTER    :: whm_pl1P
          END SUBROUTINE whm_user_getacch
     END INTERFACE

     INTERFACE
          SUBROUTINE whm_user_getacch_tp(t, ntp, whm_tp1P)
               USE module_parameters
               USE module_whm
               IMPLICIT NONE
               INTEGER(I4B), INTENT(IN) :: ntp
               REAL(DP), INTENT(IN)     :: t
               TYPE(whm_tp), POINTER    :: whm_tp1P
          END SUBROUTINE whm_user_getacch_tp
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
