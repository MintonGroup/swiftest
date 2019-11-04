!**********************************************************************************************************************************
!
!  Unit Name   : symba_discard_all
!  Unit Type   : subroutine
!  Project     : Swiftest
!  Package     : symba
!  Language    : Fortran 90/95
!
!  Description : Check to see if planets should be discarded based on their positions or because they are unbound
!
!  Input
!    Arguments : t           : time
!                npl         : number of planets
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : npl         : number of planets
!    Terminal  : none
!    File      : none
!
!  Invocation  : CALL symba_discard_pl(t, npl, nplmax, nsp, symba_pl1P, symba_pld1P, rmin, rmax, rmaxu, qmin, qmin_coord,
!                                      qmin_alo, qmin_ahi, j2rp2, j4rp4, eoffset)
!
!  Notes       : Adapted from Hal Levison's Swift routine discard_massive5.f
!
!**********************************************************************************************************************************
SUBROUTINE symba_discard_all(t, npl, ntp, symba_plA, symba_tpA, symba_pldA, symba_tpdA)

DO i = 2, npl
               !symba_plspP => symba_plP
               !symba_plP => symba_plP%nextP
               !swifter_plspP => symba_plspP%helio%swifter
               IF (swifter_plspP%status /= ACTIVE) CALL symba_discard_spill_pl(npl, nsp, symba_pld1P, symba_plspP)
          END DO