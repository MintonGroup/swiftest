!**********************************************************************************************************************************
!
!  Unit Name   : symba_rearray
!  Unit Type   : subroutine
!  Project     : Swiftest
!  Package     : symba
!  Language    : Fortran 90/95
!
!  Description : Redo array of pl and tp based on discarded and added pl and tp
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
SUBROUTINE symba_discard_all(t, npl, ntp, nsppl, nsptp, symba_plA, symba_tpA, nmergeadd, mergeadd_list, nmergesub, mergesub_list, symba_pldA, symba_tpdA)

! Modules
     USE module_parameters
     USE module_swiftest
     USE module_helio
     USE module_symba
     USE module_interfaces, EXCEPT_THIS_ONE => symba_discard_pl
     IMPLICIT NONE

! Arguments
     INTEGER(I4B), INTENT(IN)                     :: npl, ntp, nsppl, nsptp, nmergeadd, nmergesub
     REAL(DP), INTENT(IN)                         :: t
     TYPE(symba_pl), INTENT(INOUT)                :: symba_plA, symba_pldA
     TYPE(symba_tp), INTENT(INOUT)                :: symba_tpdA, symba_tpdA
     TYPE(symba_merger), DIMENSION(:), INTENT(IN) :: mergeadd_list, mergesub_list

! Internals
     INTEGER(I4B)              :: i, index, j, ncomp, ierr, nplm

    CALL symba_energy(npl, nplmax, symba_plA%helio%swiftest, j2rp2, j4rp4, ke, pe, tei, htot)

    symba_pldA = PACK(symba_plA,symba_plA%helio%swiftest%status /= ACTIVE)
    symba_plA = PACK(symba_plA,symba_plA%helio%swiftest%status = ACTIVE)
    nsppl = SIZE(symba_pldA, 1)

    npl = SIZE(symba_plA, 1)

    symba_tpdA = PACK(symba_tpA,symba_tpA%helio%swiftest%status /= ACTIVE)
    symba_tpA = PACK(symba_tpA,symba_tpA%helio%swiftest%status = ACTIVE)
    nsptp = SIZE(symba_tpdA, 1)

    ntp = SIZE(symba_tpA, 1)









    !do the discarding of pl and tp 


    CALL symba_energy(npl, nplmax, swifter_pl1P, j2rp2, j4rp4, ke, pe, tef, htot)