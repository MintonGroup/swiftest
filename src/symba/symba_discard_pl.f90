!**********************************************************************************************************************************
!
!  Unit Name   : symba_discard_pl
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
!                nplmax      : maximum allowed number of planets
!                nsp         : number of spilled planets
!                symba_pl1P  : pointer to head of SyMBA planet structure linked-list
!                symba_pld1P : pointer to head of discard SyMBA planet structure linked-list
!                rmin        : minimum allowed heliocentric radius
!                rmax        : maximum allowed heliocentric radius
!                rmaxu       : maximum allowed heliocentric radius for unbound planets
!                qmin        : minimum allowed pericenter distance
!                qmin_coord  : coordinate frame for qmin
!                qmin_alo    : minimum semimajor axis for qmin
!                qmin_ahi    : maximum semimajor axis for qmin
!                j2rp2       : J2 * R**2 for the Sun
!                j4rp4       : J4 * R**4 for the Sun
!                eoffset     : energy offset
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : npl         : number of planets
!                nsp         : number of spilled planets
!                symba_pl1P  : pointer to head of SyMBA planet structure linked-list
!                symba_pld1P : pointer to head of discard SyMBA planet structure linked-list
!                eoffset     : energy offset
!    Terminal  : none
!    File      : none
!
!  Invocation  : CALL symba_discard_pl(t, npl, nplmax, nsp, symba_pl1P, symba_pld1P, rmin, rmax, rmaxu, qmin, qmin_coord,
!                                      qmin_alo, qmin_ahi, j2rp2, j4rp4, eoffset)
!
!  Notes       : Adapted from Hal Levison's Swift routine discard_massive5.f
!
!**********************************************************************************************************************************
SUBROUTINE symba_discard_pl(t, npl, nplmax, nsp, symba_plA, rmin, rmax, rmaxu, qmin, qmin_coord, qmin_alo,          &
     qmin_ahi, j2rp2, j4rp4, eoffset)

! Modules
     USE swiftest
     USE module_helio
     USE module_symba
     USE module_interfaces, EXCEPT_THIS_ONE => symba_discard_pl
     IMPLICIT NONE

! Arguments
     INTEGER(I4B), INTENT(IN)    :: nplmax
     INTEGER(I4B), INTENT(INOUT) :: npl, nsp
     REAL(DP), INTENT(IN)        :: t, rmin, rmax, rmaxu, qmin, qmin_alo, qmin_ahi, j2rp2, j4rp4
     REAL(DP), INTENT(INOUT)     :: eoffset
     CHARACTER(*), INTENT(IN)    :: qmin_coord
     TYPE(symba_pl)              :: symba_plA

! Internals
     LOGICAL(LGT)              :: ldiscards
     INTEGER(I4B)              :: i
     REAL(DP)                  :: msys, ke, pe, tei, tef
     REAL(DP), DIMENSION(NDIM) :: htot

! Executable code
     ldiscards = .FALSE.
     IF ((rmin >= 0.0_DP) .OR. (rmax >= 0.0_DP) .OR. (rmaxu >= 0.0_DP) .OR. ((qmin >= 0.0_DP) .AND. (qmin_coord == "BARY")))      &
          CALL coord_h2b(npl, symba_plA%helio%swiftest, msys)
     IF ((rmin >= 0.0_DP) .OR. (rmax >= 0.0_DP) .OR. (rmaxu >= 0.0_DP))                                                           &
          CALL symba_discard_sun_pl(t, npl, msys, symba_plA%helio%swiftest, rmin, rmax, rmaxu, ldiscards)
     IF (qmin >= 0.0_DP) CALL symba_discard_peri_pl(t, npl, symba_plA, msys, qmin, qmin_alo, qmin_ahi, qmin_coord, ldiscards)

     RETURN

END SUBROUTINE symba_discard_pl
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
