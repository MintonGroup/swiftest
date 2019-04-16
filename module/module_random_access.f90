!**********************************************************************************************************************************
!
!  Unit Name   : module_random_access
!  Unit Type   : module
!  Project     : SWIFTER
!  Package     : module
!  Language    : Fortran 90/95
!
!  Description : Definition of data and procedures to enable random (indexed) access to planet and test particle arrays
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
MODULE module_random_access

     USE module_parameters
     USE module_swifter
     USE module_bs
     USE module_helio
     USE module_ra15
     USE module_tu4
     USE module_whm
     USE module_rmvs
     USE module_symba
     USE module_interfaces
     IMPLICIT NONE

     INTEGER(I4B)                            :: istruct_pl, istruct_tp
     TYPE(swifter_pl), DIMENSION(:), POINTER :: lswifter_pl1P
     TYPE(swifter_tp), DIMENSION(:), POINTER :: lswifter_tp1P
     TYPE(bs_pl), DIMENSION(:), POINTER      :: lbs_pl1P
     TYPE(bs_tp), DIMENSION(:), POINTER      :: lbs_tp1P
     TYPE(helio_pl), DIMENSION(:), POINTER   :: lhelio_pl1P
     TYPE(helio_tp), DIMENSION(:), POINTER   :: lhelio_tp1P
     TYPE(ra15_pl), DIMENSION(:), POINTER    :: lra15_pl1P
     TYPE(ra15_tp), DIMENSION(:), POINTER    :: lra15_tp1P
     TYPE(tu4_pl), DIMENSION(:), POINTER     :: ltu4_pl1P
     TYPE(tu4_tp), DIMENSION(:), POINTER     :: ltu4_tp1P
     TYPE(whm_pl), DIMENSION(:), POINTER     :: lwhm_pl1P
     TYPE(whm_tp), DIMENSION(:), POINTER     :: lwhm_tp1P
     TYPE(rmvs_pl), DIMENSION(:), POINTER    :: lrmvs_pl1P
     TYPE(rmvs_tp), DIMENSION(:), POINTER    :: lrmvs_tp1P
     TYPE(symba_pl), DIMENSION(:), POINTER   :: lsymba_pl1P
     TYPE(symba_tp), DIMENSION(:), POINTER   :: lsymba_tp1P

     INTERFACE set_point
          MODULE PROCEDURE set_point_swifter_pl
          MODULE PROCEDURE set_point_swifter_tp
          MODULE PROCEDURE set_point_bs_pl
          MODULE PROCEDURE set_point_bs_tp
          MODULE PROCEDURE set_point_helio_pl
          MODULE PROCEDURE set_point_helio_tp
          MODULE PROCEDURE set_point_ra15_pl
          MODULE PROCEDURE set_point_ra15_tp
          MODULE PROCEDURE set_point_tu4_pl
          MODULE PROCEDURE set_point_tu4_tp
          MODULE PROCEDURE set_point_whm_pl
          MODULE PROCEDURE set_point_whm_tp
          MODULE PROCEDURE set_point_rmvs_pl
          MODULE PROCEDURE set_point_rmvs_tp
          MODULE PROCEDURE set_point_symba_pl
          MODULE PROCEDURE set_point_symba_tp
     END INTERFACE

     INTERFACE get_point
          MODULE PROCEDURE get_point_swifter_pl
          MODULE PROCEDURE get_point_swifter_tp
          MODULE PROCEDURE get_point_bs_pl
          MODULE PROCEDURE get_point_bs_tp
          MODULE PROCEDURE get_point_helio_pl
          MODULE PROCEDURE get_point_helio_tp
          MODULE PROCEDURE get_point_ra15_pl
          MODULE PROCEDURE get_point_ra15_tp
          MODULE PROCEDURE get_point_tu4_pl
          MODULE PROCEDURE get_point_tu4_tp
          MODULE PROCEDURE get_point_whm_pl
          MODULE PROCEDURE get_point_whm_tp
          MODULE PROCEDURE get_point_rmvs_pl
          MODULE PROCEDURE get_point_rmvs_tp
          MODULE PROCEDURE get_point_symba_pl
          MODULE PROCEDURE get_point_symba_tp
     END INTERFACE

CONTAINS

     SUBROUTINE set_point_swifter_pl(s_in)
          TYPE(swifter_pl), DIMENSION(:), TARGET :: s_in
          istruct_pl = SWIFTER
          lswifter_pl1P => s_in
          RETURN
     END SUBROUTINE set_point_swifter_pl

     SUBROUTINE set_point_swifter_tp(s_in)
          TYPE(swifter_tp), DIMENSION(:), TARGET :: s_in
          istruct_tp = SWIFTER
          lswifter_tp1P => s_in
          RETURN
     END SUBROUTINE set_point_swifter_tp

     SUBROUTINE set_point_bs_pl(s_in)
          TYPE(bs_pl), DIMENSION(:), TARGET :: s_in
          istruct_pl = BS
          lbs_pl1P => s_in
          RETURN
     END SUBROUTINE set_point_bs_pl

     SUBROUTINE set_point_bs_tp(s_in)
          TYPE(bs_tp), DIMENSION(:), TARGET :: s_in
          istruct_tp = BS
          lbs_tp1P => s_in
          RETURN
     END SUBROUTINE set_point_bs_tp

     SUBROUTINE set_point_helio_pl(s_in)
          TYPE(helio_pl), DIMENSION(:), TARGET :: s_in
          istruct_pl = HELIO
          lhelio_pl1P => s_in
          RETURN
     END SUBROUTINE set_point_helio_pl

     SUBROUTINE set_point_helio_tp(s_in)
          TYPE(helio_tp), DIMENSION(:), TARGET :: s_in
          istruct_tp = HELIO
          lhelio_tp1P => s_in
          RETURN
     END SUBROUTINE set_point_helio_tp

     SUBROUTINE set_point_ra15_pl(s_in)
          TYPE(ra15_pl), DIMENSION(:), TARGET :: s_in
          istruct_pl = RA15
          lra15_pl1P => s_in
          RETURN
     END SUBROUTINE set_point_ra15_pl

     SUBROUTINE set_point_ra15_tp(s_in)
          TYPE(ra15_tp), DIMENSION(:), TARGET :: s_in
          istruct_tp = RA15
          lra15_tp1P => s_in
          RETURN
     END SUBROUTINE set_point_ra15_tp

     SUBROUTINE set_point_tu4_pl(s_in)
          TYPE(tu4_pl), DIMENSION(:), TARGET :: s_in
          istruct_pl = TU4
          ltu4_pl1P => s_in
          RETURN
     END SUBROUTINE set_point_tu4_pl

     SUBROUTINE set_point_tu4_tp(s_in)
          TYPE(tu4_tp), DIMENSION(:), TARGET :: s_in
          istruct_tp = TU4
          ltu4_tp1P => s_in
          RETURN
     END SUBROUTINE set_point_tu4_tp

     SUBROUTINE set_point_whm_pl(s_in)
          TYPE(whm_pl), DIMENSION(:), TARGET :: s_in
          istruct_pl = WHM
          lwhm_pl1P => s_in
          RETURN
     END SUBROUTINE set_point_whm_pl

     SUBROUTINE set_point_whm_tp(s_in)
          TYPE(whm_tp), DIMENSION(:), TARGET :: s_in
          istruct_tp = WHM
          lwhm_tp1P => s_in
          RETURN
     END SUBROUTINE set_point_whm_tp

     SUBROUTINE set_point_rmvs_pl(s_in)
          TYPE(rmvs_pl), DIMENSION(:), TARGET :: s_in
          istruct_pl = RMVS
          lrmvs_pl1P => s_in
          RETURN
     END SUBROUTINE set_point_rmvs_pl

     SUBROUTINE set_point_rmvs_tp(s_in)
          TYPE(rmvs_tp), DIMENSION(:), TARGET :: s_in
          istruct_tp = RMVS
          lrmvs_tp1P => s_in
          RETURN
     END SUBROUTINE set_point_rmvs_tp

     SUBROUTINE set_point_symba_pl(s_in)
          TYPE(symba_pl), DIMENSION(:), TARGET :: s_in
          istruct_pl = SYMBA
          lsymba_pl1P => s_in
          RETURN
     END SUBROUTINE set_point_symba_pl

     SUBROUTINE set_point_symba_tp(s_in)
          TYPE(symba_tp), DIMENSION(:), TARGET :: s_in
          istruct_tp = SYMBA
          lsymba_tp1P => s_in
          RETURN
     END SUBROUTINE set_point_symba_tp

     SUBROUTINE get_point_swifter_pl(i, swifter_plP)
          INTEGER(I4B), INTENT(IN)  :: i
          TYPE(swifter_pl), POINTER :: swifter_plP
          SELECT CASE (istruct_pl)
               CASE (SWIFTER)
                    IF (ASSOCIATED(lswifter_pl1P)) THEN
                         swifter_plP => lswifter_pl1P(i)
                    ELSE
                         WRITE(*, *) "Error in get_point_swifter_pl"
                         CALL util_exit(FAILURE)
                    END IF
               CASE (BS)
                    IF (ASSOCIATED(lbs_pl1P)) THEN
                         swifter_plP => lbs_pl1P(i)%swifter
                    ELSE
                         WRITE(*, *) "Error in get_point_swifter_pl"
                         CALL util_exit(FAILURE)
                    END IF
               CASE (HELIO)
                    IF (ASSOCIATED(lhelio_pl1P)) THEN
                         swifter_plP => lhelio_pl1P(i)%swifter
                    ELSE
                         WRITE(*, *) "Error in get_point_swifter_pl"
                         CALL util_exit(FAILURE)
                    END IF
               CASE (RA15)
                    IF (ASSOCIATED(lra15_pl1P)) THEN
                         swifter_plP => lra15_pl1P(i)%swifter
                    ELSE
                         WRITE(*, *) "Error in get_point_swifter_pl"
                         CALL util_exit(FAILURE)
                    END IF
               CASE (TU4)
                    IF (ASSOCIATED(ltu4_pl1P)) THEN
                         swifter_plP => ltu4_pl1P(i)%swifter
                    ELSE
                         WRITE(*, *) "Error in get_point_swifter_pl"
                         CALL util_exit(FAILURE)
                    END IF
               CASE (WHM)
                    IF (ASSOCIATED(lwhm_pl1P)) THEN
                         swifter_plP => lwhm_pl1P(i)%swifter
                    ELSE
                         WRITE(*, *) "Error in get_point_swifter_pl"
                         CALL util_exit(FAILURE)
                    END IF
               CASE (RMVS)
                    IF (ASSOCIATED(lrmvs_pl1P)) THEN
                         swifter_plP => lrmvs_pl1P(i)%whm%swifter
                    ELSE
                         WRITE(*, *) "Error in get_point_swifter_pl"
                         CALL util_exit(FAILURE)
                    END IF
               CASE (SYMBA)
                    IF (ASSOCIATED(lsymba_pl1P)) THEN
                         swifter_plP => lsymba_pl1P(i)%helio%swifter
                    ELSE
                         WRITE(*, *) "Error in get_point_swifter_pl"
                         CALL util_exit(FAILURE)
                    END IF
! DEK - improve error handler
               CASE DEFAULT
                    WRITE(*, *) "Error in get_point_swifter_pl"
                    CALL util_exit(FAILURE)
          END SELECT
          RETURN
     END SUBROUTINE get_point_swifter_pl

     SUBROUTINE get_point_swifter_tp(i, swifter_tpP)
          INTEGER(I4B), INTENT(IN)  :: i
          TYPE(swifter_tp), POINTER :: swifter_tpP
          SELECT CASE (istruct_tp)
               CASE (SWIFTER)
                    IF (ASSOCIATED(lswifter_tp1P)) THEN
                         swifter_tpP => lswifter_tp1P(i)
                    ELSE
                         WRITE(*, *) "Error in get_point_swifter_tp"
                         CALL util_exit(FAILURE)
                    END IF
               CASE (BS)
                    IF (ASSOCIATED(lbs_tp1P)) THEN
                         swifter_tpP => lbs_tp1P(i)%swifter
                    ELSE
                         WRITE(*, *) "Error in get_point_swifter_tp"
                         CALL util_exit(FAILURE)
                    END IF
               CASE (HELIO)
                    IF (ASSOCIATED(lhelio_tp1P)) THEN
                         swifter_tpP => lhelio_tp1P(i)%swifter
                    ELSE
                         WRITE(*, *) "Error in get_point_swifter_tp"
                         CALL util_exit(FAILURE)
                    END IF
               CASE (RA15)
                    IF (ASSOCIATED(lra15_tp1P)) THEN
                         swifter_tpP => lra15_tp1P(i)%swifter
                    ELSE
                         WRITE(*, *) "Error in get_point_swifter_tp"
                         CALL util_exit(FAILURE)
                    END IF
               CASE (TU4)
                    IF (ASSOCIATED(ltu4_tp1P)) THEN
                         swifter_tpP => ltu4_tp1P(i)%swifter
                    ELSE
                         WRITE(*, *) "Error in get_point_swifter_tp"
                         CALL util_exit(FAILURE)
                    END IF
               CASE (WHM)
                    IF (ASSOCIATED(lwhm_tp1P)) THEN
                         swifter_tpP => lwhm_tp1P(i)%swifter
                    ELSE
                         WRITE(*, *) "Error in get_point_swifter_tp"
                         CALL util_exit(FAILURE)
                    END IF
               CASE (RMVS)
                    IF (ASSOCIATED(lrmvs_tp1P)) THEN
                         swifter_tpP => lrmvs_tp1P(i)%whm%swifter
                    ELSE
                         WRITE(*, *) "Error in get_point_swifter_tp"
                         CALL util_exit(FAILURE)
                    END IF
               CASE (SYMBA)
                    IF (ASSOCIATED(lsymba_tp1P)) THEN
                         swifter_tpP => lsymba_tp1P(i)%helio%swifter
                    ELSE
                         WRITE(*, *) "Error in get_point_swifter_tp"
                         CALL util_exit(FAILURE)
                    END IF
! DEK - improve error handler
               CASE DEFAULT
                    WRITE(*, *) "Error in get_point_swifter_tp"
                    CALL util_exit(FAILURE)
          END SELECT
          RETURN
     END SUBROUTINE get_point_swifter_tp

     SUBROUTINE get_point_bs_pl(i, bs_plP)
          INTEGER(I4B), INTENT(IN) :: i
          TYPE(bs_pl), POINTER     :: bs_plP
          SELECT CASE (istruct_pl)
               CASE (BS)
                    IF (ASSOCIATED(lbs_pl1P)) THEN
                         bs_plP => lbs_pl1P(i)
                    ELSE
                         WRITE(*, *) "Error in get_point_bs_pl"
                         CALL util_exit(FAILURE)
                    END IF
! DEK - improve error handler
               CASE DEFAULT
                    WRITE(*, *) "Error in get_point_bs_pl"
                    CALL util_exit(FAILURE)
          END SELECT
          RETURN
     END SUBROUTINE get_point_bs_pl

     SUBROUTINE get_point_bs_tp(i, bs_tpP)
          INTEGER(I4B), INTENT(IN) :: i
          TYPE(bs_tp), POINTER     :: bs_tpP
          SELECT CASE (istruct_tp)
               CASE (BS)
                    IF (ASSOCIATED(lbs_tp1P)) THEN
                         bs_tpP => lbs_tp1P(i)
                    ELSE
                         WRITE(*, *) "Error in get_point_bs_tp"
                         CALL util_exit(FAILURE)
                    END IF
! DEK - improve error handler
               CASE DEFAULT
                    WRITE(*, *) "Error in get_point_bs_tp"
                    CALL util_exit(FAILURE)
          END SELECT
          RETURN
     END SUBROUTINE get_point_bs_tp

     SUBROUTINE get_point_helio_pl(i, helio_plP)
          INTEGER(I4B), INTENT(IN) :: i
          TYPE(helio_pl), POINTER  :: helio_plP
          SELECT CASE (istruct_pl)
               CASE (HELIO)
                    IF (ASSOCIATED(lhelio_pl1P)) THEN
                         helio_plP => lhelio_pl1P(i)
                    ELSE
                         WRITE(*, *) "Error in get_point_helio_pl"
                         CALL util_exit(FAILURE)
                    END IF
               CASE (SYMBA)
                    IF (ASSOCIATED(lsymba_pl1P)) THEN
                         helio_plP => lsymba_pl1P(i)%helio
                    ELSE
                         WRITE(*, *) "Error in get_point_helio_pl"
                         CALL util_exit(FAILURE)
                    END IF
! DEK - improve error handler
               CASE DEFAULT
                    WRITE(*, *) "Error in get_point_helio_pl"
                    CALL util_exit(FAILURE)
          END SELECT
          RETURN
     END SUBROUTINE get_point_helio_pl

     SUBROUTINE get_point_helio_tp(i, helio_tpP)
          INTEGER(I4B), INTENT(IN) :: i
          TYPE(helio_tp), POINTER  :: helio_tpP
          SELECT CASE (istruct_tp)
               CASE (HELIO)
                    IF (ASSOCIATED(lhelio_tp1P)) THEN
                         helio_tpP => lhelio_tp1P(i)
                    ELSE
                         WRITE(*, *) "Error in get_point_helio_tp"
                         CALL util_exit(FAILURE)
                    END IF
               CASE (SYMBA)
                    IF (ASSOCIATED(lsymba_tp1P)) THEN
                         helio_tpP => lsymba_tp1P(i)%helio
                    ELSE
                         WRITE(*, *) "Error in get_point_helio_tp"
                         CALL util_exit(FAILURE)
                    END IF
! DEK - improve error handler
               CASE DEFAULT
                    WRITE(*, *) "Error in get_point_helio_tp"
                    CALL util_exit(FAILURE)
          END SELECT
          RETURN
     END SUBROUTINE get_point_helio_tp

     SUBROUTINE get_point_ra15_pl(i, ra15_plP)
          INTEGER(I4B), INTENT(IN) :: i
          TYPE(ra15_pl), POINTER   :: ra15_plP
          SELECT CASE (istruct_pl)
               CASE (RA15)
                    IF (ASSOCIATED(lra15_pl1P)) THEN
                         ra15_plP => lra15_pl1P(i)
                    ELSE
                         WRITE(*, *) "Error in get_point_ra15_pl"
                         CALL util_exit(FAILURE)
                    END IF
! DEK - improve error handler
               CASE DEFAULT
                    WRITE(*, *) "Error in get_point_ra15_pl"
                    CALL util_exit(FAILURE)
          END SELECT
          RETURN
     END SUBROUTINE get_point_ra15_pl

     SUBROUTINE get_point_ra15_tp(i, ra15_tpP)
          INTEGER(I4B), INTENT(IN) :: i
          TYPE(ra15_tp), POINTER   :: ra15_tpP
          SELECT CASE (istruct_tp)
               CASE (RA15)
                    IF (ASSOCIATED(lra15_tp1P)) THEN
                         ra15_tpP => lra15_tp1P(i)
                    ELSE
                         WRITE(*, *) "Error in get_point_ra15_tp"
                         CALL util_exit(FAILURE)
                    END IF
! DEK - improve error handler
               CASE DEFAULT
                    WRITE(*, *) "Error in get_point_ra15_tp"
                    CALL util_exit(FAILURE)
          END SELECT
          RETURN
     END SUBROUTINE get_point_ra15_tp

     SUBROUTINE get_point_tu4_pl(i, tu4_plP)
          INTEGER(I4B), INTENT(IN) :: i
          TYPE(tu4_pl), POINTER    :: tu4_plP
          SELECT CASE (istruct_pl)
               CASE (TU4)
                    IF (ASSOCIATED(ltu4_pl1P)) THEN
                         tu4_plP => ltu4_pl1P(i)
                    ELSE
                         WRITE(*, *) "Error in get_point_tu4_pl"
                         CALL util_exit(FAILURE)
                    END IF
! DEK - improve error handler
               CASE DEFAULT
                    WRITE(*, *) "Error in get_point_tu4_pl"
                    CALL util_exit(FAILURE)
          END SELECT
          RETURN
     END SUBROUTINE get_point_tu4_pl

     SUBROUTINE get_point_tu4_tp(i, tu4_tpP)
          INTEGER(I4B), INTENT(IN) :: i
          TYPE(tu4_tp), POINTER    :: tu4_tpP
          SELECT CASE (istruct_tp)
               CASE (TU4)
                    IF (ASSOCIATED(ltu4_tp1P)) THEN
                         tu4_tpP => ltu4_tp1P(i)
                    ELSE
                         WRITE(*, *) "Error in get_point_tu4_tp"
                         CALL util_exit(FAILURE)
                    END IF
! DEK - improve error handler
               CASE DEFAULT
                    WRITE(*, *) "Error in get_point_tu4_tp"
                    CALL util_exit(FAILURE)
          END SELECT
          RETURN
     END SUBROUTINE get_point_tu4_tp

     SUBROUTINE get_point_whm_pl(i, whm_plP)
          INTEGER(I4B), INTENT(IN) :: i
          TYPE(whm_pl), POINTER    :: whm_plP
          SELECT CASE (istruct_pl)
               CASE (WHM)
                    IF (ASSOCIATED(lwhm_pl1P)) THEN
                         whm_plP => lwhm_pl1P(i)
                    ELSE
                         WRITE(*, *) "Error in get_point_whm_pl"
                         CALL util_exit(FAILURE)
                    END IF
               CASE (RMVS)
                    IF (ASSOCIATED(lrmvs_pl1P)) THEN
                         whm_plP => lrmvs_pl1P(i)%whm
                    ELSE
                         WRITE(*, *) "Error in get_point_whm_pl"
                         CALL util_exit(FAILURE)
                    END IF
! DEK - improve error handler
               CASE DEFAULT
                    WRITE(*, *) "Error in get_point_whm_pl"
                    CALL util_exit(FAILURE)
          END SELECT
          RETURN
     END SUBROUTINE get_point_whm_pl

     SUBROUTINE get_point_whm_tp(i, whm_tpP)
          INTEGER(I4B), INTENT(IN) :: i
          TYPE(whm_tp), POINTER    :: whm_tpP
          SELECT CASE (istruct_tp)
               CASE (WHM)
                    IF (ASSOCIATED(lwhm_tp1P)) THEN
                         whm_tpP => lwhm_tp1P(i)
                    ELSE
                         WRITE(*, *) "Error in get_point_whm_tp"
                         CALL util_exit(FAILURE)
                    END IF
               CASE (RMVS)
                    IF (ASSOCIATED(lrmvs_tp1P)) THEN
                         whm_tpP => lrmvs_tp1P(i)%whm
                    ELSE
                         WRITE(*, *) "Error in get_point_whm_tp"
                         CALL util_exit(FAILURE)
                    END IF
! DEK - improve error handler
               CASE DEFAULT
                    WRITE(*, *) "Error in get_point_whm_tp"
                    CALL util_exit(FAILURE)
          END SELECT
          RETURN
     END SUBROUTINE get_point_whm_tp

     SUBROUTINE get_point_rmvs_pl(i, rmvs_plP)
          INTEGER(I4B), INTENT(IN) :: i
          TYPE(rmvs_pl), POINTER   :: rmvs_plP
          SELECT CASE (istruct_pl)
               CASE (RMVS)
                    IF (ASSOCIATED(lrmvs_pl1P)) THEN
                         rmvs_plP => lrmvs_pl1P(i)
                    ELSE
                         WRITE(*, *) "Error in get_point_rmvs_pl"
                         CALL util_exit(FAILURE)
                    END IF
! DEK - improve error handler
               CASE DEFAULT
                    WRITE(*, *) "Error in get_point_rmvs_pl"
                    CALL util_exit(FAILURE)
          END SELECT
          RETURN
     END SUBROUTINE get_point_rmvs_pl

     SUBROUTINE get_point_rmvs_tp(i, rmvs_tpP)
          INTEGER(I4B), INTENT(IN) :: i
          TYPE(rmvs_tp), POINTER   :: rmvs_tpP
          SELECT CASE (istruct_tp)
               CASE (RMVS)
                    IF (ASSOCIATED(lrmvs_tp1P)) THEN
                         rmvs_tpP => lrmvs_tp1P(i)
                    ELSE
                         WRITE(*, *) "Error in get_point_rmvs_tp"
                         CALL util_exit(FAILURE)
                    END IF
! DEK - improve error handler
               CASE DEFAULT
                    WRITE(*, *) "Error in get_point_rmvs_tp"
                    CALL util_exit(FAILURE)
          END SELECT
          RETURN
     END SUBROUTINE get_point_rmvs_tp

     SUBROUTINE get_point_symba_pl(i, symba_plP)
          INTEGER(I4B), INTENT(IN) :: i
          TYPE(symba_pl), POINTER  :: symba_plP
          SELECT CASE (istruct_pl)
               CASE (SYMBA)
                    IF (ASSOCIATED(lsymba_pl1P)) THEN
                         symba_plP => lsymba_pl1P(i)
                    ELSE
                         WRITE(*, *) "Error in get_point_symba_pl"
                         CALL util_exit(FAILURE)
                    END IF
! DEK - improve error handler
               CASE DEFAULT
                    WRITE(*, *) "Error in get_point_symba_pl"
                    CALL util_exit(FAILURE)
          END SELECT
          RETURN
     END SUBROUTINE get_point_symba_pl

     SUBROUTINE get_point_symba_tp(i, symba_tpP)
          INTEGER(I4B), INTENT(IN) :: i
          TYPE(symba_tp), POINTER  :: symba_tpP
          SELECT CASE (istruct_tp)
               CASE (SYMBA)
                    IF (ASSOCIATED(lsymba_tp1P)) THEN
                         symba_tpP => lsymba_tp1P(i)
                    ELSE
                         WRITE(*, *) "Error in get_point_symba_tp"
                         CALL util_exit(FAILURE)
                    END IF
! DEK - improve error handler
               CASE DEFAULT
                    WRITE(*, *) "Error in get_point_symba_tp"
                    CALL util_exit(FAILURE)
          END SELECT
          RETURN
     END SUBROUTINE get_point_symba_tp

END MODULE module_random_access
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
