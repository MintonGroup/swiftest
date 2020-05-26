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
SUBROUTINE symba_rearray(t, npl, ntp, nsppl, nsptp, symba_plA, symba_tpA, nmergeadd, mergeadd_list, discard_plA, &
    discard_tpA, NPLMAX, j2rp2, j4rp4)

! Modules
     USE module_parameters
     USE module_swiftestalloc 
     USE module_swiftest
     USE module_helio
     USE module_symba
     USE module_interfaces, EXCEPT_THIS_ONE => symba_rearray
     IMPLICIT NONE

! Arguments
     INTEGER(I4B), INTENT(INOUT)                      :: npl, ntp, nsppl, nsptp, nmergeadd, NPLMAX !change to fragadd
     REAL(DP), INTENT(IN)                             :: t, j2rp2, j4rp4
     TYPE(symba_pl), INTENT(INOUT)                    :: symba_plA
     TYPE(symba_tp), INTENT(INOUT)                    :: symba_tpA
     TYPE(swiftest_tp), INTENT(INOUT)                 :: discard_tpA
     TYPE(swiftest_pl), INTENT(INOUT)                 :: discard_plA
     TYPE(symba_merger), INTENT(INOUT)                :: mergeadd_list !change to fragadd_list

! Internals
     INTEGER(I4B)                                   :: i, index, j, ncomp, ierr, nplm, nkpl, nktp, k, nfrag
     REAL(DP)                                       :: ke, pe, tei, tef
     REAL(DP), DIMENSION(NDIM)                      :: htot
     LOGICAL, DIMENSION(npl)                        :: discard_l_pl, frag_l_add
     LOGICAL, DIMENSION(ntp)                        :: discard_l_tp

! Executable code

    IF (ldiscard .eqv. .TRUE.) THEN 
        nsppl = 0
        nkpl = 0
        discard_l_pl(1:npl) = (symba_plA%helio%swiftest%status(1:npl) /= ACTIVE) 
        nsppl = COUNT(discard_l_pl)
        nkpl = npl - nsppl
        IF (lfragmentation) THEN
            frag_l_add(1:npl) = (mergeadd_list%status(1:npl) == (DISRUPTION .OR. HIT_AND_RUN .OR. SUPERCATASTROPHIC))
            nfrag = COUNT(frag_l_add)
        END IF

        CALL swiftest_pl_allocate(discard_plA,nsppl)

        discard_plA%name(1:nsppl) = PACK(symba_plA%helio%swiftest%name(1:npl), discard_l_pl)
        discard_plA%status(1:nsppl) = PACK(symba_plA%helio%swiftest%status(1:npl), discard_l_pl)
        discard_plA%mass(1:nsppl) = PACK(symba_plA%helio%swiftest%mass(1:npl), discard_l_pl)
        discard_plA%radius(1:nsppl) = PACK(symba_plA%helio%swiftest%radius(1:npl), discard_l_pl)
        discard_plA%xh(1,1:nsppl) = PACK(symba_plA%helio%swiftest%xh(1,1:npl), discard_l_pl)
        discard_plA%xh(2,1:nsppl) = PACK(symba_plA%helio%swiftest%xh(2,1:npl), discard_l_pl)
        discard_plA%xh(3,1:nsppl) = PACK(symba_plA%helio%swiftest%xh(3,1:npl), discard_l_pl)
        discard_plA%vh(1,1:nsppl) = PACK(symba_plA%helio%swiftest%vh(1,1:npl), discard_l_pl)
        discard_plA%vh(2,1:nsppl) = PACK(symba_plA%helio%swiftest%vh(2,1:npl), discard_l_pl)
        discard_plA%vh(3,1:nsppl) = PACK(symba_plA%helio%swiftest%vh(3,1:npl), discard_l_pl)
        discard_plA%rhill(1:nsppl) = PACK(symba_plA%helio%swiftest%rhill(1:npl), discard_l_pl)
        discard_plA%xb(1,1:nsppl) = PACK(symba_plA%helio%swiftest%xb(1,1:npl), discard_l_pl)
        discard_plA%xb(2,1:nsppl) = PACK(symba_plA%helio%swiftest%xb(2,1:npl), discard_l_pl)
        discard_plA%xb(3,1:nsppl) = PACK(symba_plA%helio%swiftest%xb(3,1:npl), discard_l_pl)
        discard_plA%vb(1,1:nsppl) = PACK(symba_plA%helio%swiftest%vb(1,1:npl), discard_l_pl)
        discard_plA%vb(2,1:nsppl) = PACK(symba_plA%helio%swiftest%vb(2,1:npl), discard_l_pl)
        discard_plA%vb(3,1:nsppl) = PACK(symba_plA%helio%swiftest%vb(3,1:npl), discard_l_pl)
        IF (lfragmentation .AND. (nkpl + nfrag > npl)) THEN 
            symba_plA%helio%swiftest%name(1:nkpl) = PACK(symba_plA%helio%swiftest%name(1:npl), .NOT. discard_l_pl)
            symba_plA%helio%swiftest%status(1:nkpl) = PACK(symba_plA%helio%swiftest%status(1:npl), .NOT. discard_l_pl)
            symba_plA%helio%swiftest%mass(1:nkpl) = PACK(symba_plA%helio%swiftest%mass(1:npl), .NOT. discard_l_pl)
            symba_plA%helio%swiftest%radius(1:nkpl) = PACK(symba_plA%helio%swiftest%radius(1:npl), .NOT. discard_l_pl)
            symba_plA%helio%swiftest%xh(1,1:nkpl) = PACK(symba_plA%helio%swiftest%xh(1,1:npl), .NOT. discard_l_pl)
            symba_plA%helio%swiftest%xh(2,1:nkpl) = PACK(symba_plA%helio%swiftest%xh(2,1:npl), .NOT. discard_l_pl)
            symba_plA%helio%swiftest%xh(3,1:nkpl) = PACK(symba_plA%helio%swiftest%xh(3,1:npl), .NOT. discard_l_pl)
            symba_plA%helio%swiftest%vh(1,1:nkpl) = PACK(symba_plA%helio%swiftest%vh(1,1:npl), .NOT. discard_l_pl)
            symba_plA%helio%swiftest%vh(2,1:nkpl) = PACK(symba_plA%helio%swiftest%vh(2,1:npl), .NOT. discard_l_pl)
            symba_plA%helio%swiftest%vh(3,1:nkpl) = PACK(symba_plA%helio%swiftest%vh(3,1:npl), .NOT. discard_l_pl)
            symba_plA%helio%swiftest%rhill(1:nkpl) = PACK(symba_plA%helio%swiftest%rhill(1:npl), .NOT. discard_l_pl)
            symba_plA%helio%swiftest%xb(1,1:nkpl) = PACK(symba_plA%helio%swiftest%xb(1,1:npl), .NOT. discard_l_pl)
            symba_plA%helio%swiftest%xb(2,1:nkpl) = PACK(symba_plA%helio%swiftest%xb(2,1:npl), .NOT. discard_l_pl)
            symba_plA%helio%swiftest%xb(3,1:nkpl) = PACK(symba_plA%helio%swiftest%xb(3,1:npl), .NOT. discard_l_pl)
            symba_plA%helio%swiftest%vb(1,1:nkpl) = PACK(symba_plA%helio%swiftest%vb(1,1:npl), .NOT. discard_l_pl)
            symba_plA%helio%swiftest%vb(2,1:nkpl) = PACK(symba_plA%helio%swiftest%vb(2,1:npl), .NOT. discard_l_pl)
            symba_plA%helio%swiftest%vb(3,1:nkpl) = PACK(symba_plA%helio%swiftest%vb(3,1:npl), .NOT. discard_l_pl)
            symba_plA%helio%ah(1,1:nkpl) = PACK(symba_plA%helio%ah(1,1:npl), .NOT. discard_l_pl)
            symba_plA%helio%ah(2,1:nkpl) = PACK(symba_plA%helio%ah(2,1:npl), .NOT. discard_l_pl)
            symba_plA%helio%ah(3,1:nkpl) =PACK(symba_plA%helio%ah(3,1:npl), .NOT. discard_l_pl)
            CALL util_resize_pl(symba_plA, nkpl+nfrag, npl)
            npl = nkpl  + nfrag
            !add fragments 
            symba_plA%helio%swiftest%name(nkpl+1:npl) = PACK(mergeadd_list%name(:), frag_l_add)
            symba_plA%helio%swiftest%status(nkpl+1:npl) = [(ACTIVE,i=1,nfrag)]!array of ACTIVE STATUS 
            symba_plA%helio%swiftest%mass(nkpl+1:npl) = PACK(mergeadd_list%mass(:), frag_l_add)
            symba_plA%helio%swiftest%radius(nkpl+1:npl) = PACK(mergeadd_list%radius(:), frag_l_add)
            symba_plA%helio%swiftest%xh(1,nkpl+1:npl) = PACK(mergeadd_list%xh(1,:), frag_l_add)
            symba_plA%helio%swiftest%xh(2,nkpl+1:npl) = PACK(mergeadd_list%xh(2,:), frag_l_add)
            symba_plA%helio%swiftest%xh(3,nkpl+1:npl) = PACK(mergeadd_list%xh(3,:), frag_l_add)
            symba_plA%helio%swiftest%vh(1,nkpl+1:npl) = PACK(mergeadd_list%vh(1,:), frag_l_add)
            symba_plA%helio%swiftest%vh(2,nkpl+1:npl) = PACK(mergeadd_list%vh(2,:), frag_l_add)
            symba_plA%helio%swiftest%vh(3,nkpl+1:npl) = PACK(mergeadd_list%vh(3,:), frag_l_add)

    
        ELSE
            symba_plA%helio%swiftest%name(1:nkpl) = PACK(symba_plA%helio%swiftest%name(1:npl), .NOT. discard_l_pl)
            symba_plA%helio%swiftest%status(1:nkpl) = PACK(symba_plA%helio%swiftest%status(1:npl), .NOT. discard_l_pl)
            symba_plA%helio%swiftest%mass(1:nkpl) = PACK(symba_plA%helio%swiftest%mass(1:npl), .NOT. discard_l_pl)
            symba_plA%helio%swiftest%radius(1:nkpl) = PACK(symba_plA%helio%swiftest%radius(1:npl), .NOT. discard_l_pl)
            symba_plA%helio%swiftest%xh(1,1:nkpl) = PACK(symba_plA%helio%swiftest%xh(1,1:npl), .NOT. discard_l_pl)
            symba_plA%helio%swiftest%xh(2,1:nkpl) = PACK(symba_plA%helio%swiftest%xh(2,1:npl), .NOT. discard_l_pl)
            symba_plA%helio%swiftest%xh(3,1:nkpl) = PACK(symba_plA%helio%swiftest%xh(3,1:npl), .NOT. discard_l_pl)
            symba_plA%helio%swiftest%vh(1,1:nkpl) = PACK(symba_plA%helio%swiftest%vh(1,1:npl), .NOT. discard_l_pl)
            symba_plA%helio%swiftest%vh(2,1:nkpl) = PACK(symba_plA%helio%swiftest%vh(2,1:npl), .NOT. discard_l_pl)
            symba_plA%helio%swiftest%vh(3,1:nkpl) = PACK(symba_plA%helio%swiftest%vh(3,1:npl), .NOT. discard_l_pl)
            symba_plA%helio%swiftest%rhill(1:nkpl) = PACK(symba_plA%helio%swiftest%rhill(1:npl), .NOT. discard_l_pl)
            symba_plA%helio%swiftest%xb(1,1:nkpl) = PACK(symba_plA%helio%swiftest%xb(1,1:npl), .NOT. discard_l_pl)
            symba_plA%helio%swiftest%xb(2,1:nkpl) = PACK(symba_plA%helio%swiftest%xb(2,1:npl), .NOT. discard_l_pl)
            symba_plA%helio%swiftest%xb(3,1:nkpl) = PACK(symba_plA%helio%swiftest%xb(3,1:npl), .NOT. discard_l_pl)
            symba_plA%helio%swiftest%vb(1,1:nkpl) = PACK(symba_plA%helio%swiftest%vb(1,1:npl), .NOT. discard_l_pl)
            symba_plA%helio%swiftest%vb(2,1:nkpl) = PACK(symba_plA%helio%swiftest%vb(2,1:npl), .NOT. discard_l_pl)
            symba_plA%helio%swiftest%vb(3,1:nkpl) = PACK(symba_plA%helio%swiftest%vb(3,1:npl), .NOT. discard_l_pl)
            symba_plA%helio%ah(1,1:nkpl) = PACK(symba_plA%helio%ah(1,1:npl), .NOT. discard_l_pl)
            symba_plA%helio%ah(2,1:nkpl) = PACK(symba_plA%helio%ah(2,1:npl), .NOT. discard_l_pl)
            symba_plA%helio%ah(3,1:nkpl) = PACK(symba_plA%helio%ah(3,1:npl), .NOT. discard_l_pl)
            npl = nkpl
        END IF
    END IF 

    IF (ldiscard_tp .eqv. .TRUE.) THEN 
        nktp = 0
        nsptp = 0  

        discard_l_tp(1:ntp) = (symba_tpA%helio%swiftest%status(1:ntp) /= ACTIVE)
        nsptp = COUNT(discard_l_tp)
        nktp = ntp - nsptp

        CALL swiftest_tp_allocate(discard_tpA,nsptp)

        discard_tpA%name(1:nsptp) = PACK(symba_tpA%helio%swiftest%name(1:ntp), discard_l_tp)
        discard_tpA%status(1:nsptp) = PACK(symba_tpA%helio%swiftest%status(1:ntp), discard_l_tp)
        discard_tpA%xh(1,1:nsptp) = PACK(symba_tpA%helio%swiftest%xh(1,1:ntp), discard_l_tp)
        discard_tpA%xh(2,1:nsptp) = PACK(symba_tpA%helio%swiftest%xh(2,1:ntp), discard_l_tp)
        discard_tpA%xh(3,1:nsptp) = PACK(symba_tpA%helio%swiftest%xh(3,1:ntp), discard_l_tp)
        discard_tpA%vh(1,1:nsptp) = PACK(symba_tpA%helio%swiftest%vh(1,1:ntp), discard_l_tp)
        discard_tpA%vh(2,1:nsptp) = PACK(symba_tpA%helio%swiftest%vh(2,1:ntp), discard_l_tp)
        discard_tpA%vh(3,1:nsptp) = PACK(symba_tpA%helio%swiftest%vh(3,1:ntp), discard_l_tp)
        discard_tpA%isperi(1:nsptp) = PACK(symba_tpA%helio%swiftest%isperi(1:ntp), discard_l_tp)
        discard_tpA%peri(1:nsptp) = PACK(symba_tpA%helio%swiftest%peri(1:ntp), discard_l_tp)
        discard_tpA%atp(1:nsptp) = PACK(symba_tpA%helio%swiftest%atp(1:ntp), discard_l_tp)
        discard_tpA%xb(1,1:nsptp) = PACK(symba_tpA%helio%swiftest%xb(1,1:ntp), discard_l_tp)
        discard_tpA%xb(2,1:nsptp) = PACK(symba_tpA%helio%swiftest%xb(2,1:ntp), discard_l_tp)
        discard_tpA%xb(3,1:nsptp) = PACK(symba_tpA%helio%swiftest%xb(3,1:ntp), discard_l_tp)
        discard_tpA%vb(1,1:nsptp) = PACK(symba_tpA%helio%swiftest%vb(1,1:ntp), discard_l_tp)
        discard_tpA%vb(2,1:nsptp) = PACK(symba_tpA%helio%swiftest%vb(2,1:ntp), discard_l_tp)
        discard_tpA%vb(3,1:nsptp) = PACK(symba_tpA%helio%swiftest%vb(3,1:ntp), discard_l_tp)

        symba_tpA%helio%swiftest%name(1:nktp) = PACK(symba_tpA%helio%swiftest%name(1:ntp), .NOT. discard_l_tp)
        symba_tpA%helio%swiftest%status(1:nktp) = PACK(symba_tpA%helio%swiftest%status(1:ntp), .NOT. discard_l_tp)
        symba_tpA%helio%swiftest%xh(1,1:nktp) = PACK(symba_tpA%helio%swiftest%xh(1,1:ntp), .NOT. discard_l_tp)
        symba_tpA%helio%swiftest%xh(2,1:nktp) = PACK(symba_tpA%helio%swiftest%xh(2,1:ntp), .NOT. discard_l_tp)
        symba_tpA%helio%swiftest%xh(3,1:nktp) = PACK(symba_tpA%helio%swiftest%xh(3,1:ntp), .NOT. discard_l_tp)
        symba_tpA%helio%swiftest%vh(1,1:nktp) = PACK(symba_tpA%helio%swiftest%vh(1,1:ntp), .NOT. discard_l_tp)
        symba_tpA%helio%swiftest%vh(2,1:nktp) = PACK(symba_tpA%helio%swiftest%vh(2,1:ntp), .NOT. discard_l_tp)
        symba_tpA%helio%swiftest%vh(3,1:nktp) = PACK(symba_tpA%helio%swiftest%vh(3,1:ntp), .NOT. discard_l_tp)
        symba_tpA%helio%swiftest%xb(1,1:nktp) = PACK(symba_tpA%helio%swiftest%xb(1,1:ntp), .NOT. discard_l_tp)
        symba_tpA%helio%swiftest%xb(2,1:nktp) = PACK(symba_tpA%helio%swiftest%xb(2,1:ntp), .NOT. discard_l_tp)
        symba_tpA%helio%swiftest%xb(3,1:nktp) = PACK(symba_tpA%helio%swiftest%xb(3,1:ntp), .NOT. discard_l_tp)
        symba_tpA%helio%swiftest%vb(1,1:nktp) = PACK(symba_tpA%helio%swiftest%vb(1,1:ntp), .NOT. discard_l_tp)
        symba_tpA%helio%swiftest%vb(2,1:nktp) = PACK(symba_tpA%helio%swiftest%vb(2,1:ntp), .NOT. discard_l_tp)
        symba_tpA%helio%swiftest%vb(3,1:nktp) = PACK(symba_tpA%helio%swiftest%vb(3,1:ntp), .NOT. discard_l_tp)
        symba_tpA%helio%swiftest%isperi(1:nktp) = PACK(symba_tpA%helio%swiftest%isperi(1:ntp), .NOT. discard_l_tp)
        symba_tpA%helio%swiftest%peri(1:nktp) = PACK(symba_tpA%helio%swiftest%peri(1:ntp), .NOT. discard_l_tp)
        symba_tpA%helio%swiftest%atp(1:nktp) = PACK(symba_tpA%helio%swiftest%atp(1:ntp), .NOT. discard_l_tp)
        symba_tpA%helio%ah(1,1:nktp) = PACK(symba_tpA%helio%ah(1,1:ntp), .NOT. discard_l_tp)
        symba_tpA%helio%ah(2,1:nktp) = PACK(symba_tpA%helio%ah(2,1:ntp), .NOT. discard_l_tp)
        symba_tpA%helio%ah(3,1:nktp) = PACK(symba_tpA%helio%ah(3,1:ntp), .NOT. discard_l_tp)
        ntp = nktp
    END IF 

END SUBROUTINE symba_rearray


    


