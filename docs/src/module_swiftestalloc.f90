!**********************************************************************************************************************************
!
!  Unit Name   : module_swiftestalloc
!  Unit Type   : module
!  Project     : SWIFTEST
!  Package     : module
!  Language    : Fortran 2003
!
!  Description : 
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
!
!  Author(s)   : Jennifer Pouplin & Carlisle Wisahrd
!
!**********************************************************************************************************************************
MODULE module_swiftestalloc

    USE swiftest
    IMPLICIT NONE

    CONTAINS 

        SUBROUTINE swiftest_pl_allocate(swiftest_plA, npl)
            USE swiftest
            USE module_swiftest
            IMPLICIT NONE

            ! Arguments
            INTEGER(I4B), INTENT(IN)            :: npl
            TYPE(swiftest_pl), INTENT(INOUT)    :: swiftest_plA

            ALLOCATE(swiftest_plA%name(npl))
            ALLOCATE(swiftest_plA%status(npl))
            ALLOCATE(swiftest_plA%mass(npl))
            ALLOCATE(swiftest_plA%radius(npl))
            ALLOCATE(swiftest_plA%rhill(npl))
            ALLOCATE(swiftest_plA%xh(NDIM,npl))
            ALLOCATE(swiftest_plA%vh(NDIM,npl))
            ALLOCATE(swiftest_plA%xb(NDIM,npl))
            ALLOCATE(swiftest_plA%vb(NDIM,npl))
            RETURN
        END SUBROUTINE swiftest_pl_allocate


        SUBROUTINE helio_pl_allocate(helio_plA, npl)
            USE swiftest
            USE module_helio
            IMPLICIT NONE

            ! Arguments
            INTEGER(I4B), INTENT(IN)            :: npl
            TYPE(helio_pl), INTENT(INOUT)        :: helio_plA

            ALLOCATE(helio_plA%ah(NDIM,npl))
             ALLOCATE(helio_plA%ahi(NDIM,npl))
             CALL swiftest_pl_allocate(helio_plA%swiftest,npl)
            return
        END SUBROUTINE helio_pl_allocate


        SUBROUTINE symba_pl_allocate(symba_plA, npl)
            USE swiftest
            USE module_symba
            IMPLICIT NONE

            ! Arguments
            INTEGER(I4B), INTENT(IN)            :: npl
            TYPE(symba_pl), INTENT(INOUT)        :: symba_plA
            ALLOCATE(symba_plA%lmerged(npl))
            ALLOCATE(symba_plA%nplenc(npl))
            ALLOCATE(symba_plA%ntpenc(npl))
            ALLOCATE(symba_plA%levelg(npl))
            ALLOCATE(symba_plA%levelm(npl))
            ALLOCATE(symba_plA%nchild(npl))
            ALLOCATE(symba_plA%isperi(npl))
            ALLOCATE(symba_plA%peri(npl))
            ALLOCATE(symba_plA%atp(npl))
            ALLOCATE(symba_plA%index_parent(npl))
            ALLOCATE(symba_plA%index_child(npl,npl))
            CALL helio_pl_allocate(symba_plA%helio,npl)
            return
        END SUBROUTINE symba_pl_allocate

        SUBROUTINE symba_plplenc_allocate(plplenc_list, nplplenc)
            USE swiftest
            USE module_symba
            IMPLICIT NONE

            ! Arguments
            INTEGER(I4B), INTENT(IN)                :: nplplenc
            TYPE(symba_plplenc), INTENT(INOUT)        :: plplenc_list

            ALLOCATE(plplenc_list%lvdotr(nplplenc))
            ALLOCATE(plplenc_list%status(nplplenc))
            ALLOCATE(plplenc_list%level(nplplenc))
            ALLOCATE(plplenc_list%index1(nplplenc))
            ALLOCATE(plplenc_list%index2(nplplenc))
            ALLOCATE(plplenc_list%enc_child(nplplenc))
            ALLOCATE(plplenc_list%enc_parent(nplplenc))
            return
        END SUBROUTINE symba_plplenc_allocate

        SUBROUTINE symba_merger_allocate(mergeadd_list, nmergeadd)
            USE swiftest
            USE module_symba
            IMPLICIT NONE

            ! Arguments
            INTEGER(I4B), INTENT(IN)                :: nmergeadd
            TYPE(symba_merger), INTENT(INOUT)        :: mergeadd_list

            ALLOCATE(mergeadd_list%name(nmergeadd))
            ALLOCATE(mergeadd_list%index_ps(nmergeadd))
            ALLOCATE(mergeadd_list%status(nmergeadd))
            ALLOCATE(mergeadd_list%ncomp(nmergeadd))
            ALLOCATE(mergeadd_list%xh(NDIM,nmergeadd))
            ALLOCATE(mergeadd_list%vh(NDIM,nmergeadd))
            ALLOCATE(mergeadd_list%mass(nmergeadd))
            ALLOCATE(mergeadd_list%radius(nmergeadd))
            return
        END SUBROUTINE symba_merger_allocate

        SUBROUTINE swiftest_tp_allocate(swiftest_tpA, ntp)
            USE swiftest
            USE module_swiftest
            IMPLICIT NONE

            ! Arguments
            INTEGER(I4B), INTENT(IN)            :: ntp
            TYPE(swiftest_tp), INTENT(INOUT)    :: swiftest_tpA

            ALLOCATE(swiftest_tpA%name(ntp))
             ALLOCATE(swiftest_tpA%status(ntp))
             ALLOCATE(swiftest_tpA%peri(ntp))
             ALLOCATE(swiftest_tpA%atp(ntp))
             ALLOCATE(swiftest_tpA%isperi(ntp))
             ALLOCATE(swiftest_tpA%xh(NDIM,ntp))
             ALLOCATE(swiftest_tpA%vh(NDIM,ntp))
             ALLOCATE(swiftest_tpA%xb(NDIM,ntp))
             ALLOCATE(swiftest_tpA%vb(NDIM,ntp))
            return
        END SUBROUTINE swiftest_tp_allocate


        SUBROUTINE helio_tp_allocate(helio_tpA, ntp)
            USE swiftest
            USE module_helio
            IMPLICIT NONE

            ! Arguments
            INTEGER(I4B), INTENT(IN)            :: ntp
            TYPE(helio_tp), INTENT(INOUT)        :: helio_tpA

            ALLOCATE(helio_tpA%ah(NDIM,ntp))
             ALLOCATE(helio_tpA%ahi(NDIM,ntp))
             CALL swiftest_tp_allocate(helio_tpA%swiftest,ntp)

            return
        END SUBROUTINE helio_tp_allocate


        SUBROUTINE symba_tp_allocate(symba_tpA, ntp)
            USE swiftest
            USE module_symba
            IMPLICIT NONE

            ! Arguments
            INTEGER(I4B), INTENT(IN)            :: ntp
            TYPE(symba_tp), INTENT(INOUT)        :: symba_tpA

            ALLOCATE(symba_tpA%nplenc(ntp))
            ALLOCATE(symba_tpA%levelg(ntp))
            ALLOCATE(symba_tpA%levelm(ntp))
            CALL helio_tp_allocate(symba_tpA%helio,ntp)
            return
        END SUBROUTINE symba_tp_allocate

        SUBROUTINE symba_pltpenc_allocate(pltpenc_list, npltpenc)
            USE swiftest
            USE module_symba
            IMPLICIT NONE

            ! Arguments
            INTEGER(I4B), INTENT(IN)                :: npltpenc
            TYPE(symba_pltpenc), INTENT(INOUT)        :: pltpenc_list

            ALLOCATE(pltpenc_list%lvdotr(npltpenc))
            ALLOCATE(pltpenc_list%status(npltpenc))
            ALLOCATE(pltpenc_list%level(npltpenc))
            ALLOCATE(pltpenc_list%indexpl(npltpenc))
            ALLOCATE(pltpenc_list%indextp(npltpenc))
            return
        END SUBROUTINE symba_pltpenc_allocate

!___________________________


        SUBROUTINE swiftest_pl_deallocate(swiftest_plA)
            USE swiftest
            USE module_swiftest
            IMPLICIT NONE

            ! Arguments
            TYPE(swiftest_pl), INTENT(INOUT)    :: swiftest_plA

            DEALLOCATE(swiftest_plA%name)
             DEALLOCATE(swiftest_plA%status)
             DEALLOCATE(swiftest_plA%mass)
             DEALLOCATE(swiftest_plA%radius)
             DEALLOCATE(swiftest_plA%rhill)
             DEALLOCATE(swiftest_plA%xh)
             DEALLOCATE(swiftest_plA%vh)
             DEALLOCATE(swiftest_plA%xb)
             DEALLOCATE(swiftest_plA%vb)
            return
        END SUBROUTINE swiftest_pl_deallocate


        SUBROUTINE helio_pl_deallocate(helio_plA)
            USE swiftest
            USE module_helio
            IMPLICIT NONE

            ! Arguments
            TYPE(helio_pl), INTENT(INOUT)        :: helio_plA

            DEALLOCATE(helio_plA%ah)
             DEALLOCATE(helio_plA%ahi)
             CALL swiftest_pl_deallocate(helio_plA%swiftest)
            return
        END SUBROUTINE helio_pl_deallocate


        SUBROUTINE symba_pl_deallocate(symba_plA)
            USE swiftest
            USE module_symba
            IMPLICIT NONE

            ! Arguments
            TYPE(symba_pl), INTENT(INOUT)        :: symba_plA

            DEALLOCATE(symba_plA%lmerged)
            DEALLOCATE(symba_plA%nplenc)
            DEALLOCATE(symba_plA%ntpenc)
            DEALLOCATE(symba_plA%levelg)
            DEALLOCATE(symba_plA%levelm)
            DEALLOCATE(symba_plA%nchild)
            DEALLOCATE(symba_plA%isperi)
            DEALLOCATE(symba_plA%peri)
            DEALLOCATE(symba_plA%atp)
            DEALLOCATE(symba_plA%index_parent)
            DEALLOCATE(symba_plA%index_child)
            CALL helio_pl_deallocate(symba_plA%helio)
            return
        END SUBROUTINE symba_pl_deallocate

        SUBROUTINE symba_plplenc_deallocate(plplenc_list)
            USE swiftest
            USE module_symba
            IMPLICIT NONE

            ! Arguments
            TYPE(symba_plplenc), INTENT(INOUT)        :: plplenc_list

            DEALLOCATE(plplenc_list%lvdotr)
            DEALLOCATE(plplenc_list%status)
            DEALLOCATE(plplenc_list%level)
            DEALLOCATE(plplenc_list%index1)
            DEALLOCATE(plplenc_list%index2)
            DEALLOCATE(plplenc_list%enc_child)
            DEALLOCATE(plplenc_list%enc_parent)
            return
        END SUBROUTINE symba_plplenc_deallocate

        SUBROUTINE symba_merger_deallocate(mergeadd_list)
            USE swiftest
            USE module_symba
            IMPLICIT NONE

            ! Arguments
            TYPE(symba_merger), INTENT(INOUT)        :: mergeadd_list

            DEALLOCATE(mergeadd_list%name)
            DEALLOCATE(mergeadd_list%index_ps)
            DEALLOCATE(mergeadd_list%status)
            DEALLOCATE(mergeadd_list%ncomp)
            DEALLOCATE(mergeadd_list%xh)
            DEALLOCATE(mergeadd_list%vh)
            DEALLOCATE(mergeadd_list%mass)
            DEALLOCATE(mergeadd_list%radius)
            return
        END SUBROUTINE symba_merger_deallocate

        SUBROUTINE swiftest_tp_deallocate(swiftest_tpA)
            USE swiftest
            USE module_swiftest
            IMPLICIT NONE

            ! Arguments
            TYPE(swiftest_tp), INTENT(INOUT)    :: swiftest_tpA

            DEALLOCATE(swiftest_tpA%name)
             DEALLOCATE(swiftest_tpA%status)
             DEALLOCATE(swiftest_tpA%peri)
             DEALLOCATE(swiftest_tpA%atp)
             DEALLOCATE(swiftest_tpA%isperi)
             DEALLOCATE(swiftest_tpA%xh)
             DEALLOCATE(swiftest_tpA%vh)
             DEALLOCATE(swiftest_tpA%xb)
             DEALLOCATE(swiftest_tpA%vb)
            return
        END SUBROUTINE swiftest_tp_deallocate


        SUBROUTINE helio_tp_deallocate(helio_tpA)
            USE swiftest
            USE module_helio
            IMPLICIT NONE

            ! Arguments
            TYPE(helio_tp), INTENT(INOUT)        :: helio_tpA

            DEALLOCATE(helio_tpA%ah)
             DEALLOCATE(helio_tpA%ahi)
             CALL swiftest_tp_deallocate(helio_tpA%swiftest)
            return
        END SUBROUTINE helio_tp_deallocate


        SUBROUTINE symba_tp_deallocate(symba_tpA)
            USE swiftest
            USE module_symba
            IMPLICIT NONE

            ! Arguments
            TYPE(symba_tp), INTENT(INOUT)        :: symba_tpA

            DEALLOCATE(symba_tpA%nplenc)
            DEALLOCATE(symba_tpA%levelg)
            DEALLOCATE(symba_tpA%levelm)
            
            return
        END SUBROUTINE symba_tp_deallocate

        SUBROUTINE symba_pltpenc_deallocate(pltpenc_list)
            USE swiftest
            USE module_symba
            IMPLICIT NONE

            ! Arguments
            TYPE(symba_pltpenc), INTENT(INOUT)        :: pltpenc_list

            DEALLOCATE(pltpenc_list%lvdotr)
            DEALLOCATE(pltpenc_list%status)
            DEALLOCATE(pltpenc_list%level)
            DEALLOCATE(pltpenc_list%indexpl)
            DEALLOCATE(pltpenc_list%indextp)
            return
        END SUBROUTINE symba_pltpenc_deallocate

        
END MODULE module_swiftestalloc




