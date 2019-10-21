!**********************************************************************************************************************************
!
!  Unit Name   : module_swiftest_allocation
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
MODULE module_swiftest_allocation

	USE module_parameters
	IMPLICIT NONE

	CONTAINS 

		SUBROUTINE swiftest_pl_allocate(swiftest_plA, npl)
			USE module_parameters
			USE module_swiftest
			IMPLICIT NONE

			! Arguments
			INTEGER(I4B), INTENT(IN)			:: npl
			TYPE(swiftest_pl), INTENT(INOUT)	:: swiftest_plA

			ALLOCATE(swiftest_plA%id(npl))
         	ALLOCATE(swiftest_plA%status(npl))
         	ALLOCATE(swiftest_plA%mass(npl))
         	ALLOCATE(swiftest_plA%radius(npl))
         	ALLOCATE(swiftest_plA%rhill(npl))
         	ALLOCATE(swiftest_plA%xh(NDIM,npl))
         	ALLOCATE(swiftest_plA%vh(NDIM,npl))
         	ALLOCATE(swiftest_plA%xb(NDIM,npl))
         	ALLOCATE(swiftest_plA%vb(NDIM,npl))
        	return
        END SUBROUTINE


        SUBROUTINE helio_pl_allocate(helio_plA, npl)
			USE module_parameters
			USE module_helio
			IMPLICIT NONE

			! Arguments
			INTEGER(I4B), INTENT(IN)			:: npl
			TYPE(helio_pl), INTENT(INOUT)		:: helio_plA

			ALLOCATE(helio_plA%ah(npl))
         	ALLOCATE(helio_plA%ahi(npl))
         	ALLOCATE(helio_plA%swiftest(npl))
        	return
        END SUBROUTINE


        SUBROUTINE symba_pl_allocate(symba_plA, npl)
			USE module_parameters
			USE module_symba
			IMPLICIT NONE

			! Arguments
			INTEGER(I4B), INTENT(IN)			:: npl
			TYPE(symba_pl), INTENT(INOUT)		:: symba_plA

			ALLOCATE(symba_plA%lmerged(npl))
			ALLOCATE(symba_plA%nplenc(npl))
			ALLOCATE(symba_plA%ntpenc(npl))
			ALLOCATE(symba_plA%levelg(npl))
			ALLOCATE(symba_plA%levelm(npl))
			ALLOCATE(symba_plA%nchild(npl))
			ALLOCATE(symba_plA%isperi(npl))
			ALLOCATE(symba_plA%peri(npl))
			ALLOCATE(symba_plA%atp(npl))
			ALLOCATE(symba_plA%helio(npl))
			return
		END SUBROUTINE

		SUBROUTINE symba_plplenc_allocate(plplenc_list, nplplenc)
			USE module_parameters
			USE module_symba
			IMPLICIT NONE

			! Arguments
			INTEGER(I4B), INTENT(IN)				:: nplplenc
			TYPE(symba_plplenc), INTENT(INOUT)		:: plplenc_list

			ALLOCATE(plplenc_list%lvdotr(nplplenc))
			ALLOCATE(plplenc_list%status(nplplenc))
			ALLOCATE(plplenc_list%level(nplplenc))
			return
		END SUBROUTINE

		SUBROUTINE symba_merger_allocate(mergeadd_list, nmergeadd)
			USE module_parameters
			USE module_symba
			IMPLICIT NONE

			! Arguments
			INTEGER(I4B), INTENT(IN)				:: nmergeadd
			TYPE(symba_plplenc), INTENT(INOUT)		:: mergeadd_list

			ALLOCATE(mergeadd_list%id(nmergeadd))
			ALLOCATE(mergeadd_list%status(nmergeadd))
			ALLOCATE(mergeadd_list%ncomp(nmergeadd))
			ALLOCATE(mergeadd_list%xh(nmergeadd))
			ALLOCATE(mergeadd_list%vh(nmergeadd))
			ALLOCATE(mergeadd_list%mass(nmergeadd))
			ALLOCATE(mergeadd_list%radius(nmergeadd))
			return
		END SUBROUTINE

		SUBROUTINE swiftest_tp_allocate(swiftest_tpA, ntp)
			USE module_parameters
			USE module_swiftest
			IMPLICIT NONE

			! Arguments
			INTEGER(I4B), INTENT(IN)			:: ntp
			TYPE(swiftest_tp), INTENT(INOUT)	:: swiftest_tpA

			ALLOCATE(swiftest_tpA%id(ntp))
         	ALLOCATE(swiftest_tpA%status(ntp))
         	ALLOCATE(swiftest_tpA%mass(ntp))
         	ALLOCATE(swiftest_tpA%radius(ntp))
         	ALLOCATE(swiftest_tpA%rhill(ntp))
         	ALLOCATE(swiftest_tpA%xh(NDIM,ntp))
         	ALLOCATE(swiftest_tpA%vh(NDIM,ntp))
         	ALLOCATE(swiftest_tpA%xb(NDIM,ntp))
         	ALLOCATE(swiftest_tpA%vb(NDIM,ntp))
        	return
        END SUBROUTINE


        SUBROUTINE helio_tp_allocate(helio_tpA, ntp)
			USE module_parameters
			USE module_helio
			IMPLICIT NONE

			! Arguments
			INTEGER(I4B), INTENT(IN)			:: ntp
			TYPE(helio_tp), INTENT(INOUT)		:: helio_tpA

			ALLOCATE(helio_tpA%ah(ntp))
         	ALLOCATE(helio_tpA%ahi(ntp))
         	ALLOCATE(helio_tpA%swiftest(ntp))
        	return
        END SUBROUTINE


        SUBROUTINE symba_tp_allocate(symba_tpA, ntp)
			USE module_parameters
			USE module_symba
			IMPLICIT NONE

			! Arguments
			INTEGER(I4B), INTENT(IN)			:: ntp
			TYPE(symba_tp), INTENT(INOUT)		:: symba_tpA

			ALLOCATE(symba_tpA%lmerged(ntp))
			ALLOCATE(symba_tpA%nplenc(ntp))
			ALLOCATE(symba_tpA%ntpenc(ntp))
			ALLOCATE(symba_tpA%levelg(ntp))
			ALLOCATE(symba_tpA%levelm(ntp))
			ALLOCATE(symba_tpA%nchild(ntp))
			ALLOCATE(symba_tpA%isperi(ntp))
			ALLOCATE(symba_tpA%peri(ntp))
			ALLOCATE(symba_tpA%atp(ntp))
			ALLOCATE(symba_tpA%helio(ntp))
			return
		END SUBROUTINE

		SUBROUTINE symba_pltpenc_allocate(pltpenc_list, npltpenc)
			USE module_parameters
			USE module_symba
			IMPLICIT NONE

			! Arguments
			INTEGER(I4B), INTENT(IN)				:: npltpenc
			TYPE(symba_pltpenc), INTENT(INOUT)		:: pltpenc_list

			ALLOCATE(pltpenc_list%lvdotr(npltpenc))
			ALLOCATE(pltpenc_list%status(npltpenc))
			ALLOCATE(pltpenc_list%level(npltpenc))
			return
		END SUBROUTINE

!___________________________

		SUBROUTINE swiftest_pl_deallocate(swiftest_plA, npl)
			USE module_parameters
			USE module_swiftest
			IMPLICIT NONE

			! Arguments
			INTEGER(I4B), INTENT(IN)			:: npl
			TYPE(swiftest_pl), INTENT(INOUT)	:: swiftest_plA

			DEALLOCATE(swiftest_plA%id(npl))
         	DEALLOCATE(swiftest_plA%status(npl))
         	DEALLOCATE(swiftest_plA%mass(npl))
         	DEALLOCATE(swiftest_plA%radius(npl))
         	DEALLOCATE(swiftest_plA%rhill(npl))
         	DEALLOCATE(swiftest_plA%xh(NDIM,npl))
         	DEALLOCATE(swiftest_plA%vh(NDIM,npl))
         	DEALLOCATE(swiftest_plA%xb(NDIM,npl))
         	DEALLOCATE(swiftest_plA%vb(NDIM,npl))
        	return
        END SUBROUTINE


        SUBROUTINE helio_pl_deallocate(helio_plA, npl)
			USE module_parameters
			USE module_helio
			IMPLICIT NONE

			! Arguments
			INTEGER(I4B), INTENT(IN)			:: npl
			TYPE(helio_pl), INTENT(INOUT)		:: helio_plA

			DEALLOCATE(helio_plA%ah(npl))
         	DEALLOCATE(helio_plA%ahi(npl))
         	DEALLOCATE(helio_plA%swiftest(npl))
        	return
        END SUBROUTINE


        SUBROUTINE symba_pl_deallocate(symba_plA, npl)
			USE module_parameters
			USE module_symba
			IMPLICIT NONE

			! Arguments
			INTEGER(I4B), INTENT(IN)			:: npl
			TYPE(symba_pl), INTENT(INOUT)		:: symba_plA

			DEALLOCATE(symba_plA%lmerged(npl))
			DEALLOCATE(symba_plA%nplenc(npl))
			DEALLOCATE(symba_plA%ntpenc(npl))
			DEALLOCATE(symba_plA%levelg(npl))
			DEALLOCATE(symba_plA%levelm(npl))
			DEALLOCATE(symba_plA%nchild(npl))
			DEALLOCATE(symba_plA%isperi(npl))
			DEALLOCATE(symba_plA%peri(npl))
			DEALLOCATE(symba_plA%atp(npl))
			DEALLOCATE(symba_plA%helio(npl))
			return
		END SUBROUTINE

		SUBROUTINE symba_plplenc_deallocate(plplenc_list, nplplenc)
			USE module_parameters
			USE module_symba
			IMPLICIT NONE

			! Arguments
			INTEGER(I4B), INTENT(IN)				:: nplplenc
			TYPE(symba_plplenc), INTENT(INOUT)		:: plplenc_list

			DEALLOCATE(plplenc_list%lvdotr(nplplenc))
			DEALLOCATE(plplenc_list%status(nplplenc))
			DEALLOCATE(plplenc_list%level(nplplenc))
			return
		END SUBROUTINE

		SUBROUTINE symba_merger_deallocate(mergeadd_list, nmergeadd)
			USE module_parameters
			USE module_symba
			IMPLICIT NONE

			! Arguments
			INTEGER(I4B), INTENT(IN)				:: nmergeadd
			TYPE(symba_plplenc), INTENT(INOUT)		:: mergeadd_list

			DEALLOCATE(mergeadd_list%id(nmergeadd))
			DEALLOCATE(mergeadd_list%status(nmergeadd))
			DEALLOCATE(mergeadd_list%ncomp(nmergeadd))
			DEALLOCATE(mergeadd_list%xh(nmergeadd))
			DEALLOCATE(mergeadd_list%vh(nmergeadd))
			DEALLOCATE(mergeadd_list%mass(nmergeadd))
			DEALLOCATE(mergeadd_list%radius(nmergeadd))
			return
		END SUBROUTINE

		SUBROUTINE swiftest_tp_deallocate(swiftest_tpA, ntp)
			USE module_parameters
			USE module_swiftest
			IMPLICIT NONE

			! Arguments
			INTEGER(I4B), INTENT(IN)			:: ntp
			TYPE(swiftest_tp), INTENT(INOUT)	:: swiftest_tpA

			DEALLOCATE(swiftest_tpA%id(ntp))
         	DEALLOCATE(swiftest_tpA%status(ntp))
         	DEALLOCATE(swiftest_tpA%mass(ntp))
         	DEALLOCATE(swiftest_tpA%radius(ntp))
         	DEALLOCATE(swiftest_tpA%rhill(ntp))
         	DEALLOCATE(swiftest_tpA%xh(NDIM,ntp))
         	DEALLOCATE(swiftest_tpA%vh(NDIM,ntp))
         	DEALLOCATE(swiftest_tpA%xb(NDIM,ntp))
         	DEALLOCATE(swiftest_tpA%vb(NDIM,ntp))
        	return
        END SUBROUTINE


        SUBROUTINE helio_tp_deallocate(helio_tpA, ntp)
			USE module_parameters
			USE module_helio
			IMPLICIT NONE

			! Arguments
			INTEGER(I4B), INTENT(IN)			:: ntp
			TYPE(helio_tp), INTENT(INOUT)		:: helio_tpA

			DEALLOCATE(helio_tpA%ah(ntp))
         	DEALLOCATE(helio_tpA%ahi(ntp))
         	DEALLOCATE(helio_tpA%swiftest(ntp))
        	return
        END SUBROUTINE


        SUBROUTINE symba_tp_deallocate(symba_tpA, ntp)
			USE module_parameters
			USE module_symba
			IMPLICIT NONE

			! Arguments
			INTEGER(I4B), INTENT(IN)			:: ntp
			TYPE(symba_tp), INTENT(INOUT)		:: symba_tpA

			DEALLOCATE(symba_tpA%lmerged(ntp))
			DEALLOCATE(symba_tpA%nplenc(ntp))
			DEALLOCATE(symba_tpA%ntpenc(ntp))
			DEALLOCATE(symba_tpA%levelg(ntp))
			DEALLOCATE(symba_tpA%levelm(ntp))
			DEALLOCATE(symba_tpA%nchild(ntp))
			DEALLOCATE(symba_tpA%isperi(ntp))
			DEALLOCATE(symba_tpA%peri(ntp))
			DEALLOCATE(symba_tpA%atp(ntp))
			DEALLOCATE(symba_tpA%helio(ntp))
			return
		END SUBROUTINE

		SUBROUTINE symba_pltpenc_deallocate(pltpenc_list, npltpenc)
			USE module_parameters
			USE module_symba
			IMPLICIT NONE

			! Arguments
			INTEGER(I4B), INTENT(IN)				:: npltpenc
			TYPE(symba_pltpenc), INTENT(INOUT)		:: pltpenc_list

			DEALLOCATE(pltpenc_list%lvdotr(npltpenc))
			DEALLOCATE(pltpenc_list%status(npltpenc))
			DEALLOCATE(pltpenc_list%level(npltpenc))
			return
		END SUBROUTINE
		
END PROGRAM module_swiftest_allocation




