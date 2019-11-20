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
SUBROUTINE symba_rearray(t, npl, ntp, nsppl, nsptp, symba_plA, symba_tpA, nmergeadd, mergeadd_list, discard_plA, discard_tpA, discard_plA_id_status,discard_tpA_id_status)

! Modules
     USE module_parameters
     USE module_swiftestalloc 
     USE module_swiftest
     USE module_helio
     USE module_symba
     USE module_interfaces, EXCEPT_THIS_ONE => symba_rearray
     IMPLICIT NONE

! Arguments
     INTEGER(I4B), INTENT(INOUT)                  	:: npl, ntp, nsppl, nsptp, nmergeadd
     REAL(DP), INTENT(IN)                         	:: t
     TYPE(symba_pl), INTENT(INOUT)                	:: symba_plA
     TYPE(symba_tp), INTENT(INOUT)                	:: symba_tpA
     TYPE(symba_merger), DIMENSION(:), INTENT(IN) 	:: mergeadd_list
     REAL(DP), DIMENSION(8,NPLMAX), INTENT(OUT)  	:: discard_plA
     REAL(DP), DIMENSION(8,ntp), INTENT(OUT)     	:: discard_tpA
     INTEGER(I4B), DIMENSION(2,NPLMAX), INTENT(OUT) :: discard_plA_id_status
     INTEGER(I4B), DIMENSION(2,ntp), INTENT(OUT) 	:: discard_tpA_id_status

! Internals
     INTEGER(I4B)                   				:: i, index, j, ncomp, ierr, nplm, nkpl, nktp
     REAL(DP), DIMENSION(8,NPLMAX) 				:: keep_plA
     REAL(DP), DIMENSION(8,ntp)    				:: keep_tpA
     INTEGER(I4B), DIMENSION(2,NPLMAX), INTENT(OUT) :: keep_plA_id_status
     INTEGER(I4B), DIMENSION(2,ntp), INTENT(OUT) 	:: keep_tpA_id_status					


    IF (ldiscard = .TRUE.) THEN 
    	CALL symba_energy(npl, nplmax, symba_plA%helio%swiftest, j2rp2, j4rp4, ke, pe, tei, htot)
    	nsppl = 1
    	nkpl = 1

    	DO i = 1, npl 
        	IF (symba_plA%helio%swiftest%status(i) /= ACTIVE) THEN 
            	discard_plA_id_status(1,nsppl) = symba_plA%helio%swiftest%id(i)
            	discard_plA_id_status(2,nsppl) = symba_plA%helio%swiftest%status(i)
            	discard_plA(1,nsppl) = symba_plA%helio%swiftest%mass(i)
            	discard_plA(2,nsppl) = symba_plA%helio%swiftest%radius(i)
            	discard_plA(3,nsppl) = symba_plA%helio%swiftest%xh(1,i)
            	discard_plA(4,nsppl) = symba_plA%helio%swiftest%xh(2,i)
            	discard_plA(5,nsppl) = symba_plA%helio%swiftest%xh(3,i)
            	discard_plA(6,nsppl) = symba_plA%helio%swiftest%vh(1,i)
            	discard_plA(7,nsppl) = symba_plA%helio%swiftest%vh(2,i)
            	discard_plA(8,nsppl) = symba_plA%helio%swiftest%vh(3,i)
            	nsppl = nsppl+1
        	ELSE
            	keep_plA_id_status(1,nkpl) = symba_plA%helio%swiftest%id(i)
            	keep_plA_id_status(2,nkpl) = symba_plA%helio%swiftest%status(i)
            	keep_plA(1,nkpl) = symba_plA%helio%swiftest%mass(i)
            	keep_plA(2,nkpl) = symba_plA%helio%swiftest%radius(i)
            	keep_plA(3,nkpl) = symba_plA%helio%swiftest%xh(1,i)
            	keep_plA(4,nkpl) = symba_plA%helio%swiftest%xh(2,i)
            	keep_plA(5,nkpl) = symba_plA%helio%swiftest%xh(3,i)
            	keep_plA(6,nkpl) = symba_plA%helio%swiftest%vh(1,i)
            	keep_plA(7,nkpl) = symba_plA%helio%swiftest%vh(2,i)
            	keep_plA(8,nkpl) = symba_plA%helio%swiftest%vh(3,i)
            	nkpl = nkpl+1
        	END IF 
    	END DO

    	IF (nmergeadd = 0) THEN 
    		CALL symba_pl_deallocate(symba_plA,npl)
    		CALL symba_pl_allocate(symba_plA, nkpl)! + nmergeadd)
    		symba_plA%helio%swiftest%id(1:nkpl) = keep_plA_id_status(1,:)
    		symba_plA%helio%swiftest%mass(1:nkpl) = keep_plA(1,:)
    		symba_plA%helio%swiftest%radius(1:nkpl) = keep_plA(2,:)
    		symba_plA%helio%swiftest%xh(1,1:nkpl) = keep_plA(3,:)
    		symba_plA%helio%swiftest%xh(2,1:nkpl) = keep_plA(4,:)
    		symba_plA%helio%swiftest%xh(3,1:nkpl) = keep_plA(5,:)
    		symba_plA%helio%swiftest%vh(1,1:nkpl) = keep_plA(6,:)
    		symba_plA%helio%swiftest%vh(2,1:nkpl) = keep_plA(7,:)
    		symba_plA%helio%swiftest%vh(3,1:nkpl) = keep_plA(8,:)

    	ELSE 
    		CALL symba_pl_deallocate(symba_plA,npl)
    		CALL symba_pl_allocate(symba_plA, nkpl)! + nmergeadd)
    		symba_plA%helio%swiftest%id(1:nkpl) = keep_plA_id_status(1,:)
    		!symba_plA%helio%swiftest%id(nkpl+1:nkpl+nmergeadd) = mergeadd_list%id(:)

    		symba_plA%helio%swiftest%mass(1:nkpl) = keep_plA(1,:)
    		!symba_plA%helio%swiftest%mass(nkpl+1:nkpl+nmergeadd) = mergeadd_list%mass(:)

    		symba_plA%helio%swiftest%radius(1:nkpl) = keep_plA(2,:)
    		!symba_plA%helio%swiftest%radius(nkpl+1:nkpl+nmergeadd) = mergeadd_list%radius(:)

    		symba_plA%helio%swiftest%xh(1,1:nkpl) = keep_plA(3,:)
    		!symba_plA%helio%swiftest%xh(1,nkpl+1:nkpl+nmergeadd) = mergeadd_list%xh(1,:)

    		symba_plA%helio%swiftest%xh(2,1:nkpl) = keep_plA(4,:)
    		!symba_plA%helio%swiftest%xh(2,nkpl+1:nkpl+nmergeadd) = mergeadd_list%xh(2,:)

    		symba_plA%helio%swiftest%xh(3,1:nkpl) = keep_plA(5,:)
    		!symba_plA%helio%swiftest%xh(3,nkpl+1:nkpl+nmergeadd) = mergeadd_list%xh(3,:)

    		symba_plA%helio%swiftest%vh(1,1:nkpl) = keep_plA(6,:)
    		!symba_plA%helio%swiftest%vh(1,nkpl+1:nkpl+nmergeadd) = mergeadd_list%vh(1,:)

    		symba_plA%helio%swiftest%vh(2,1:nkpl) = keep_plA(7,:)
    		!symba_plA%helio%swiftest%vh(2,nkpl+1:nkpl+nmergeadd) = mergeadd_list%vh(2,:)

    		symba_plA%helio%swiftest%vh(3,1:nkpl) = keep_plA(8,:)
    		!symba_plA%helio%swiftest%vh(3,nkpl+1:nkpl+nmergeadd) = mergeadd_list%vh(3,:)
    		CALL symba_energy(npl, nplmax, symba_plA%helio%swiftest, j2rp2, j4rp4, ke, pe, tei, htot
		END IF
	END IF 


	IF (ldiscard_tp = .TRUE.) THEN 
    	nktp = 1
    	nsptp = 1  
    	DO i = 1, ntp
        	IF (symba_tpA%helio%swiftest%status(i) /= ACTIVE) THEN 
            	discard_tpA_id_status(1,nsptp) = symba_tpA%helio%swiftest%id(i)
            	discard_tpA_id_status(2,nsptp) = symba_tpA%helio%swiftest%status(i)
            	discard_tpA(1,nsptp) = symba_tpA%helio%swiftest%mass(i)
            	discard_tpA(2,nsptp) = symba_tpA%helio%swiftest%radius(i)
            	discard_tpA(3,nsptp) = symba_tpA%helio%swiftest%xh(1,i)
            	discard_tpA(4,nsptp) = symba_tpA%helio%swiftest%xh(2,i)
            	discard_tpA(5,nsptp) = symba_tpA%helio%swiftest%xh(3,i)
            	discard_tpA(6,nsptp) = symba_tpA%helio%swiftest%vh(1,i)
            	discard_tpA(7,nsptp) = symba_tpA%helio%swiftest%vh(2,i)
            	discard_tpA(8,nsptp) = symba_tpA%helio%swiftest%vh(3,i)
            	nsptp = nsptp+1
        	ELSE 
            	keep_tpA_id_status(1,nktp) = symba_tpA%helio%swiftest%id(i)
            	keep_tpA_id_status(2,nktp) = symba_tpA%helio%swiftest%status(i)
            	keep_tpA(1,nktp) = symba_tpA%helio%swiftest%mass(i)
            	keep_tpA(2,nktp) = symba_tpA%helio%swiftest%radius(i)
            	keep_tpA(3,nktp) = symba_tpA%helio%swiftest%xh(1,i)
            	keep_tpA(4,nktp) = symba_tpA%helio%swiftest%xh(2,i)
            	keep_tpA(5,nktp) = symba_tpA%helio%swiftest%xh(3,i)
            	keep_tpA(6,nktp) = symba_tpA%helio%swiftest%vh(1,i)
            	keep_tpA(7,nktp) = symba_tpA%helio%swiftest%vh(2,i)
            	keep_tpA(8,nktp) = symba_tpA%helio%swiftest%vh(3,i)
            	nktp = nktp+1
        	END IF 
    	END DO 

    	CALL symba_tp_deallocate(symba_tpA,ntp)
    	CALL symba_tp_allocate(symba_tpA, nktp)


    	symba_tpA%helio%swiftest%id(1:nktp) = keep_tpA_id_status(1,:)
    	symba_tpA%helio%swiftest%mass(1:nktp) = keep_tpA(1,:)
    	symba_tpA%helio%swiftest%radius(1:nktp) = keep_tpA(2,:)
    	symba_tpA%helio%swiftest%xh(1,1:nktp) = keep_tpA(3,:)
    	symba_tpA%helio%swiftest%xh(2,1:nktp) = keep_tpA(4,:)
    	symba_tpA%helio%swiftest%xh(3,1:nktp) = keep_tpA(5,:)
    	symba_tpA%helio%swiftest%vh(1,1:nktp) = keep_tpA(6,:)
    	symba_tpA%helio%swiftest%vh(2,1:nktp) = keep_tpA(7,:)
    	symba_tpA%helio%swiftest%vh(3,1:nktp) = keep_tpA(8,:)

    END IF 


END SUBROUTINE symba_rearray


    


