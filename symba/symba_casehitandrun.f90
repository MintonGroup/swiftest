!**********************************************************************************************************************************
!
!  Unit Name   : symba_casehitandrun
!  Unit Type   : subroutine
!  Project     : Swiftest
!  Package     : symba
!  Language    : Fortran 90/95
!
!  Description : Merge planets
!
!  Input
!    Arguments : t            : time
!                npl          : number of planets
!                nsppl        : number of spilled planets
!                symba_pl1P   : pointer to head of SyMBA planet structure linked-list
!                symba_pld1P  : pointer to head of discard SyMBA planet structure linked-list
!                nplplenc     : number of planet-planet encounters
!                plplenc_list : array of planet-planet encounter structures
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : npl          : number of planets
!                nsppl        : number of spilled planets
!                symba_pl1P   : pointer to head of SyMBA planet structure linked-list
!                symba_pld1P  : pointer to head of discard SyMBA planet structure linked-list
!    Terminal  : none
!    File      : none
!
!  Invocation  : CALL symba_casehitandrun(t, npl, nsppl, symba_pl1P, symba_pld1P, nplplenc, plplenc_list)
!
!  Notes       : Adapted from Hal Levison's Swift routine discard_mass_merge.f
!
!**********************************************************************************************************************************
SUBROUTINE symba_casehitandrun (t, index_enc, nmergeadd, nmergesub, mergeadd_list, mergesub_list, eoffset, vbs, & 
     npl, symba_plA, nplplenc, plplenc_list, &
     nplmax, ntpmax, fragmax, mres, rres, m1, m2, rad1, rad2, x1, x2, v1, v2)

! Modules
     USE module_parameters
     USE module_swiftest
     USE module_helio
     USE module_symba
     USE module_interfaces, EXCEPT_THIS_ONE => symba_casehitandrun
     IMPLICIT NONE

! Arguments
     INTEGER(I4B), INTENT(IN)                         :: index_enc, nplmax, ntpmax
     INTEGER(I4B), INTENT(INOUT)                      :: npl, nmergeadd, nmergesub, nplplenc, fragmax
     REAL(DP), INTENT(IN)                             :: t
     REAL(DP), INTENT(INOUT)                          :: eoffset, m1, m2, rad1, rad2
     REAL(DP), DIMENSION(3), INTENT(INOUT)            :: mres, rres
     REAL(DP), DIMENSION(NDIM), INTENT(IN)            :: vbs
     REAL(DP), DIMENSION(NDIM), INTENT(INOUT)         :: x1, x2, v1, v2
     TYPE(symba_plplenc), INTENT(INOUT)               :: plplenc_list
     TYPE(symba_merger), INTENT(INOUT)                :: mergeadd_list, mergesub_list
     TYPE(symba_pl), INTENT(INOUT)                    :: symba_plA

! Internals
 
     INTEGER(I4B)                                     :: nfrag, i, k, index1, index2
     INTEGER(I4B)                                     :: index1_parent, index2_parent
     INTEGER(I4B)                                     :: name1, name2, index_keep, index_rm
     REAL(DP)                                         :: mtot, d_rm, m_rm, r_rm, vx_rm, vy_rm
     REAL(DP)                                         :: rhill_keep, r_circle, theta
     REAL(DP)                                         :: m_rem, m_test, mass1, mass2, enew, eold
     REAL(DP)                                         :: x_com, y_com, z_com, vx_com, vy_com, vz_com
     REAL(DP)                                         :: x_frag, y_frag, z_frag, vx_frag, vy_frag, vz_frag
     REAL(DP), DIMENSION(NDIM)                        :: vnew, xr, mv


! Executable code

     nfrag = 4
     index1 = plplenc_list%index1(index_enc)
     index2 = plplenc_list%index2(index_enc)
     index1_parent = symba_plA%index_parent(index1)
     index2_parent = symba_plA%index_parent(index2)
     name1 = symba_plA%helio%swiftest%name(index1)
     name2 = symba_plA%helio%swiftest%name(index2)
     mass1 = symba_plA%helio%swiftest%mass(index1)
     mass2 = symba_plA%helio%swiftest%mass(index2)
     radius1 = symba_plA%helio%swiftest%radius(index1)
     radius2 = symba_plA%helio%swiftest%radius(index2)

     IF (m2 > m1) THEN
          index_keep = index2
          index_rm = index1
     ELSE
          index_keep = index1
          index_rm = index2
     END IF

     ! Find COM
     x_com = ((x1(1) * m1) + (x2(1) * m2)) / (m1 + m2)
     y_com = ((x1(2) * m1) + (x2(2) * m2)) / (m1 + m2)
     z_com = ((x1(3) * m1) + (x2(3) * m2)) / (m1 + m2)

     vx_com = ((v1(1) * m1) + (v2(1) * m2)) / (m1 + m2)
     vy_com = ((v1(2) * m1) + (v2(2) * m2)) / (m1 + m2)
     vz_com = ((v1(3) * m1) + (v2(3) * m2)) / (m1 + m2)

     ! Find energy pre-frag
     eold = 0.5_DP*(m1*DOT_PRODUCT(v1(:), v1(:)) + m2*DOT_PRODUCT(v2(:), v2(:)))
     xr(:) = x2(:) - x1(:)
     eold = eold - (m1*m2/(SQRT(DOT_PRODUCT(xr(:), xr(:)))))

     WRITE(*, *) "Hit and run between particles ", name1, " and ", name2, " at time t = ",t
     WRITE(*, *) "Particle ", index_keep, " survives; Particle ", index_rm, " is fragmented."
     WRITE(*, *) "Number of fragments added: ", (nfrag - 1)

     ! Add both parents to mergesub_list
     nmergesub = nmergesub + 1
     mergesub_list%name(nmergesub) = name1
     mergesub_list%status(nmergesub) = HIT_AND_RUN ! possibly change to disruption for new flag in discard.out
     mergesub_list%xh(:,nmergesub) = x1(:)
     mergesub_list%vh(:,nmergesub) = v1(:) - vbs(:)
     mergesub_list%mass(nmergesub) = mass1
     mergesub_list%radius(nmergesub) = rad1
     nmergesub = nmergesub + 1
     mergesub_list%name(nmergesub) = name2
     mergesub_list%status(nmergesub) = HIT_AND_RUN
     mergesub_list%xh(:,nmergesub) = x2(:)
     mergesub_list%vh(:,nmergesub) = v2(:) - vbs(:)
     mergesub_list%mass(nmergesub) = mass2
     mergesub_list%radius(nmergesub) = rad2

    ! go through the encounter list and for particles actively encoutering
    ! prevent them from having further encounters in this timestep by setting status to MERGED
     DO k = 1, nplplenc 
          IF ((plplenc_list%status(k) == ACTIVE) .AND. &
             ((index1 == plplenc_list%index1(k) .OR. index2 == plplenc_list%index2(k)) .OR. &
             (index2 == plplenc_list%index1(k) .OR. index1 == plplenc_list%index2(k)))) THEN
                    plplenc_list%status(k) = MERGED
          END IF
     END DO

     symba_plA%helio%swiftest%status(index1) = HIT_AND_RUN
     symba_plA%helio%swiftest%status(index2) = HIT_AND_RUN

     ! Calculate the positions of the new fragments
     rhill_keep = symba_plA%helio%swiftest%rhill(index_keep)
     r_circle = rhill_keep
     theta = (2.0_DP * PI) / (nfrag - 1)

     ! Add new fragments to mergeadd_list
     mtot = 0.0_DP ! running total mass of new fragments
     mv = 0.0_DP   ! running sum of m*v of new fragments to be used in COM calculation

     DO i = 1, nfrag
         nmergeadd = nmergeadd + 1
         mergeadd_list%status(nmergeadd) = HIT_AND_RUN
         mergeadd_list%ncomp(nmergeadd) = 2
         m_rm = symba_plA%helio%swiftest%mass(index_rm)
         r_rm = symba_plA%helio%swiftest%radius(index_rm)
         vx_rm = symba_plA%helio%swiftest%vh(1,index_rm)
         vy_rm = symba_plA%helio%swiftest%vh(2,index_rm)
         d_rm = (3.0_DP * m_rm) / (4.0_DP * PI * (r_rm ** 3.0_DP))

         IF (i == 1) THEN
             ! first largest particle equal to index_keep
             mergeadd_list%name(nmergeadd) = symba_plA%helio%swiftest%name(index_keep)
             mergeadd_list%mass(nmergeadd) = symba_plA%helio%swiftest%mass(index_keep)
             mergeadd_list%radius(nmergeadd) = symba_plA%helio%swiftest%radius(index_keep)
             mergeadd_list%xh(:,nmergeadd) = symba_plA%helio%swiftest%xh(:,index_keep)
             mergeadd_list%vh(:,nmergeadd) = symba_plA%helio%swiftest%vh(:,index_keep)
             mtot = mtot + mergeadd_list%mass(nmergeadd)                             
         END IF
         IF (i == 2) THEN
             ! second largest particle from collresolve mres[1] rres[1]
             mergeadd_list%name(nmergeadd) = nplmax + ntpmax + fragmax + i - 1
             mergeadd_list%mass(nmergeadd) = mres(2)
             mergeadd_list%radius(nmergeadd) = rres(2) 
             mtot = mtot + mergeadd_list%mass(nmergeadd)              
         END IF
         IF (i > 2) THEN
             ! FIXME all other particles implement eq. 31 LS12
             ! FIXME current equation taken from Durda et al 2007 Figure 2 Supercatastrophic: N = (1.5e5)e(-1.3*D)
             mergeadd_list%name(nmergeadd) = nplmax + ntpmax + fragmax + i - 1
             m_rem = m_rm - mres(2)
             m_test = (((- 1.0_DP / 2.6_DP) * log(i / (1.5_DP * 10.0_DP ** 5))) ** 3.0_DP) * ((4.0_DP / 3.0_DP) * PI * d_rm)
             
             IF (m_test < m_rem) THEN
                mergeadd_list%mass(nmergeadd) = m_test
             ELSE
                mergeadd_list%mass(nmergeadd) = (m1 + m2) - mtot 
             END IF 
             mergeadd_list%radius(nmergeadd) = ((3.0_DP * mergeadd_list%mass(nmergeadd)) / (4.0_DP * PI * d_rm))  & 
                ** (1.0_DP / 3.0_DP) 
             mtot = mtot + mergeadd_list%mass(nmergeadd)
         END IF
         IF (i > 1) THEN
             x_frag = (r_circle * cos(theta * i)) + x_com
             y_frag = (r_circle * sin(theta * i)) + y_com
             z_frag = z_com
             vx_frag = ((1.0_DP / (nfrag-1)) * (1.0_DP / mergeadd_list%mass(nmergeadd)) * ((m_rm * vx_rm))) - vbs(1)
             vy_frag = ((1.0_DP / (nfrag-1)) * (1.0_DP / mergeadd_list%mass(nmergeadd)) * ((m_rm * vy_rm))) - vbs(2)
             vz_frag = vz_com - vbs(3)
             mergeadd_list%xh(1,nmergeadd) = x_frag
             mergeadd_list%xh(2,nmergeadd) = y_frag 
             mergeadd_list%xh(3,nmergeadd) = z_frag                                                    
             mergeadd_list%vh(1,nmergeadd) = vx_frag
             mergeadd_list%vh(2,nmergeadd) = vy_frag
             mergeadd_list%vh(3,nmergeadd) = vz_frag   
          END IF                                             
          mv = mv + (mergeadd_list%mass(nmergeadd) * mergeadd_list%vh(:,nmergeadd))
     END DO

     ! Calculate energy after frag                                                                           
     vnew(:) = mv / mtot    ! COM of new fragments                               
     enew = 0.5_DP*mtot*DOT_PRODUCT(vnew(:), vnew(:))
     eoffset = eoffset + eold - enew

     ! Update fragmax to account for new fragments
     fragmax = fragmax + (nfrag - 1)

     RETURN 
END SUBROUTINE symba_casehitandrun