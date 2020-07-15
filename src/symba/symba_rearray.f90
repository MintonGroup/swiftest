submodule (symba) s_symba_rearray
contains
   module procedure symba_rearray
   !! author: Jennifer L. L. Pouplin, Carlisle A. Wishard, and David A. Minton
   !!
   !! Redo array of pl and tp based on discarded and added pl and tp
   use swiftest
   implicit none
   integer(I4B)                     :: i, nkpl, nktp, nfrag
   real(DP)                         :: mu, energy, ap, r, v2
   logical, dimension(npl)                :: discard_l_pl, frag_l_add
   logical, dimension(ntp)                :: discard_l_tp

! executable code

    if (ldiscard) then 
      nsppl = 0
      nkpl = 0
      discard_l_pl(1:npl) = (symba_plA%status(1:npl) /= ACTIVE) 
      nsppl = count(discard_l_pl)
      nkpl = npl - nsppl
      frag_l_add = [(.false.,i=1,npl)]
      if (config%lfragmentation) then
        do i = 1, npl
          if (mergeadd_list%status(i) == DISRUPTION) then
            frag_l_add(i) = .true.
          else if (mergeadd_list%status(i) == HIT_AND_RUN) then
            frag_l_add(i) = .true.
          else if (mergeadd_list%status(i) == SUPERCATASTROPHIC) then
            frag_l_add(i) = .true.
          else
            frag_l_add(i) = .false.
          end if
        end do
      end if
      nfrag = count(frag_l_add)

      call discard_plA%alloc(nsppl)

      discard_plA%name(1:nsppl) = pack(symba_plA%name(1:npl), discard_l_pl)
      discard_plA%status(1:nsppl) = pack(symba_plA%status(1:npl), discard_l_pl)
      discard_plA%mass(1:nsppl) = pack(symba_plA%mass(1:npl), discard_l_pl)
      discard_plA%radius(1:nsppl) = pack(symba_plA%radius(1:npl), discard_l_pl)
      discard_plA%xh(1:nsppl, 1) = pack(symba_plA%xh(1:npl, 1), discard_l_pl)
      discard_plA%xh(1:nsppl, 2) = pack(symba_plA%xh(1:npl, 2), discard_l_pl)
      discard_plA%xh(1:nsppl, 3) = pack(symba_plA%xh(1:npl, 3), discard_l_pl)
      discard_plA%vh(1:nsppl, 1) = pack(symba_plA%vh(1:npl, 1), discard_l_pl)
      discard_plA%vh(1:nsppl, 2) = pack(symba_plA%vh(1:npl, 2), discard_l_pl)
      discard_plA%vh(1:nsppl, 3) = pack(symba_plA%vh(1:npl, 3), discard_l_pl)
      discard_plA%rhill(1:nsppl) = pack(symba_plA%rhill(1:npl), discard_l_pl)
      discard_plA%xb(1:nsppl, 1) = pack(symba_plA%xb(1:npl, 1), discard_l_pl)
      discard_plA%xb(1:nsppl, 2) = pack(symba_plA%xb(1:npl, 2), discard_l_pl)
      discard_plA%xb(1:nsppl, 3) = pack(symba_plA%xb(1:npl, 3), discard_l_pl)
      discard_plA%vb(1:nsppl, 1) = pack(symba_plA%vb(1:npl, 1), discard_l_pl)
      discard_plA%vb(1:nsppl, 2) = pack(symba_plA%vb(1:npl, 2), discard_l_pl)
      discard_plA%vb(1:nsppl, 3) = pack(symba_plA%vb(1:npl, 3), discard_l_pl)
      if (config%lfragmentation .and. (nkpl + nfrag > npl)) then 
        symba_plA%name(1:nkpl) = pack(symba_plA%name(1:npl), .not. discard_l_pl)
        symba_plA%status(1:nkpl) = pack(symba_plA%status(1:npl), .not. discard_l_pl)
        symba_plA%mass(1:nkpl) = pack(symba_plA%mass(1:npl), .not. discard_l_pl)
        symba_plA%radius(1:nkpl) = pack(symba_plA%radius(1:npl), .not. discard_l_pl)
        symba_plA%xh(1:nkpl, 1) = pack(symba_plA%xh(1:npl, 1), .not. discard_l_pl)
        symba_plA%xh(1:nkpl, 2) = pack(symba_plA%xh(1:npl, 2), .not. discard_l_pl)
        symba_plA%xh(1:nkpl, 3) = pack(symba_plA%xh(1:npl, 3), .not. discard_l_pl)
        symba_plA%vh(1:nkpl, 1) = pack(symba_plA%vh(1:npl, 1), .not. discard_l_pl)
        symba_plA%vh(1:nkpl, 2) = pack(symba_plA%vh(1:npl, 2), .not. discard_l_pl)
        symba_plA%vh(1:nkpl, 3) = pack(symba_plA%vh(1:npl, 3), .not. discard_l_pl)
        symba_plA%rhill(1:nkpl) = pack(symba_plA%rhill(1:npl), .not. discard_l_pl)
        symba_plA%xb(1:nkpl, 1) = pack(symba_plA%xb(1:npl, 1), .not. discard_l_pl)
        symba_plA%xb(1:nkpl, 2) = pack(symba_plA%xb(1:npl, 2), .not. discard_l_pl)
        symba_plA%xb(1:nkpl, 3) = pack(symba_plA%xb(1:npl, 3), .not. discard_l_pl)
        symba_plA%vb(1:nkpl, 1) = pack(symba_plA%vb(1:npl, 1), .not. discard_l_pl)
        symba_plA%vb(1:nkpl, 2) = pack(symba_plA%vb(1:npl, 2), .not. discard_l_pl)
        symba_plA%vb(1:nkpl, 3) = pack(symba_plA%vb(1:npl, 3), .not. discard_l_pl)
        symba_plA%ah(1:nkpl, 1) = pack(symba_plA%ah(1:npl, 1), .not. discard_l_pl)
        symba_plA%ah(1:nkpl, 2) = pack(symba_plA%ah(1:npl, 2), .not. discard_l_pl)
        symba_plA%ah(1:nkpl, 3) =pack(symba_plA%ah(1:npl, 3), .not. discard_l_pl)

        call util_resize_pl(symba_plA, nkpl+nfrag, npl)

        npl = nkpl  + nfrag
        !add fragments 
        symba_plA%name(nkpl+1:npl) = pack(mergeadd_list%name(1:nmergeadd), frag_l_add)
        symba_plA%status(nkpl+1:npl) = [(ACTIVE,i=1,nfrag)]!array of ACTIVE status 
        symba_plA%mass(nkpl+1:npl) = pack(mergeadd_list%mass(1:nmergeadd), frag_l_add)
        symba_plA%radius(nkpl+1:npl) = pack(mergeadd_list%radius(1:nmergeadd), frag_l_add)
        symba_plA%xh(1,nkpl+1:npl) = pack(mergeadd_list%xh(1:n, 1mergeadd), frag_l_add)
        symba_plA%xh(2,nkpl+1:npl) = pack(mergeadd_list%xh(1:n, 2mergeadd), frag_l_add)
        symba_plA%xh(3,nkpl+1:npl) = pack(mergeadd_list%xh(1:n, 3mergeadd), frag_l_add)
        symba_plA%vh(1,nkpl+1:npl) = pack(mergeadd_list%vh(1:n, 1mergeadd), frag_l_add)
        symba_plA%vh(2,nkpl+1:npl) = pack(mergeadd_list%vh(1:n, 2mergeadd), frag_l_add)
        symba_plA%vh(3,nkpl+1:npl) = pack(mergeadd_list%vh(1:n, 3mergeadd), frag_l_add)

        do i = nkpl+1, npl
          mu = symba_plA%mass(1) + symba_plA%mass(i)
          r = norm2(symba_plA%xh(:,i))
          v2 = dot_product(symba_plA%vh(:,i), symba_plA%vh(:,i))
          energy = 0.5_DP*v2 - mu/r
          ap = -0.5_DP*mu/energy
         symba_plA%rhill(i) = ap*(((symba_plA%mass(i)/mu)/3.0_DP)**(1.0_DP/3.0_DP))
        end do

      else
        symba_plA%name(1:nkpl) = pack(symba_plA%name(1:npl), .not. discard_l_pl)
        symba_plA%status(1:nkpl) = pack(symba_plA%status(1:npl), .not. discard_l_pl)
        symba_plA%mass(1:nkpl) = pack(symba_plA%mass(1:npl), .not. discard_l_pl)
        symba_plA%radius(1:nkpl) = pack(symba_plA%radius(1:npl), .not. discard_l_pl)
        symba_plA%xh(1:nkpl, 1) = pack(symba_plA%xh(1:npl, 1), .not. discard_l_pl)
        symba_plA%xh(1:nkpl, 2) = pack(symba_plA%xh(1:npl, 2), .not. discard_l_pl)
        symba_plA%xh(1:nkpl, 3) = pack(symba_plA%xh(1:npl, 3), .not. discard_l_pl)
        symba_plA%vh(1:nkpl, 1) = pack(symba_plA%vh(1:npl, 1), .not. discard_l_pl)
        symba_plA%vh(1:nkpl, 2) = pack(symba_plA%vh(1:npl, 2), .not. discard_l_pl)
        symba_plA%vh(1:nkpl, 3) = pack(symba_plA%vh(1:npl, 3), .not. discard_l_pl)
        symba_plA%rhill(1:nkpl) = pack(symba_plA%rhill(1:npl), .not. discard_l_pl)
        symba_plA%xb(1:nkpl, 1) = pack(symba_plA%xb(1:npl, 1), .not. discard_l_pl)
        symba_plA%xb(1:nkpl, 2) = pack(symba_plA%xb(1:npl, 2), .not. discard_l_pl)
        symba_plA%xb(1:nkpl, 3) = pack(symba_plA%xb(1:npl, 3), .not. discard_l_pl)
        symba_plA%vb(1:nkpl, 1) = pack(symba_plA%vb(1:npl, 1), .not. discard_l_pl)
        symba_plA%vb(1:nkpl, 2) = pack(symba_plA%vb(1:npl, 2), .not. discard_l_pl)
        symba_plA%vb(1:nkpl, 3) = pack(symba_plA%vb(1:npl, 3), .not. discard_l_pl)
        symba_plA%ah(1:nkpl, 1) = pack(symba_plA%ah(1:npl, 1), .not. discard_l_pl)
        symba_plA%ah(1:nkpl, 2) = pack(symba_plA%ah(1:npl, 2), .not. discard_l_pl)
        symba_plA%ah(1:nkpl, 3) = pack(symba_plA%ah(1:npl, 3), .not. discard_l_pl)
        npl = nkpl
        symba_plA%nbody = npl
      end if
    end if 

    if (ldiscard_tp) then 
      nktp = 0
      nsptp = 0  

      discard_l_tp(1:ntp) = (symba_tpA%status(1:ntp) /= ACTIVE)
      nsptp = count(discard_l_tp)
      nktp = ntp - nsptp

      call discard_tpA%alloc(nsptp) 

      discard_tpA%name(1:nsptp) = pack(symba_tpA%name(1:ntp), discard_l_tp)
      discard_tpA%status(1:nsptp) = pack(symba_tpA%status(1:ntp), discard_l_tp)
      discard_tpA%xh(1:nsptp, 1) = pack(symba_tpA%xh(1:ntp, 1), discard_l_tp)
      discard_tpA%xh(1:nsptp, 2) = pack(symba_tpA%xh(1:ntp, 2), discard_l_tp)
      discard_tpA%xh(1:nsptp, 3) = pack(symba_tpA%xh(1:ntp, 3), discard_l_tp)
      discard_tpA%vh(1:nsptp, 1) = pack(symba_tpA%vh(1:ntp, 1), discard_l_tp)
      discard_tpA%vh(1:nsptp, 2) = pack(symba_tpA%vh(1:ntp, 2), discard_l_tp)
      discard_tpA%vh(1:nsptp, 3) = pack(symba_tpA%vh(1:ntp, 3), discard_l_tp)
      discard_tpA%isperi(1:nsptp) = pack(symba_tpA%isperi(1:ntp), discard_l_tp)
      discard_tpA%peri(1:nsptp) = pack(symba_tpA%peri(1:ntp), discard_l_tp)
      discard_tpA%atp(1:nsptp) = pack(symba_tpA%atp(1:ntp), discard_l_tp)
      discard_tpA%xb(1:nsptp, 1) = pack(symba_tpA%xb(1:ntp, 1), discard_l_tp)
      discard_tpA%xb(1:nsptp, 2) = pack(symba_tpA%xb(1:ntp, 2), discard_l_tp)
      discard_tpA%xb(1:nsptp, 3) = pack(symba_tpA%xb(1:ntp, 3), discard_l_tp)
      discard_tpA%vb(1:nsptp, 1) = pack(symba_tpA%vb(1:ntp, 1), discard_l_tp)
      discard_tpA%vb(1:nsptp, 2) = pack(symba_tpA%vb(1:ntp, 2), discard_l_tp)
      discard_tpA%vb(1:nsptp, 3) = pack(symba_tpA%vb(1:ntp, 3), discard_l_tp)

      symba_tpA%name(1:nktp) = pack(symba_tpA%name(1:ntp), .not. discard_l_tp)
      symba_tpA%status(1:nktp) = pack(symba_tpA%status(1:ntp), .not. discard_l_tp)
      symba_tpA%xh(1:nktp, 1) = pack(symba_tpA%xh(1:ntp, 1), .not. discard_l_tp)
      symba_tpA%xh(1:nktp, 2) = pack(symba_tpA%xh(1:ntp, 2), .not. discard_l_tp)
      symba_tpA%xh(1:nktp, 3) = pack(symba_tpA%xh(1:ntp, 3), .not. discard_l_tp)
      symba_tpA%vh(1:nktp, 1) = pack(symba_tpA%vh(1:ntp, 1), .not. discard_l_tp)
      symba_tpA%vh(1:nktp, 2) = pack(symba_tpA%vh(1:ntp, 2), .not. discard_l_tp)
      symba_tpA%vh(1:nktp, 3) = pack(symba_tpA%vh(1:ntp, 3), .not. discard_l_tp)
      symba_tpA%xb(1:nktp, 1) = pack(symba_tpA%xb(1:ntp, 1), .not. discard_l_tp)
      symba_tpA%xb(1:nktp, 2) = pack(symba_tpA%xb(1:ntp, 2), .not. discard_l_tp)
      symba_tpA%xb(1:nktp, 3) = pack(symba_tpA%xb(1:ntp, 3), .not. discard_l_tp)
      symba_tpA%vb(1:nktp, 1) = pack(symba_tpA%vb(1:ntp, 1), .not. discard_l_tp)
      symba_tpA%vb(1:nktp, 2) = pack(symba_tpA%vb(1:ntp, 2), .not. discard_l_tp)
      symba_tpA%vb(1:nktp, 3) = pack(symba_tpA%vb(1:ntp, 3), .not. discard_l_tp)
      symba_tpA%isperi(1:nktp) = pack(symba_tpA%isperi(1:ntp), .not. discard_l_tp)
      symba_tpA%peri(1:nktp) = pack(symba_tpA%peri(1:ntp), .not. discard_l_tp)
      symba_tpA%atp(1:nktp) = pack(symba_tpA%atp(1:ntp), .not. discard_l_tp)
      symba_tpA%ah(1:nktp, 1) = pack(symba_tpA%ah(1:ntp, 1), .not. discard_l_tp)
      symba_tpA%ah(1:nktp, 2) = pack(symba_tpA%ah(1:ntp, 2), .not. discard_l_tp)
      symba_tpA%ah(1:nktp, 3) = pack(symba_tpA%ah(1:ntp, 3), .not. discard_l_tp)
      ntp = nktp
      symba_tpA%nbody = ntp
    end if 

  end procedure symba_rearray
end submodule s_symba_rearray
