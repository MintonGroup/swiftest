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
      discard_l_pl(1:npl) = (symba_pla%helio%swiftest%status(1:npl) /= active) 
      nsppl = count(discard_l_pl)
      nkpl = npl - nsppl
      frag_l_add = [(.false.,i=1,npl)]
      if (config%lfragmentation) then
        do i = 1, npl
          if (mergeadd_list%status(i) == disruption) then
            frag_l_add(i) = .true.
          else if (mergeadd_list%status(i) == hit_and_run) then
            frag_l_add(i) = .true.
          else if (mergeadd_list%status(i) == supercatastrophic) then
            frag_l_add(i) = .true.
          else
            frag_l_add(i) = .false.
          end if
        end do
      end if
      nfrag = count(frag_l_add)

      call discard_pla%alloc(nsppl)

      discard_pla%name(1:nsppl) = pack(symba_pla%helio%swiftest%name(1:npl), discard_l_pl)
      discard_pla%status(1:nsppl) = pack(symba_pla%helio%swiftest%status(1:npl), discard_l_pl)
      discard_pla%mass(1:nsppl) = pack(symba_pla%helio%swiftest%mass(1:npl), discard_l_pl)
      discard_pla%radius(1:nsppl) = pack(symba_pla%helio%swiftest%radius(1:npl), discard_l_pl)
      discard_pla%xh(1,1:nsppl) = pack(symba_pla%helio%swiftest%xh(1,1:npl), discard_l_pl)
      discard_pla%xh(2,1:nsppl) = pack(symba_pla%helio%swiftest%xh(2,1:npl), discard_l_pl)
      discard_pla%xh(3,1:nsppl) = pack(symba_pla%helio%swiftest%xh(3,1:npl), discard_l_pl)
      discard_pla%vh(1,1:nsppl) = pack(symba_pla%helio%swiftest%vh(1,1:npl), discard_l_pl)
      discard_pla%vh(2,1:nsppl) = pack(symba_pla%helio%swiftest%vh(2,1:npl), discard_l_pl)
      discard_pla%vh(3,1:nsppl) = pack(symba_pla%helio%swiftest%vh(3,1:npl), discard_l_pl)
      discard_pla%rhill(1:nsppl) = pack(symba_pla%helio%swiftest%rhill(1:npl), discard_l_pl)
      discard_pla%xb(1,1:nsppl) = pack(symba_pla%helio%swiftest%xb(1,1:npl), discard_l_pl)
      discard_pla%xb(2,1:nsppl) = pack(symba_pla%helio%swiftest%xb(2,1:npl), discard_l_pl)
      discard_pla%xb(3,1:nsppl) = pack(symba_pla%helio%swiftest%xb(3,1:npl), discard_l_pl)
      discard_pla%vb(1,1:nsppl) = pack(symba_pla%helio%swiftest%vb(1,1:npl), discard_l_pl)
      discard_pla%vb(2,1:nsppl) = pack(symba_pla%helio%swiftest%vb(2,1:npl), discard_l_pl)
      discard_pla%vb(3,1:nsppl) = pack(symba_pla%helio%swiftest%vb(3,1:npl), discard_l_pl)
      if (config%lfragmentation .and. (nkpl + nfrag > npl)) then 
        symba_pla%helio%swiftest%name(1:nkpl) = pack(symba_pla%helio%swiftest%name(1:npl), .not. discard_l_pl)
        symba_pla%helio%swiftest%status(1:nkpl) = pack(symba_pla%helio%swiftest%status(1:npl), .not. discard_l_pl)
        symba_pla%helio%swiftest%mass(1:nkpl) = pack(symba_pla%helio%swiftest%mass(1:npl), .not. discard_l_pl)
        symba_pla%helio%swiftest%radius(1:nkpl) = pack(symba_pla%helio%swiftest%radius(1:npl), .not. discard_l_pl)
        symba_pla%helio%swiftest%xh(1,1:nkpl) = pack(symba_pla%helio%swiftest%xh(1,1:npl), .not. discard_l_pl)
        symba_pla%helio%swiftest%xh(2,1:nkpl) = pack(symba_pla%helio%swiftest%xh(2,1:npl), .not. discard_l_pl)
        symba_pla%helio%swiftest%xh(3,1:nkpl) = pack(symba_pla%helio%swiftest%xh(3,1:npl), .not. discard_l_pl)
        symba_pla%helio%swiftest%vh(1,1:nkpl) = pack(symba_pla%helio%swiftest%vh(1,1:npl), .not. discard_l_pl)
        symba_pla%helio%swiftest%vh(2,1:nkpl) = pack(symba_pla%helio%swiftest%vh(2,1:npl), .not. discard_l_pl)
        symba_pla%helio%swiftest%vh(3,1:nkpl) = pack(symba_pla%helio%swiftest%vh(3,1:npl), .not. discard_l_pl)
        symba_pla%helio%swiftest%rhill(1:nkpl) = pack(symba_pla%helio%swiftest%rhill(1:npl), .not. discard_l_pl)
        symba_pla%helio%swiftest%xb(1,1:nkpl) = pack(symba_pla%helio%swiftest%xb(1,1:npl), .not. discard_l_pl)
        symba_pla%helio%swiftest%xb(2,1:nkpl) = pack(symba_pla%helio%swiftest%xb(2,1:npl), .not. discard_l_pl)
        symba_pla%helio%swiftest%xb(3,1:nkpl) = pack(symba_pla%helio%swiftest%xb(3,1:npl), .not. discard_l_pl)
        symba_pla%helio%swiftest%vb(1,1:nkpl) = pack(symba_pla%helio%swiftest%vb(1,1:npl), .not. discard_l_pl)
        symba_pla%helio%swiftest%vb(2,1:nkpl) = pack(symba_pla%helio%swiftest%vb(2,1:npl), .not. discard_l_pl)
        symba_pla%helio%swiftest%vb(3,1:nkpl) = pack(symba_pla%helio%swiftest%vb(3,1:npl), .not. discard_l_pl)
        symba_pla%helio%ah(1,1:nkpl) = pack(symba_pla%helio%ah(1,1:npl), .not. discard_l_pl)
        symba_pla%helio%ah(2,1:nkpl) = pack(symba_pla%helio%ah(2,1:npl), .not. discard_l_pl)
        symba_pla%helio%ah(3,1:nkpl) =pack(symba_pla%helio%ah(3,1:npl), .not. discard_l_pl)

        call util_resize_pl(symba_pla, nkpl+nfrag, npl)

        npl = nkpl  + nfrag
        !add fragments 
        symba_pla%helio%swiftest%name(nkpl+1:npl) = pack(mergeadd_list%name(1:nmergeadd), frag_l_add)
        symba_pla%helio%swiftest%status(nkpl+1:npl) = [(active,i=1,nfrag)]!array of active status 
        symba_pla%helio%swiftest%mass(nkpl+1:npl) = pack(mergeadd_list%mass(1:nmergeadd), frag_l_add)
        symba_pla%helio%swiftest%radius(nkpl+1:npl) = pack(mergeadd_list%radius(1:nmergeadd), frag_l_add)
        symba_pla%helio%swiftest%xh(1,nkpl+1:npl) = pack(mergeadd_list%xh(1,1:nmergeadd), frag_l_add)
        symba_pla%helio%swiftest%xh(2,nkpl+1:npl) = pack(mergeadd_list%xh(2,1:nmergeadd), frag_l_add)
        symba_pla%helio%swiftest%xh(3,nkpl+1:npl) = pack(mergeadd_list%xh(3,1:nmergeadd), frag_l_add)
        symba_pla%helio%swiftest%vh(1,nkpl+1:npl) = pack(mergeadd_list%vh(1,1:nmergeadd), frag_l_add)
        symba_pla%helio%swiftest%vh(2,nkpl+1:npl) = pack(mergeadd_list%vh(2,1:nmergeadd), frag_l_add)
        symba_pla%helio%swiftest%vh(3,nkpl+1:npl) = pack(mergeadd_list%vh(3,1:nmergeadd), frag_l_add)

        do i = nkpl+1, npl
          mu = symba_pla%helio%swiftest%mass(1) + symba_pla%helio%swiftest%mass(i)
          r = sqrt(dot_product(symba_pla%helio%swiftest%xh(:,i), symba_pla%helio%swiftest%xh(:,i)))
          v2 = dot_product(symba_pla%helio%swiftest%vh(:,i), symba_pla%helio%swiftest%vh(:,i))
          energy = 0.5_DP*v2 - mu/r
          ap = -0.5_DP*mu/energy
         symba_pla%helio%swiftest%rhill(i) = ap*(((symba_pla%helio%swiftest%mass(i)/mu)/3.0_DP)**(1.0_DP/3.0_DP))
        end do

      else
        symba_pla%helio%swiftest%name(1:nkpl) = pack(symba_pla%helio%swiftest%name(1:npl), .not. discard_l_pl)
        symba_pla%helio%swiftest%status(1:nkpl) = pack(symba_pla%helio%swiftest%status(1:npl), .not. discard_l_pl)
        symba_pla%helio%swiftest%mass(1:nkpl) = pack(symba_pla%helio%swiftest%mass(1:npl), .not. discard_l_pl)
        symba_pla%helio%swiftest%radius(1:nkpl) = pack(symba_pla%helio%swiftest%radius(1:npl), .not. discard_l_pl)
        symba_pla%helio%swiftest%xh(1,1:nkpl) = pack(symba_pla%helio%swiftest%xh(1,1:npl), .not. discard_l_pl)
        symba_pla%helio%swiftest%xh(2,1:nkpl) = pack(symba_pla%helio%swiftest%xh(2,1:npl), .not. discard_l_pl)
        symba_pla%helio%swiftest%xh(3,1:nkpl) = pack(symba_pla%helio%swiftest%xh(3,1:npl), .not. discard_l_pl)
        symba_pla%helio%swiftest%vh(1,1:nkpl) = pack(symba_pla%helio%swiftest%vh(1,1:npl), .not. discard_l_pl)
        symba_pla%helio%swiftest%vh(2,1:nkpl) = pack(symba_pla%helio%swiftest%vh(2,1:npl), .not. discard_l_pl)
        symba_pla%helio%swiftest%vh(3,1:nkpl) = pack(symba_pla%helio%swiftest%vh(3,1:npl), .not. discard_l_pl)
        symba_pla%helio%swiftest%rhill(1:nkpl) = pack(symba_pla%helio%swiftest%rhill(1:npl), .not. discard_l_pl)
        symba_pla%helio%swiftest%xb(1,1:nkpl) = pack(symba_pla%helio%swiftest%xb(1,1:npl), .not. discard_l_pl)
        symba_pla%helio%swiftest%xb(2,1:nkpl) = pack(symba_pla%helio%swiftest%xb(2,1:npl), .not. discard_l_pl)
        symba_pla%helio%swiftest%xb(3,1:nkpl) = pack(symba_pla%helio%swiftest%xb(3,1:npl), .not. discard_l_pl)
        symba_pla%helio%swiftest%vb(1,1:nkpl) = pack(symba_pla%helio%swiftest%vb(1,1:npl), .not. discard_l_pl)
        symba_pla%helio%swiftest%vb(2,1:nkpl) = pack(symba_pla%helio%swiftest%vb(2,1:npl), .not. discard_l_pl)
        symba_pla%helio%swiftest%vb(3,1:nkpl) = pack(symba_pla%helio%swiftest%vb(3,1:npl), .not. discard_l_pl)
        symba_pla%helio%ah(1,1:nkpl) = pack(symba_pla%helio%ah(1,1:npl), .not. discard_l_pl)
        symba_pla%helio%ah(2,1:nkpl) = pack(symba_pla%helio%ah(2,1:npl), .not. discard_l_pl)
        symba_pla%helio%ah(3,1:nkpl) = pack(symba_pla%helio%ah(3,1:npl), .not. discard_l_pl)
        npl = nkpl
        symba_pla%helio%swiftest%nbody = npl
      end if
    end if 

    if (ldiscard_tp) then 
      nktp = 0
      nsptp = 0  

      discard_l_tp(1:ntp) = (symba_tpa%helio%swiftest%status(1:ntp) /= active)
      nsptp = count(discard_l_tp)
      nktp = ntp - nsptp

      call discard_tpa%alloc(nsptp) 

      discard_tpa%name(1:nsptp) = pack(symba_tpa%helio%swiftest%name(1:ntp), discard_l_tp)
      discard_tpa%status(1:nsptp) = pack(symba_tpa%helio%swiftest%status(1:ntp), discard_l_tp)
      discard_tpa%xh(1,1:nsptp) = pack(symba_tpa%helio%swiftest%xh(1,1:ntp), discard_l_tp)
      discard_tpa%xh(2,1:nsptp) = pack(symba_tpa%helio%swiftest%xh(2,1:ntp), discard_l_tp)
      discard_tpa%xh(3,1:nsptp) = pack(symba_tpa%helio%swiftest%xh(3,1:ntp), discard_l_tp)
      discard_tpa%vh(1,1:nsptp) = pack(symba_tpa%helio%swiftest%vh(1,1:ntp), discard_l_tp)
      discard_tpa%vh(2,1:nsptp) = pack(symba_tpa%helio%swiftest%vh(2,1:ntp), discard_l_tp)
      discard_tpa%vh(3,1:nsptp) = pack(symba_tpa%helio%swiftest%vh(3,1:ntp), discard_l_tp)
      discard_tpa%isperi(1:nsptp) = pack(symba_tpa%helio%swiftest%isperi(1:ntp), discard_l_tp)
      discard_tpa%peri(1:nsptp) = pack(symba_tpa%helio%swiftest%peri(1:ntp), discard_l_tp)
      discard_tpa%atp(1:nsptp) = pack(symba_tpa%helio%swiftest%atp(1:ntp), discard_l_tp)
      discard_tpa%xb(1,1:nsptp) = pack(symba_tpa%helio%swiftest%xb(1,1:ntp), discard_l_tp)
      discard_tpa%xb(2,1:nsptp) = pack(symba_tpa%helio%swiftest%xb(2,1:ntp), discard_l_tp)
      discard_tpa%xb(3,1:nsptp) = pack(symba_tpa%helio%swiftest%xb(3,1:ntp), discard_l_tp)
      discard_tpa%vb(1,1:nsptp) = pack(symba_tpa%helio%swiftest%vb(1,1:ntp), discard_l_tp)
      discard_tpa%vb(2,1:nsptp) = pack(symba_tpa%helio%swiftest%vb(2,1:ntp), discard_l_tp)
      discard_tpa%vb(3,1:nsptp) = pack(symba_tpa%helio%swiftest%vb(3,1:ntp), discard_l_tp)

      symba_tpa%helio%swiftest%name(1:nktp) = pack(symba_tpa%helio%swiftest%name(1:ntp), .not. discard_l_tp)
      symba_tpa%helio%swiftest%status(1:nktp) = pack(symba_tpa%helio%swiftest%status(1:ntp), .not. discard_l_tp)
      symba_tpa%helio%swiftest%xh(1,1:nktp) = pack(symba_tpa%helio%swiftest%xh(1,1:ntp), .not. discard_l_tp)
      symba_tpa%helio%swiftest%xh(2,1:nktp) = pack(symba_tpa%helio%swiftest%xh(2,1:ntp), .not. discard_l_tp)
      symba_tpa%helio%swiftest%xh(3,1:nktp) = pack(symba_tpa%helio%swiftest%xh(3,1:ntp), .not. discard_l_tp)
      symba_tpa%helio%swiftest%vh(1,1:nktp) = pack(symba_tpa%helio%swiftest%vh(1,1:ntp), .not. discard_l_tp)
      symba_tpa%helio%swiftest%vh(2,1:nktp) = pack(symba_tpa%helio%swiftest%vh(2,1:ntp), .not. discard_l_tp)
      symba_tpa%helio%swiftest%vh(3,1:nktp) = pack(symba_tpa%helio%swiftest%vh(3,1:ntp), .not. discard_l_tp)
      symba_tpa%helio%swiftest%xb(1,1:nktp) = pack(symba_tpa%helio%swiftest%xb(1,1:ntp), .not. discard_l_tp)
      symba_tpa%helio%swiftest%xb(2,1:nktp) = pack(symba_tpa%helio%swiftest%xb(2,1:ntp), .not. discard_l_tp)
      symba_tpa%helio%swiftest%xb(3,1:nktp) = pack(symba_tpa%helio%swiftest%xb(3,1:ntp), .not. discard_l_tp)
      symba_tpa%helio%swiftest%vb(1,1:nktp) = pack(symba_tpa%helio%swiftest%vb(1,1:ntp), .not. discard_l_tp)
      symba_tpa%helio%swiftest%vb(2,1:nktp) = pack(symba_tpa%helio%swiftest%vb(2,1:ntp), .not. discard_l_tp)
      symba_tpa%helio%swiftest%vb(3,1:nktp) = pack(symba_tpa%helio%swiftest%vb(3,1:ntp), .not. discard_l_tp)
      symba_tpa%helio%swiftest%isperi(1:nktp) = pack(symba_tpa%helio%swiftest%isperi(1:ntp), .not. discard_l_tp)
      symba_tpa%helio%swiftest%peri(1:nktp) = pack(symba_tpa%helio%swiftest%peri(1:ntp), .not. discard_l_tp)
      symba_tpa%helio%swiftest%atp(1:nktp) = pack(symba_tpa%helio%swiftest%atp(1:ntp), .not. discard_l_tp)
      symba_tpa%helio%ah(1,1:nktp) = pack(symba_tpa%helio%ah(1,1:ntp), .not. discard_l_tp)
      symba_tpa%helio%ah(2,1:nktp) = pack(symba_tpa%helio%ah(2,1:ntp), .not. discard_l_tp)
      symba_tpa%helio%ah(3,1:nktp) = pack(symba_tpa%helio%ah(3,1:ntp), .not. discard_l_tp)
      ntp = nktp
      symba_tpa%helio%swiftest%nbody = ntp
    end if 

  end procedure symba_rearray
end submodule s_symba_rearray
