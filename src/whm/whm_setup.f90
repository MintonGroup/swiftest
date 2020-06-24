submodule(whm) s_whm_setup
contains
   module procedure whm_setup(npl, ntp, whm_pla, whm_tpa, whm_pl1p, whm_tp1p, swifter_pl1p, swifter_tp1p)
   !! author: David A. Minton
   !!
   !! Set up pointers within WHM and SWIFTER planet and test particle structure linked-lists
   !!
   !! 
   !! Adapted from David E. Kaufmann's Swifter routine whm_setup.f90
   use swiftest
   implicit none
   integer(I4B)          :: i
   type(whm_pl), pointer   :: whm_plp
   type(whm_tp), pointer   :: whm_tpp
   type(swifter_pl), pointer :: swifter_plp
   type(swifter_tp), pointer :: swifter_tpp

! executable code
   whm_pl1p => whm_pla(1)
   swifter_pl1p => whm_pla(1)%swifter
   nullify(whm_pl1p%prevp)
   nullify(swifter_pl1p%prevp)
   if (npl == 1) then
      nullify(whm_pl1p%nextp)
      nullify(swifter_pl1p%nextp)
   else
      whm_pl1p%nextp => whm_pla(2)
      swifter_pl1p%nextp => whm_pla(2)%swifter
      do i = 2, npl - 1
         whm_pla(i)%prevp => whm_pla(i-1)
         whm_pla(i)%nextp => whm_pla(i+1)
         swifter_plp => whm_pla(i)%swifter
         swifter_plp%prevp => whm_pla(i-1)%swifter
         swifter_plp%nextp => whm_pla(i+1)%swifter
      end do
      whm_pla(npl)%prevp => whm_pla(npl-1)
      whm_plp => whm_pla(npl)
      nullify(whm_plp%nextp)
      swifter_plp => whm_pla(npl)%swifter
      swifter_plp%prevp => whm_pla(npl-1)%swifter
      nullify(swifter_plp%nextp)
   end if
   nullify(whm_tp1p)
   nullify(swifter_tp1p)
   if (ntp > 0) then
      whm_tp1p => whm_tpa(1)
      swifter_tp1p => whm_tpa(1)%swifter
      nullify(whm_tp1p%prevp)
      nullify(swifter_tp1p%prevp)
      if (ntp == 1) then
         nullify(whm_tp1p%nextp)
         nullify(swifter_tp1p%nextp)
      else
         whm_tp1p%nextp => whm_tpa(2)
         swifter_tp1p%nextp => whm_tpa(2)%swifter
         do i = 2, ntp - 1
            whm_tpa(i)%prevp => whm_tpa(i-1)
            whm_tpa(i)%nextp => whm_tpa(i+1)
            swifter_tpp => whm_tpa(i)%swifter
            swifter_tpp%prevp => whm_tpa(i-1)%swifter
            swifter_tpp%nextp => whm_tpa(i+1)%swifter
         end do
         whm_tpa(ntp)%prevp => whm_tpa(ntp-1)
         whm_tpp => whm_tpa(ntp)
         nullify(whm_tpp%nextp)
         swifter_tpp => whm_tpa(ntp)%swifter
         swifter_tpp%prevp => whm_tpa(ntp-1)%swifter
         nullify(swifter_tpp%nextp)
      end if
   end if

   return

   end procedure whm_setup
end submodule s_whm_setup
