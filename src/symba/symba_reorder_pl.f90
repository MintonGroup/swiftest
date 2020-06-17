submodule (symba) s_symba_reorder_pl
contains
   module procedure symba_reorder_pl
   !! author: David A. Minton
   !!
   !! Rearrange SyMBA planet arrays in order of decreasing mass
   !!
   !! Adapted from David E. Kaufmann's Swifter modules: symba_reorder_pl.f90
   use swiftest
   implicit none
   integer(I4B)                  :: i
   integer(I4B), dimension(:), allocatable   :: index
   real(DP), dimension(:), allocatable     :: mass
   real(DP), dimension(:,:), allocatable   :: symba_plwkspa
   integer(I4B), dimension(:,:), allocatable :: symba_plwkspa_id_status

! executable code
   allocate(index(npl), mass(npl))
   allocate(symba_plwkspa(12,npl))
   allocate(symba_plwkspa_id_status(2,npl))

   do i = 1, npl
      mass(i) = symba_pla%helio%swiftest%mass(i)
      symba_plwkspa_id_status(1,i) = symba_pla%helio%swiftest%name(i)
      symba_plwkspa_id_status(2,i) = symba_pla%helio%swiftest%status(i)
      symba_plwkspa(1,i) = symba_pla%helio%swiftest%mass(i)
      symba_plwkspa(2,i) = symba_pla%helio%swiftest%radius(i)
      symba_plwkspa(3,i) = symba_pla%helio%swiftest%xh(1,i)
      symba_plwkspa(4,i) = symba_pla%helio%swiftest%xh(2,i)
      symba_plwkspa(5,i) = symba_pla%helio%swiftest%xh(3,i)
      symba_plwkspa(6,i) = symba_pla%helio%swiftest%vh(1,i)
      symba_plwkspa(7,i) = symba_pla%helio%swiftest%vh(2,i)
      symba_plwkspa(8,i) = symba_pla%helio%swiftest%vh(3,i)
      symba_plwkspa(9,i) = symba_pla%helio%swiftest%rhill(i)
      symba_plwkspa(10,i) = symba_pla%helio%ah(1,i)
      symba_plwkspa(11,i) = symba_pla%helio%ah(2,i)
      symba_plwkspa(12,i) = symba_pla%helio%ah(3,i)
   end do
   call util_index(mass, index)
   write(*,*) "************ Reorder ***************"
   do i = 1, npl
      symba_pla%helio%swiftest%name(i) = symba_plwkspa_id_status(1,index(npl-i+1))
      symba_pla%helio%swiftest%status(i) = symba_plwkspa_id_status(2,index(npl-i+1))
      symba_pla%helio%swiftest%mass(i) = symba_plwkspa(1,index(npl-i+1))
      symba_pla%helio%swiftest%radius(i) = symba_plwkspa(2,index(npl-i+1))
      symba_pla%helio%swiftest%xh(1,i) = symba_plwkspa(3,index(npl-i+1))
      symba_pla%helio%swiftest%xh(2,i) = symba_plwkspa(4,index(npl-i+1))
      symba_pla%helio%swiftest%xh(3,i) = symba_plwkspa(5,index(npl-i+1))
      symba_pla%helio%swiftest%vh(1,i) = symba_plwkspa(6,index(npl-i+1))
      symba_pla%helio%swiftest%vh(2,i) = symba_plwkspa(7,index(npl-i+1))
      symba_pla%helio%swiftest%vh(3,i) = symba_plwkspa(8,index(npl-i+1))
      symba_pla%helio%swiftest%rhill(i) = symba_plwkspa(9,index(npl-i+1))
      symba_pla%helio%ah(1,i) = symba_plwkspa(10,index(npl-i+1))
      symba_pla%helio%ah(2,i) = symba_plwkspa(11,index(npl-i+1))
      symba_pla%helio%ah(3,i) = symba_plwkspa(12,index(npl-i+1))


   end do
   if (allocated(symba_plwkspa)) deallocate(symba_plwkspa)
   if (allocated(symba_plwkspa_id_status)) deallocate(symba_plwkspa_id_status)
   if (allocated(mass)) deallocate(mass)
   if (allocated(index)) deallocate(index)

   return

   end procedure symba_reorder_pl
end submodule s_symba_reorder_pl
