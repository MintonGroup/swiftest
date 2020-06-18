submodule (util) s_util_sort_dp
contains
   module procedure util_sort_dp
   !! author: David A. Minton
   !!
   !! Sort input double precision array into ascending numerical order using Quicksort algorithm
   !!
   !! Adapted from David E. Kaufmann's Swifter modules: util_sort_dp.f90
   !! Adapted from Numerical Recipes in Fortran 90: The Art of Parallel Scientific Computing, by Press, Teukolsky,
   !!        Vetterling, and Flannery, 2nd ed., pp. 1169-70
   use swiftest
   integer(I4B), parameter :: NN = 15, NSTACK = 50
   real(DP)                :: a, dum
   integer(I4B)            :: n, k, i, j, jstack, l, r
   integer(I4B), dimension(NSTACK) :: istack

! executable code
   n = size(arr)
   jstack = 0
   l = 1
   r = n
   do
      if ((r - l) < NN) then
         do j = l + 1, r
            a = arr(j)
            do i = j - 1, l, -1
               if (arr(i) <= a) exit
               arr(i+1) = arr(i)
            end do
            arr(i+1) = a
         end do
         if (jstack == 0) return
         r = istack(jstack)
         l = istack(jstack-1)
         jstack = jstack - 2
      else
         k = (l + r)/2
         dum = arr(k); arr(k) = arr(l+1); arr(l+1) = dum
         if (arr(l) > arr(r)) then
            dum = arr(l); arr(l) = arr(r); arr(r) = dum
         end if
         if (arr(l+1) > arr(r)) then
            dum = arr(l+1); arr(l+1) = arr(r); arr(r) = dum
         end if
         if (arr(l) > arr(l+1)) then
            dum = arr(l); arr(l) = arr(l+1); arr(l+1) = dum
         end if
         i = l + 1
         j = r
         a = arr(l+1)
         do
            do
               i = i + 1
               if (arr(i) >= a) exit
            end do
            do
               j = j - 1
               if (arr(j) <= a) exit
            end do
            if (j < i) exit
            dum = arr(i); arr(i) = arr(j); arr(j) = dum
         end do
         arr(l+1) = arr(j)
         arr(j) = a
         jstack = jstack + 2
         if (jstack > NSTACK) then
            write(*, *) "Swiftest Error:"
            write(*, *) "   NSTACK too small in util_sort_I4B"
            call util_exit(failure)
         end if
         if ((r - i + 1) >= (j - l)) then
            istack(jstack) = r
            istack(jstack-1) = i
            r = j - 1
         else
            istack(jstack) = j - 1
            istack(jstack-1) = l
            l = i
         end if
      end if
   end do

   return

   end procedure util_sort_dp
end submodule s_util_sort_dp
