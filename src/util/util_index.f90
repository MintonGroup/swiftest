submodule (util) s_util_index
contains
   module procedure util_index
   !! author: David A. Minton
   !!
   !! Index input real array into ascending numerical order using Quicksort algorithm
   !!
   !! Adapted from David E. Kaufmann's Swifter modules: util_index.f90
   !! Adapted from Numerical Recipes in Fortran 90: The Art of Parallel Scientific Computing, by Press, Teukolsky,
   !!      Vetterling, and Flannery, 2nd ed., pp. 1173-4
   use swiftest
   integer(I4B), parameter       :: nn = 15, nstack = 50
   integer(I4B)            :: n, k, i, j, indext, jstack, l, r, dum
   integer(I4B), dimension(nstack) :: istack
   real(DP)                :: a

   n = size(arr)
   if (n /= size(index)) then
      write(*, *) "Swiftest Error:"
      write(*, *) "   array size mismatch in util_index"
      call util_exit(FAILURE)
   end if
   index = arth(1, 1, n)
   jstack = 0
   ! l is the counter ie 'the one we are at'
   l = 1
   ! r is the length of the array ie 'the total number of particles'
   r = n
   do
      if ((r - l) < nn) then
         do j = l + 1, r
            indext = index(j)
            a = arr(indext)
            do i = j - 1, l, -1
               if (arr(index(i)) <= a) exit
               index(i+1) = index(i)
            end do
            index(i+1) = indext
         end do
         if (jstack == 0) return
         r = istack(jstack)
         l = istack(jstack-1)
         jstack = jstack - 2
      else
         k = (l + r)/2
         dum = index(k); index(k) = index(l+1); index(l+1) = dum
         ! if the mass of the particle we are at in our counting is greater than the mass of the last particle then put the particle we are at above the last one
         if (arr(index(l)) > arr(index(r))) then
            dum = index(l); index(l) = index(r); index(r) = dum
         end if
         ! if the mass of the particle above the one we are at in our counting is greater than the last particle then put that particle above the last one
         if (arr(index(l+1)) > arr(index(r))) then
            dum = index(l+1); index(l+1) = index(r); index(r) = dum
         end if
         ! if the mass of teh particle we are at in our counting is greater than the one above it, then put it above the one above it
         if (arr(index(l)) > arr(index(l+1))) then
            dum = index(l); index(l) = index(l+1); index(l+1) = dum
         end if
         i = l + 1
         j = r
         indext = index(l+1)
         a = arr(indext)
         do
            do
               i = i + 1
               if (arr(index(i)) >= a) exit
            end do
            do
               j = j - 1
               if (arr(index(j)) <= a) exit
            end do
            if (j < i) exit
            dum = index(i); index(i) = index(j); index(j) = dum
         end do
         index(l+1) = index(j)
         index(j) = indext
         jstack = jstack + 2
         if (jstack > nstack) then
            write(*, *) "Swiftest Error:"
            write(*, *) "   nstack too small in util_sort"
            call util_exit(FAILURE)
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

   end procedure util_index
end submodule s_util_index
