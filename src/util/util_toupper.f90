submodule (util) s_util_toupper
contains
   module procedure util_toupper
   !! author: David A. Minton
   !!
   !! Convert string to uppercase
   !!
   !! Adapted from David E. Kaufmann's Swifter routine: util_toupper.f90
   use swiftest
   integer(I4B) :: i, length, idx

   length = len(string)
   do i = 1, length
      idx = iachar(string(i:i))
      if ((idx >= lowercase_begin) .and. (idx <= lowercase_end)) then
         idx = idx + uppercase_offset
         string(i:i) = achar(idx)
      end if
   end do

   return

   end procedure util_toupper
end submodule s_util_toupper
