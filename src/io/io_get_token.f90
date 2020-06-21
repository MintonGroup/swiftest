submodule(swiftest_classes) s_io_get_token
contains
   module procedure io_get_token
   !! author: David A. Minton
   !!
   !! Retrieves a character token from an input string. Here a token is defined as any set of contiguous non-blank characters not 
   !! beginning with or containing "!". If "!" is present, any remaining part of the buffer including the "!" is ignored
   !!
   !! Adapted from David E. Kaufmann's Swifter routine io_get_token.f90
   implicit none

   integer(I4B) :: i,ilength

   ilength = len(buffer)

   if (ifirst > ilength) then
       ilast = ifirst
       ierr = -1 !! Bad input
       token = ''
       return
   end if
   do i = ifirst, ilength
       if (buffer(i:i) /= ' ') exit
   end do
   if ((i > ilength) .or. (buffer(i:i) == '!')) then
       ifirst = i
       ilast = i
       ierr = -2 !! No valid token
       token = ''
       return
   end if
   ifirst = i
   do i = ifirst, ilength
       if ((buffer(i:i) == ' ') .or. (buffer(i:i) == '!')) exit
   end do
   ilast = i - 1
   ierr = 0

   token = buffer(ifirst:ilast)

   return

   end procedure io_get_token
end submodule s_io_get_token
