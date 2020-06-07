submodule (io) s_io_getn
contains
   module procedure io_getn
   !! author: David A. Minton
   !!
   !! Read the number of planets and test particles from respective input files
   !! param%nplmax (param%ntpmax) is reset to npl (ntp) if the latter exceeds the former
   !! Adapted from David E. Kaufmann's Swifter routine io_getn.f90
   !! Adapted from Martin Duncan's Swift routine io_init_param.f
   implicit none

   integer(I4B), parameter :: LUN = 7           !! File unit number
   integer(I4B)            :: ierr              !! Error code
   integer(I4B)            :: npl = 0,ntp  = 0  !! Temporary variable to store values in

   npl = 0
   open(unit = LUN, file = param%inplfile, status = 'old', iostat = ierr)
   read(LUN, *) npl
   close(LUN)
   if (npl < 1) then
      write(*, *) "Error: the number of planets, ", swiftest_plA%npl, ","
      write(*, *) "       must be at least 1"
      call util_exit(FAILURE)
   else if ((npl > param%nplmax) .and. .not. (param%nplmax < 0)) then
      write(*, *) "Warning: the number of planets, ", npl, ","
      write(*, *) "         exceeds the specified maximum number of planets, ", param%nplmax
      write(*, *) "         ...resetting param%nplmax to ", npl
      param%nplmax = npl
   else if (param%nplmax < 0) then
      param%nplmax = npl
   end if

   ntp = 0
   if (param%intpfile /= "") then
      open(unit = LUN, file = param%intpfile, status = 'old', iostat = ierr)
      read(LUN, *) ntp
      close(LUN)
   end if
   if (ntp < 0) then
      write(*, *) "Error: the number of test particles, ", ntp, ","
      write(*, *) "       must be at least 0"
      call util_exit(FAILURE)
   else if ((ntp > param%ntpmax) .and. .not. (param%ntpmax < 0)) then
      write(*, *) "Warning: the number of test particles, ", ntp, ","
      write(*, *) "         exceeds the specified maximum number of test particles, ", param%ntpmax
      write(*, *) "         ...resetting param%ntpmax to ", ntp
      param%ntpmax = ntp
   else if (param%ntpmax < 0) then
      param%ntpmax = ntp
   end if

   swiftest_plA%npl = npl
   swiftest_tpA%ntp = ntp
   return

   end procedure io_getn
end submodule s_io_getn
