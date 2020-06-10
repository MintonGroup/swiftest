!**********************************************************************************************************************************
!
!  Unit Name   : python_io_write_frame_pl
!  Unit Type   : subroutine
!  Project     : Swiftest
!  Package     : io
!  Language    : Fortran 90/95
!
!  Description : Write a frame of Swiftest output to a binary file
!
!  Input
!    Arguments : t              : time
!                Symba_plA      : Python readable data structure
!                python_outfile : output file name
!                out_stat     : open status for output binary file (either "APPEND" or "NEW")
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : none
!    Terminal  : error message
!    File      : none
!
!  Invocation  : CALL python_io_write_frame_pl(t, ring)
!
!  Notes       : 
!
!**********************************************************************************************************************************
!**********************************************************************************************************************************
!  Author(s)   : David A. Minton, Carlisle Wishard, Jennifer Pouplin
!**********************************************************************************************************************************
subroutine python_io_write_frame_pl(t, symba_plA, npl, out_stat)

   ! modules
   use swiftest
   use module_helio
   use module_symba
   use module_interfaces, except_this_one => python_io_write_frame_pl
   implicit none

   ! arguments
   real(DP), intent(in)      :: t
   type(symba_pl),intent(in) :: symba_plA
   integer, intent(in)       :: npl
   character(*), intent(in)  :: out_stat

   ! internals
   logical(LGT), save        :: lfirst = .true.
   integer(I4B), parameter   :: lun = 88
   integer(I4B)              :: ierr 
   integer(I4B), save        :: iu = lun
   !real(DP),dimension(count(seeds%active)) :: aseeds, Gmseeds

   ! executable code
   if (lfirst) then
       if (out_stat == "APPEND") then
           call io_open(iu, pl_outfile, out_stat, "unformatted", ierr)
       else if (out_stat == "NEW") then
           call io_open(iu, pl_outfile, out_stat, "unformatted", ierr)
           if (ierr /= 0) then
              write(*, *) "SWIFTEST error:"
              write(*, *) "   Binary output file already exists"
              call util_exit(FAILURE)
           end if
       else
           call io_open(iu, pl_outfile, "REPLACE", "UNFORMATTED", ierr)
       end if
       if (ierr /= 0) then
           write(*, *) "SWIFTEST error:"
           write(*, *) "   Unable to open binary output file"
           call util_exit(FAILURE)
       end if
       lfirst = .false.
   else
       call io_open(iu, pl_outfile, "APPEND", "UNFORMATTED", ierr)
       if (ierr /= 0) then
            write(*, *) "SWIFTEST error:"
            write(*, *) "   Unable to open binary output file for append"
            call util_exit(failure)
       end if
   end if
   write(iu, iostat = ierr) t
   write(iu, iostat = ierr) npl
   write(iu, iostat = ierr) symba_plA%helio%swiftest%name(1:npl)
   write(iu, iostat = ierr) symba_plA%helio%swiftest%status(1:npl)
   write(iu, iostat = ierr) symba_plA%helio%swiftest%mass(1:npl)
   write(iu, iostat = ierr) symba_plA%helio%swiftest%radius(1:npl)
   write(iu, iostat = ierr) symba_plA%helio%swiftest%rhill(1:npl)
   write(iu, iostat = ierr) symba_plA%helio%swiftest%xh(1,1:npl)
   write(iu, iostat = ierr) symba_plA%helio%swiftest%xh(2,1:npl)
   write(iu, iostat = ierr) symba_plA%helio%swiftest%xh(3,1:npl)
   write(iu, iostat = ierr) symba_plA%helio%swiftest%vh(1,1:npl)
   write(iu, iostat = ierr) symba_plA%helio%swiftest%vh(2,1:npl)
   write(iu, iostat = ierr) symba_plA%helio%swiftest%vh(3,1:npl)
   write(iu, iostat = ierr) symba_plA%helio%swiftest%xb(1,1:npl)
   write(iu, iostat = ierr) symba_plA%helio%swiftest%xb(2,1:npl)
   write(iu, iostat = ierr) symba_plA%helio%swiftest%xb(3,1:npl)
   write(iu, iostat = ierr) symba_plA%helio%swiftest%vb(1,1:npl)
   write(iu, iostat = ierr) symba_plA%helio%swiftest%vb(2,1:npl)
   write(iu, iostat = ierr) symba_plA%helio%swiftest%vb(3,1:npl)
   close(unit = iu, iostat = ierr)
   if (ierr /= 0) then
       write(*, *) "SWIFTEST error:"
       write(*, *) "   Unable to close binary output file"
       call util_exit(FAILURE)
   end if

return

end subroutine python_io_write_frame_pl
