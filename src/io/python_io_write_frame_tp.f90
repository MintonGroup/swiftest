!**********************************************************************************************************************************
!
!  Unit Name   : python_io_write_frame_tp
!  Unit Type   : subroutine
!  Project     : Swiftest
!  Package     : io
!  Language    : Fortran 90/95
!
!  Description : Write a frame of Swiftest output to a binary file
!
!  Input
!    Arguments : t              : time
!                Symba_tpA      : Python readable data structure for test particles
!                ntp            : number of test particles 
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
!  Invocation  : CALL python_io_write_frame_tp(t, symba_tpA, ntp, tp_outfile, out_stat)
!
!  Notes       : 
!
!**********************************************************************************************************************************
!**********************************************************************************************************************************
!  Author(s)   : David A. Minton, Carlisle Wishard, Jennifer Pouplin
!**********************************************************************************************************************************
subroutine python_io_write_frame_tp(t, symba_tpA, ntp, out_stat)

   ! modules
   use module_parameters
   use module_swiftest
   use module_helio
   use module_symba
   use module_interfaces, except_this_one => python_io_write_frame_tp
   implicit none

   ! arguments
   real(DP), intent(in)      :: t
   type(symba_tp),intent(in) :: symba_tpA
   integer, intent(in)       :: ntp
   character(*), intent(in)  :: out_stat

   ! internals
   logical(LGT), save        :: lfirst = .true.
   integer(I4B), parameter   :: lun = 88
   integer(I4B)              :: ierr 
   integer(I4B), save        :: iu = lun


   ! executable code
   if (lfirst) then
       if (out_stat == "APPEND") then
           call io_open(iu, tp_outfile, out_stat, "unformatted", ierr)
       else if (out_stat == "NEW") then
           call io_open(iu, tp_outfile, out_stat, "unformatted", ierr)
           if (ierr /= 0) then
              write(*, *) "SWIFTEST error:"
              write(*, *) " TP  Binary output file already exists"
              call util_exit(FAILURE)
           end if
       else
           call io_open(iu, tp_outfile, "REPLACE", "UNFORMATTED", ierr)
       end if
       if (ierr /= 0) then
           write(*, *) "SWIFTEST error:"
           write(*, *) "  Unable to open TP binary output file"
           call util_exit(FAILURE)
       end if
       lfirst = .false.
   else
       call io_open(iu, tp_outfile, "APPEND", "UNFORMATTED", ierr)
       if (ierr /= 0) then
            write(*, *) "SWIFTEST error:"
            write(*, *) "   Unable to open TP binary output file for append"
            call util_exit(failure)
       end if
   end if
   write(iu, iostat = ierr) t
   write(iu, iostat = ierr) ntp
   write(iu, iostat = ierr) symba_tpA%helio%swiftest%name(1:ntp)
   write(iu, iostat = ierr) symba_tpA%helio%swiftest%status(1:ntp)
   write(iu, iostat = ierr) symba_tpA%helio%swiftest%peri(1:ntp)
   write(iu, iostat = ierr) symba_tpA%helio%swiftest%atp(1:ntp)
   write(iu, iostat = ierr) symba_tpA%helio%swiftest%xh(1,1:ntp)
   write(iu, iostat = ierr) symba_tpA%helio%swiftest%xh(2,1:ntp)
   write(iu, iostat = ierr) symba_tpA%helio%swiftest%xh(3,1:ntp)
   write(iu, iostat = ierr) symba_tpA%helio%swiftest%vh(1,1:ntp)
   write(iu, iostat = ierr) symba_tpA%helio%swiftest%vh(2,1:ntp)
   write(iu, iostat = ierr) symba_tpA%helio%swiftest%vh(3,1:ntp)
   write(iu, iostat = ierr) symba_tpA%helio%swiftest%xb(1,1:ntp)
   write(iu, iostat = ierr) symba_tpA%helio%swiftest%xb(2,1:ntp)
   write(iu, iostat = ierr) symba_tpA%helio%swiftest%xb(3,1:ntp)
   write(iu, iostat = ierr) symba_tpA%helio%swiftest%vb(1,1:ntp)
   write(iu, iostat = ierr) symba_tpA%helio%swiftest%vb(2,1:ntp)
   write(iu, iostat = ierr) symba_tpA%helio%swiftest%vb(3,1:ntp)
   close(unit = iu, iostat = ierr)
   if (ierr /= 0) then
       write(*, *) "SWIFTEST error:"
       write(*, *) "   Unable to close TP binary output file"
       call util_exit(FAILURE)
   end if

return

end subroutine python_io_write_frame_tp
