!**********************************************************************************************************************************
!
!  Unit Name   : ringmoons_io_write_frame
!  Unit Type   : subroutine
!  Project     : Swifter
!  Package     : io
!  Language    : Fortran 90/95
!
!  Description : Write a frame of a RINGMOONS output to a binary file
!
!  Input
!    Arguments : t            : time
!                ring         : RINGMOONS data structure
!                ring_outfile : output file name
!                out_stat     : open status for output binary file (either "APPEND" or "NEW")
!    Terminal  : none
!    File      : none
!
!  Output
!    Arguments : none
!    Terminal  : error message
!    File      : none
!
!  Invocation  : CALL ringmoons_io_write_frame(t, ring)
!
!  Notes       : 
!
!**********************************************************************************************************************************
!**********************************************************************************************************************************
!  Author(s)   : David A. Minton
!**********************************************************************************************************************************
subroutine ringmoons_io_write_frame(t, ring, seeds, ring_outfile, out_stat)

   ! modules
   use module_parameters
   use module_ringmoons
   use module_ringmoons_interfaces, except_this_one => ringmoons_io_write_frame
   implicit none

   ! arguments
   real(DP), intent(in)            :: t
   type(ringmoons_ring),intent(in) :: ring
   type(ringmoons_seeds),intent(in) :: seeds
   character(*), intent(in)  :: ring_outfile, out_stat

   ! internals
   logical(LGT), save        :: lfirst = .true.
   integer(I4B), parameter   :: lun = 88
   integer(I4B)              :: i,ierr 
   integer(I4B), save        :: iu = lun
   real(DP),dimension(count(seeds%active)) :: aseeds, Gmseeds

   ! executable code
   if (lfirst) then
       if (out_stat == "APPEND") then
           call io_open(iu, ring_outfile, out_stat, "unformatted", ierr)
       else if (out_stat == "NEW") then
           call io_open(iu, ring_outfile, out_stat, "unformatted", ierr)
           if (ierr /= 0) then
              write(*, *) "RINGMOONS error:"
              write(*, *) "   Binary output file already exists"
              call util_exit(FAILURE)
           end if
       else
           call io_open(iu, ring_outfile, "REPLACE", "UNFORMATTED", ierr)
       end if
       if (ierr /= 0) then
           write(*, *) "RINGMOONS error:"
           write(*, *) "   Unable to open binary output file"
           call util_exit(FAILURE)
       end if
       lfirst = .false.
   else
       call io_open(iu, ring_outfile, "APPEND", "UNFORMATTED", ierr)
       if (ierr /= 0) then
            write(*, *) "RINGMOONS error:"
            write(*, *) "   Unable to open binary output file for append"
            call util_exit(failure)
       end if
   end if
   write(iu, iostat = ierr) t
   write(iu, iostat = ierr) ring%N
   write(iu, iostat = ierr) ring%r
   write(iu, iostat = ierr) ring%Gsigma
   write(iu, iostat = ierr) ring%nu
   !write(iu, iostat = ierr) seeds%N
   !write(iu, iostat = ierr) seeds%a
   !write(iu, iostat = ierr) seeds%Gm
   write(iu, iostat = ierr) count(seeds%active)
   aseeds = pack(seeds%a,seeds%active) 
   Gmseeds = pack(seeds%Gm,seeds%active) 
   write(iu, iostat = ierr) aseeds
   write(iu, iostat = ierr) Gmseeds
   close(unit = iu, iostat = ierr)
   if (ierr /= 0) then
       write(*, *) "RINGMOONS error:"
       write(*, *) "   Unable to close binary output file"
       call util_exit(FAILURE)
   end if

return

end subroutine ringmoons_io_write_frame
