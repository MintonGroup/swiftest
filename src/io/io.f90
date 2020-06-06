module io
   !! Module containing all input/output subroutine interface blocks 
   use module_parameters
   implicit none

   interface
      module function io_read_param_in(inparfile) result(param)
         implicit none
         type(input_parameters) :: param
         character(*), intent(in)  :: inparfile
      end function io_read_param_in
   end interface

end module io


