module io
   !! Module containing all input/output subroutine interface blocks 
   use module_parameters
   use module_interfaces

   interface
      module function io_read_param_in(inparfile) result(param)
         type(input_parameters) :: param
         character(*), intent(in)  :: inparfile
      end function io_read_param_in

      module function io_get_token(buffer, ilength, ifirst, ilast, ierr) result(token)
         integer(I4B), intent(in)     :: ilength
         integer(I4B), intent(inout)  :: ifirst
         integer(I4B), intent(out)    :: ilast, ierr
         character(len=*), intent(in) :: buffer
         character(len=len(buffer))   :: token
      end function io_get_token
         
   end interface

end module io


