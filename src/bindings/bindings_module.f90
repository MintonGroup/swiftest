module bindings_module
   use iso_c_binding, only : c_char, c_null_char, c_int
   use swiftest, only : swiftest_driver, swiftest_io_get_args, STRMAX
   implicit none
   contains

      subroutine bindings_c2f_string(c_string, f_string)
         implicit none
         ! Arguments
         character(len=1,kind=c_char),              intent(in)  :: c_string(*)
         character(len=:),             allocatable, intent(out) :: f_string
         ! Internals
         integer :: i
         character(len=STRMAX) :: tmp_string

         i=1
         tmp_string = ''
         do while(c_string(i) /= c_null_char .and. i <= STRMAX)
            tmp_string(i:i) = c_string(i)
            i=i+1
         end do

         if (i > 1) then
            f_string = trim(tmp_string)
         else
            f_string = ""
         end if

         return
      end subroutine bindings_c2f_string

      subroutine bindings_c_driver(c_integrator, c_param_file_name, c_display_style) bind(c)
         implicit none
         character(kind=c_char), dimension(*), intent(in) :: c_integrator, c_param_file_name, c_display_style
         character(len=:), allocatable :: integrator, param_file_name, display_style

         call bindings_c2f_string(c_integrator, integrator)
         call bindings_c2f_string(c_param_file_name, param_file_name)
         call bindings_c2f_string(c_display_style, display_style)

         call swiftest_io_get_args(integrator,param_file_name,display_style,from_cli=.false.)
         call swiftest_driver(integrator,param_file_name,display_style)
   end subroutine bindings_c_driver
end module bindings_module 
