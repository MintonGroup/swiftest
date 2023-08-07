module pydriver_interface
   use iso_c_binding, only : c_char
   use pydriver_module, only : driver
   implicit none
   contains
      subroutine c_driver(integrator, param_file_name, display_style) bind(c)
      character(kind=c_char) :: integrator(:), param_file_name(:), display_style(:)
   end subroutine c_driver
end module pydriver_interface
