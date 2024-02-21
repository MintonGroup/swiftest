# cython: language_level=3, c_string_type=unicode, c_string_encoding=ascii

cdef extern from "_bindings.h":
    void bindings_c_driver(char* integrator, char* param_file_name, char* display_style) noexcept nogil

def driver(integrator, param_file_name, display_style):
    b_integrator = bytes(integrator,'ascii') + b'\x00'
    b_param_file_name = bytes(param_file_name,'ascii') + b'\x00'
    b_display_style = bytes(display_style,'ascii') + b'\x00'

    cdef:
        char* c_integrator = b_integrator 
        char* c_param_file_name = b_param_file_name 
        char* c_display_style = b_display_style 

    try:
        with nogil:
            bindings_c_driver(c_integrator, c_param_file_name, c_display_style)
    except:
        raise Warning("The Swiftest driver did not terminate normally")

    return