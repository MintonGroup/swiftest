module bindings_module
   use iso_c_binding, only : c_char, c_null_char, c_int
   use swiftest 
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

      subroutine bindings_c_driver(c_integrator, c_param_file_name, c_display_style) bind(c, name="bindings_c_driver")
         implicit none
         character(kind=c_char), dimension(*), intent(in) :: c_integrator, c_param_file_name, c_display_style
         character(len=:), allocatable :: integrator, param_file_name, display_style

         call bindings_c2f_string(c_integrator, integrator)
         call bindings_c2f_string(c_param_file_name, param_file_name)
         call bindings_c2f_string(c_display_style, display_style)

         call swiftest_io_get_args(integrator,param_file_name,display_style,from_cli=.false.)
         call swiftest_driver(integrator,param_file_name,display_style)
      end subroutine bindings_c_driver


      subroutine bindings_orbel_el2xv(nbody, c_mu, c_a, c_e, c_inc, c_capom, c_omega, c_capm, &
                                      c_rx, c_ry, c_rz, c_vx, c_vy, c_vz) &
         bind(c, name="bindings_orbel_el2xv")
         !! author: David A. Minton
         !!
         !! Implements a bindings version of swiftest_orbel_el2xv to be called from Python via Cython
         implicit none
         ! Arguments
         integer(I4B), intent(in), value :: nbody !! The number of bodies
         type(c_ptr), intent(in), value :: c_mu !! The gravitational parameter G*(M+m)
         type(c_ptr), intent(in), value :: c_a !! The semi-major axis
         type(c_ptr), intent(in), value :: c_e !! The eccentricity
         type(c_ptr), intent(in), value :: c_inc !! The inclination
         type(c_ptr), intent(in), value :: c_capom !! The longitude of the ascending node
         type(c_ptr), intent(in), value :: c_omega !! The argument of periapsis
         type(c_ptr), intent(in), value :: c_capm !! The mean anomaly
         type(c_ptr), intent(in), value :: c_rx !! The x-component of the position vector
         type(c_ptr), intent(in), value :: c_ry !! The y-component of the position vector
         type(c_ptr), intent(in), value :: c_rz !! The z-component of the position vector
         type(c_ptr), intent(in), value :: c_vx !! The x-component of the velocity vector
         type(c_ptr), intent(in), value :: c_vy !! The y-component of the velocity vector
         type(c_ptr), intent(in), value :: c_vz !! The z-component of the velocity vector
         ! Internals
         real(DP), dimension(:), pointer :: mu, a, e, inc, capom, omega, capm, rx, ry, rz, vx, vy, vz
         integer(I4B) :: i

         if (c_associated(c_mu)) then
            call c_f_pointer(c_mu, mu, shape=[nbody])
         else
            error stop "c_mu is not associated"
         end if

         if (c_associated(c_a)) then
            call c_f_pointer(c_a, a, shape=[nbody])
         else
            error stop "c_a is not associated"
         end if

         if (c_associated(c_e)) then
            call c_f_pointer(c_e, e, shape=[nbody])
         else
            error stop "c_e is not associated"
         end if

         if (c_associated(c_inc)) then
            call c_f_pointer(c_inc, inc, shape=[nbody])
         else
            error stop "c_inc is not associated"
         end if

         if (c_associated(c_capom)) then
            call c_f_pointer(c_capom, capom, shape=[nbody])
         else
            error stop "c_capom is not associated"
         end if

         if (c_associated(c_omega)) then
            call c_f_pointer(c_omega, omega, shape=[nbody])
         else
            error stop "c_omega is not associated"
         end if

         if (c_associated(c_capm)) then
            call c_f_pointer(c_capm, capm, shape=[nbody])
         else
            error stop "c_capm is not associated"
         end if

         allocate(rx(nbody))
         if (c_associated(c_rx)) then
            call c_f_pointer(c_rx, rx, shape=[nbody])
         else
            error stop "c_rx is not associated"
         end if

         allocate(ry(nbody))
         if (c_associated(c_ry)) then
            call c_f_pointer(c_ry, ry, shape=[nbody])
         else
            error stop "c_ry is not associated"
         end if

         allocate(rz(nbody))
         if (c_associated(c_rz)) then
            call c_f_pointer(c_rz, rz, shape=[nbody])
         else
            error stop "c_rz is not associated"
         end if

         allocate(vx(nbody))
         if (c_associated(c_vx)) then
            call c_f_pointer(c_vx, vx, shape=[nbody])
         else
            error stop "c_vx is not associated"
         end if   

         allocate(vy(nbody))
         if (c_associated(c_vy)) then
            call c_f_pointer(c_vy, vy, shape=[nbody])
         else
            error stop "c_vy is not associated"
         end if

         allocate(vz(nbody))
         if (c_associated(c_vz)) then
            call c_f_pointer(c_vz, vz, shape=[nbody])
         else
            error stop "c_vz is not associated"
         end if

#ifdef DOCONLOC
         do concurrent (i = 1:nbody) shared(mu, a, e, inc, capom, omega, capm, rx, ry, rz, vx, vy, vz)
#else
         do concurrent (i = 1:nbody)
#endif
            call swiftest_orbel_el2xv(mu(i), a(i), e(i), inc(i), capom(i), omega(i), capm(i), &
                                      rx(i), ry(i), rz(i), vx(i), vy(i), vz(i))
         end do

         return
      end subroutine bindings_orbel_el2xv

      subroutine bindings_orbel_xv2el(nbody, c_mu, c_rx, c_ry, c_rz, c_vx, c_vy, c_vz, &
                                      c_a, c_e, c_inc, c_capom, c_omega, c_capm, c_varpi, c_lam, c_f, c_cape, c_capf) &
         bind(c, name="bindings_orbel_xv2el")
         !! author: David A. Minton
         !!
         !! Implements a bindings version of swiftest_orbel_xv2el to be called from Python via Cython
         implicit none
         ! Arguments
         integer(I4B), intent(in), value :: nbody !! The number of bodies
         type(c_ptr), intent(in), value :: c_mu !! The gravitational parameter G*(M+m)
         type(c_ptr), intent(in), value :: c_rx !! The x-component of the position vector
         type(c_ptr), intent(in), value :: c_ry !! The y-component of the position vector
         type(c_ptr), intent(in), value :: c_rz !! The z-component of the position vector
         type(c_ptr), intent(in), value :: c_vx !! The x-component of the velocity vector
         type(c_ptr), intent(in), value :: c_vy !! The y-component of the velocity vector
         type(c_ptr), intent(in), value :: c_vz !! The z-component of the velocity vector
         type(c_ptr), intent(in), value :: c_a !! The semi-major axis
         type(c_ptr), intent(in), value :: c_e !! The eccentricity
         type(c_ptr), intent(in), value :: c_inc !! The inclination
         type(c_ptr), intent(in), value :: c_capom !! The longitude of the ascending node
         type(c_ptr), intent(in), value :: c_omega !! The argument of periapsis
         type(c_ptr), intent(in), value :: c_capm !! The mean anomaly
         type(c_ptr), intent(in), value :: c_varpi !! The longitude of periapsis
         type(c_ptr), intent(in), value :: c_lam !! The mean longitude
         type(c_ptr), intent(in), value :: c_f !! The true anomaly
         type(c_ptr), intent(in), value :: c_cape !! The eccentric anomaly (elliptical orbits)
         type(c_ptr), intent(in), value :: c_capf !! The hyperbolic anomaly (hyperbolic orbits)
         ! Internals
         real(DP), dimension(:), pointer :: mu, rx, ry, rz, vx, vy, vz, a, e, inc, capom, omega, capm, varpi, lam, f, cape, capf
         integer(I4B) :: i

         if (c_associated(c_mu)) then
            call c_f_pointer(c_mu, mu, shape=[nbody])
         else
            error stop "c_mu is not associated"
         end if

         if (c_associated(c_rx)) then
            call c_f_pointer(c_rx, rx, shape=[nbody])
         else
            error stop "c_rx is not associated"
         end if   

         if (c_associated(c_ry)) then
            call c_f_pointer(c_ry, ry, shape=[nbody])
         else
            error stop "c_ry is not associated"
         end if

         if (c_associated(c_rz)) then
            call c_f_pointer(c_rz, rz, shape=[nbody])
         else
            error stop "c_rz is not associated"
         end if   

         if (c_associated(c_vx)) then
            call c_f_pointer(c_vx, vx, shape=[nbody])
         else
            error stop "c_vx is not associated"
         end if

         if (c_associated(c_vy)) then
            call c_f_pointer(c_vy, vy, shape=[nbody])
         else
            error stop "c_vy is not associated"
         end if

         if (c_associated(c_vz)) then
            call c_f_pointer(c_vz, vz, shape=[nbody])
         else
            error stop "c_vz is not associated"
         end if

         allocate(a(nbody))
         if (c_associated(c_a)) then
            call c_f_pointer(c_a, a, shape=[nbody])
         else
            error stop "c_a is not associated"
         end if

         allocate(e(nbody))
         if (c_associated(c_e)) then
            call c_f_pointer(c_e, e, shape=[nbody])
         else
            error stop "c_e is not associated"
         end if

         allocate(inc(nbody))
         if (c_associated(c_inc)) then
            call c_f_pointer(c_inc, inc, shape=[nbody])
         else
            error stop "c_inc is not associated"
         end if   

         allocate(capom(nbody))
         if (c_associated(c_capom)) then
            call c_f_pointer(c_capom, capom, shape=[nbody])
         else
            error stop "c_capom is not associated"
         end if

         allocate(omega(nbody))
         if (c_associated(c_omega)) then
            call c_f_pointer(c_omega, omega, shape=[nbody])
         else
            error stop "c_omega is not associated"
         end if

         allocate(capm(nbody))
         if (c_associated(c_capm)) then
            call c_f_pointer(c_capm, capm, shape=[nbody])
         else
            error stop "c_capm is not associated"
         end if

         allocate(varpi(nbody))
         if (c_associated(c_varpi)) then
            call c_f_pointer(c_varpi, varpi, shape=[nbody])
         else
            error stop "c_varpi is not associated"
         end if

         allocate(lam(nbody))
         if (c_associated(c_lam)) then
            call c_f_pointer(c_lam, lam, shape=[nbody])
         else
            error stop "c_lam is not associated"
         end if

         allocate(f(nbody))
         if (c_associated(c_f)) then
            call c_f_pointer(c_f, f, shape=[nbody])
         else
            error stop "c_f is not associated"
         end if

         allocate(cape(nbody))
         if (c_associated(c_cape)) then
            call c_f_pointer(c_cape, cape, shape=[nbody])
         else
            error stop "c_cape is not associated"
         end if

         allocate(capf(nbody))
         if (c_associated(c_capf)) then
            call c_f_pointer(c_capf, capf, shape=[nbody])
         else
            error stop "c_capf is not associated"
         end if

#ifdef DOCONLOC
         do concurrent (i = 1:nbody) shared(mu, rx, ry, rz, vx, vy, vz, a, e, inc, capom, omega, capm, varpi, lam, f, cape, capf) 
#else
         do concurrent (i = 1:nbody)
#endif
            call swiftest_orbel_xv2el(mu(i), rx(i), ry(i), rz(i), vx(i), vy(i), vz(i), &
                                      a(i), e(i), inc(i), capom(i), omega(i), capm(i), varpi(i), lam(i), f(i), cape(i), capf(i))
         end do


         return
      end subroutine bindings_orbel_xv2el

end module bindings_module 
