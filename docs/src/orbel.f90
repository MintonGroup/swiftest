module orbel
    use swiftest
    private
    public :: xv2el
contains
    pure elemental subroutine xv2el(mu, px, py, pz, vx, vy, vz, a, e, inc, capom, omega, capm)
        use swiftest_classes, only : orbel_xv2el
        implicit none
        ! Arguments
        real*8, intent(in) :: mu, px, py, pz, vx, vy, vz
        real*8, intent(out) :: a, e, inc, capom, omega, capm
        !$f2py intent(in) mu
        !$f2py intent(in) px
        !$f2py intent(in) py
        !$f2py intent(in) pz
        !$f2py intent(in) vx
        !$f2py intent(in) vy
        !$f2py intent(in) vz
        !$f2py intent(out) a
        !$f2py intent(out) e
        !$f2py intent(out) inc
        !$f2py intent(out) capom
        !$f2py intent(out) omega
        !$f2py intent(out) capm
        ! Internals
        real*8, dimension(3) :: x, v
        x = [px, py, pz]
        v = [vx, vy, vz]
        call orbel_xv2el(mu, x(:), v(:), a, e, inc, capom, omega, capm)
        return
    end subroutine xv2el
end module orbel