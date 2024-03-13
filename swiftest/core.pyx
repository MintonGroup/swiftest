# cython: language_level=3, c_string_type=unicode, c_string_encoding=ascii
cimport cython
cimport numpy as cnp
import numpy as np

cdef extern from "core.h":
    void bindings_c_driver(char* integrator, char* param_file_name, char* display_style) noexcept nogil
# extern void bindings_orbel_el2xv(int nbody, double *mu, double *a, double *e, double *inc, double *capom, double *omega, double *capm, double *rx, double *ry, double *rz, double *vx, double *vy, double *vz);
#extern void bindings_orbel_xv2el(int nbody, double *mu, double *rx, double *ry, double *rz, double *vx, double *vy, double *vz, double *a, double *e, double *inc, double *capom, double *omega, double *capm);     
    void bindings_orbel_el2xv(int nbody, double *mu, double *a, double *e, double *inc, double *capom, double *omega, double *capm, double *rx, double *ry, double *rz, double *vx, double *vy, double *vz) noexcept nogil
    void bindings_orbel_xv2el(int nbody, double *mu, double *rx, double *ry, double *rz, double *vx, double *vy, double *vz, double *a, double *e, double *inc, double *capom, double *omega, double *capm, double *lam, double *f, double *cape, double *capf) noexcept nogil

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

def el2xv(cnp.ndarray[cnp.float64_t, ndim=1] mu,
          cnp.ndarray[cnp.float64_t, ndim=1] a,
          cnp.ndarray[cnp.float64_t, ndim=1] e,
            cnp.ndarray[cnp.float64_t, ndim=1] inc,
            cnp.ndarray[cnp.float64_t, ndim=1] capom,
            cnp.ndarray[cnp.float64_t, ndim=1] omega,
            cnp.ndarray[cnp.float64_t, ndim=1] capm):
    """
    Convert orbital elements to state vectors

    Parameters
    ----------
    mu : array of floats
        Gravitational parameter
    a : array of floats
        Semi-major axis
    e : array of floats
        Eccentricity
    inc : array of floats
        Inclination (degrees)
    capom : array of floats
        Longitude of the ascending node (degrees)
    omega : array of floats
        Argument of periapsis (degrees)
    capm : array of floats
        Mean anomaly (degrees)

    Returns
    -------
    rx : array of floats
        x-coordinate of position vector
    ry : array of floats
        y-coordinate of position vector
    rz : array of floats
        z-coordinate of position vector
    vx : array of floats
        x-coordinate of velocity vector
    vy : array of floats
        y-coordinate of velocity vector
    vz : array of floats
        z-coordinate of velocity vector
    """ 
    if not (mu.size == a.size == e.size == inc.size == capom.size == omega.size == capm.size):
        raise ValueError("All input arrays must have the same length")

    cdef int nbody = mu.size

    # Ensure memory-contiguous numpy arrays
    mu = np.ascontiguousarray(mu, dtype=np.float64)
    a = np.ascontiguousarray(a, dtype=np.float64)
    e = np.ascontiguousarray(e, dtype=np.float64)
    inc = np.ascontiguousarray(inc, dtype=np.float64)
    capom = np.ascontiguousarray(capom, dtype=np.float64)
    omega = np.ascontiguousarray(omega, dtype=np.float64)
    capm = np.ascontiguousarray(capm, dtype=np.float64)

    # Make memory view of the numpy arrays and convert angular quantities to radians
    cdef cnp.float64_t[::1] mu_v = mu
    cdef cnp.float64_t[::1] a_v = a
    cdef cnp.float64_t[::1] e_v = e
    cdef cnp.float64_t[::1] inc_v = np.deg2rad(inc)
    cdef cnp.float64_t[::1] capom_v = np.deg2rad(capom)
    cdef cnp.float64_t[::1] omega_v = np.deg2rad(omega)
    cdef cnp.float64_t[::1] capm_v = np.deg2rad(capm) 

    # Create arrays for outputs
    cdef cnp.float64_t[::1] rx = np.empty(nbody, dtype=np.float64)
    cdef cnp.float64_t[::1] ry = np.empty(nbody, dtype=np.float64)
    cdef cnp.float64_t[::1] rz = np.empty(nbody, dtype=np.float64)
    cdef cnp.float64_t[::1] vx = np.empty(nbody, dtype=np.float64)
    cdef cnp.float64_t[::1] vy = np.empty(nbody, dtype=np.float64)
    cdef cnp.float64_t[::1] vz = np.empty(nbody, dtype=np.float64)

    try:
        with nogil:
            bindings_orbel_el2xv(nbody, &mu_v[0], &a_v[0], &e_v[0], &inc_v[0], &capom_v[0], &omega_v[0], &capm_v[0], &rx[0], &ry[0], &rz[0], &vx[0], &vy[0], &vz[0])
    except:
        raise Warning("Failure in bindings_orbel_el2xv")

    return rx, ry, rz, vx, vy, vz

def xv2el(cnp.ndarray[cnp.float64_t, ndim=1] mu,
          cnp.ndarray[cnp.float64_t, ndim=1] rx,
          cnp.ndarray[cnp.float64_t, ndim=1] ry,
          cnp.ndarray[cnp.float64_t, ndim=1] rz,
          cnp.ndarray[cnp.float64_t, ndim=1] vx,
          cnp.ndarray[cnp.float64_t, ndim=1] vy,
          cnp.ndarray[cnp.float64_t, ndim=1] vz):
    """
    Convert state vectors to orbital elements

    Parameters
    ----------
    mu : array of floats
        Gravitational parameter
    rx : array of floats
        x-coordinate of position vector
    ry : array of floats
        y-coordinate of position vector
    rz : array of floats
        z-coordinate of position vector
    vx : array of floats
        x-coordinate of velocity vector
    vy : array of floats
        y-coordinate of velocity vector
    vz : array of floats
        z-coordinate of velocity vector

    Returns
    -------
    a : array of floats
        Semi-major axis
    e : array of floats
        Eccentricity
    inc : array of floats
        Inclination (degrees)
    capom : array of floats
        Longitude of the ascending node (degrees)
    omega : array of floats
        Argument of periapsis (degrees)
    capm : array of floats
        Mean anomaly (degrees)
    lam : array of floats
        True longitude (degrees)
    f : array of floats
        True anomaly (degrees)
    cape : array of floats
        Eccentric anomaly (degrees)
    capf : array of floats
        Eccentric true anomaly (degrees)
    """ 

    if not (mu.size == rx.size == ry.size == rz.size == vx.size == vy.size == vz.size):
        raise ValueError("All input arrays must have the same length")

    cdef int nbody = mu.size

    # Ensure memory-contiguous numpy arrays
    mu = np.ascontiguousarray(mu, dtype=np.float64)
    rx = np.ascontiguousarray(rx, dtype=np.float64)
    ry = np.ascontiguousarray(ry, dtype=np.float64)
    rz = np.ascontiguousarray(rz, dtype=np.float64)
    vx = np.ascontiguousarray(vx, dtype=np.float64)
    vy = np.ascontiguousarray(vy, dtype=np.float64)
    vz = np.ascontiguousarray(vz, dtype=np.float64)

    # Make memory view of the numpy arrays
    cdef cnp.float64_t[::1] mu_v = mu
    cdef cnp.float64_t[::1] rx_v = rx
    cdef cnp.float64_t[::1] ry_v = ry
    cdef cnp.float64_t[::1] rz_v = rz
    cdef cnp.float64_t[::1] vx_v = vx
    cdef cnp.float64_t[::1] vy_v = vy
    cdef cnp.float64_t[::1] vz_v = vz

    # Create arrays for outputs
    cdef cnp.float64_t[::1] a = np.empty(nbody, dtype=np.float64)
    cdef cnp.float64_t[::1] e = np.empty(nbody, dtype=np.float64)
    cdef cnp.float64_t[::1] inc = np.empty(nbody, dtype=np.float64)
    cdef cnp.float64_t[::1] capom = np.empty(nbody, dtype=np.float64)
    cdef cnp.float64_t[::1] omega = np.empty(nbody, dtype=np.float64)
    cdef cnp.float64_t[::1] capm = np.empty(nbody, dtype=np.float64)
    cdef cnp.float64_t[::1] lam = np.empty(nbody, dtype=np.float64)
    cdef cnp.float64_t[::1] f = np.empty(nbody, dtype=np.float64)
    cdef cnp.float64_t[::1] cape = np.empty(nbody, dtype=np.float64)
    cdef cnp.float64_t[::1] capf = np.empty(nbody, dtype=np.float64)

    try:
        with nogil:
            bindings_orbel_xv2el(nbody, &mu_v[0], &rx_v[0], &ry_v[0], &rz_v[0], &vx_v[0], &vy_v[0], &vz_v[0], &a[0], &e[0], &inc[0], &capom[0], &omega[0], &capm[0], &lam[0], &f[0], &cape[0], &capf[0])
    except:
        raise Warning("Failure in bindings_orbel_xv2el")


    inc = np.rad2deg(inc)
    capom = np.rad2deg(capom)
    omega = np.rad2deg(omega)
    capm = np.rad2deg(capm)
    lam = np.rad2deg(lam)
    f = np.rad2deg(f)
    cape = np.rad2deg(cape)
    capf = np.rad2deg(capf)

    return a, e, inc, capom, omega, capm, lam, f, cape, capf


