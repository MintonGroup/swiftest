# cython: language_level=3, c_string_type=unicode, c_string_encoding=ascii
cimport cython
cimport numpy as cnp
import numpy as np
import warnings
import traceback




cdef extern from "core.h":
    void bindings_c_driver(char* integrator, char* param_file_name, char* display_style) nogil
    void bindings_orbel_el2xv(int nbody, double *mu, double *a, double *e, double *inc, double *capom, double *omega, double *capm, double *rx, double *ry, double *rz, double *vx, double *vy, double *vz) 
    void bindings_orbel_xv2el(int nbody, double *mu, double *rx, double *ry, double *rz, double *vx, double *vy, double *vz, double *a, double *e, double *inc, double *capom, double *omega, double *capm, double *lam, double *f, double *cape, double *capf) 

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
    except Exception as e:  
        traceback_details = traceback.format_exc()
        warning_message = f"An unexpected error occurred in bindings_c_driver: {e}\n{traceback_details}"
        warnings.warn(warning_message, RuntimeWarning)

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
    rh : (n,3) array of floats
        position vector
    vh : (n,3) array of floats
    """ 
    if not (mu.size == a.size == e.size == inc.size == capom.size == omega.size == capm.size):
        raise ValueError("All input arrays must have the same length")

    cdef int nbody = mu.size

    # Ensure memory-contiguous numpy arrays
    mu = np.ascontiguousarray(mu, dtype=np.float64)
    a = np.ascontiguousarray(a, dtype=np.float64)
    e = np.ascontiguousarray(e, dtype=np.float64)
    inc_rad = np.ascontiguousarray(np.deg2rad(inc), dtype=np.float64)
    capom_rad = np.ascontiguousarray(np.deg2rad(capom), dtype=np.float64)
    omega_rad = np.ascontiguousarray(np.deg2rad(omega), dtype=np.float64)
    capm_rad = np.ascontiguousarray(np.deg2rad(capm), dtype=np.float64)

    # Make memory view of the numpy arrays and convert angular quantities to radians
    cdef cnp.float64_t[::1] mu_v = mu
    cdef cnp.float64_t[::1] a_v = a
    cdef cnp.float64_t[::1] e_v = e
    cdef cnp.float64_t[::1] inc_v = inc_rad
    cdef cnp.float64_t[::1] capom_v = capom_rad
    cdef cnp.float64_t[::1] omega_v = omega_rad
    cdef cnp.float64_t[::1] capm_v = capm_rad

    # Create arrays for outputs
    _rx = np.empty(nbody, dtype=np.float64)
    _ry = np.empty(nbody, dtype=np.float64)
    _rz = np.empty(nbody, dtype=np.float64)
    _vx = np.empty(nbody, dtype=np.float64)
    _vy = np.empty(nbody, dtype=np.float64)
    _vz = np.empty(nbody, dtype=np.float64)

    cdef cnp.float64_t[::1] rx = _rx 
    cdef cnp.float64_t[::1] ry = _ry
    cdef cnp.float64_t[::1] rz = _rz
    cdef cnp.float64_t[::1] vx = _vx
    cdef cnp.float64_t[::1] vy = _vy
    cdef cnp.float64_t[::1] vz = _vz

    try:
        bindings_orbel_el2xv(nbody, &mu_v[0], &a_v[0], &e_v[0], &inc_v[0], &capom_v[0], &omega_v[0], &capm_v[0], &rx[0], &ry[0], &rz[0], &vx[0], &vy[0], &vz[0])
    except Exception as e:  
        traceback_details = traceback.format_exc()
        warning_message = f"An unexpected error occurred in bindings_orbel_el2xv: {e}\n{traceback_details}"
        warnings.warn(warning_message, RuntimeWarning)

    rh = np.stack((rx, ry, rz), axis=1)
    vh = np.stack((vx, vy, vz), axis=1)

    return rh, vh


def xv2el(cnp.ndarray[cnp.float64_t, ndim=1] mu,
          cnp.ndarray[cnp.float64_t, ndim=2] rh,
          cnp.ndarray[cnp.float64_t, ndim=2] vh):
    """
    Convert state vectors to orbital elements

    Parameters
    ----------
    mu : array of floats
        Gravitational parameter
    rh : (N,3) array of floats
        position vector
    vh : (N,3) array of floats

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

    # Ensure the input arrays are compatible
    if not (mu.size == rh.shape[0] == vh.shape[0]):
        raise ValueError("mu, rh, and vh must have compatible dimensions.")

    cdef int nbody = mu.size

    # Ensure memory-contiguous numpy arrays
    mu = np.ascontiguousarray(mu, dtype=np.float64)

    # Extracting individual components from 2D arrays
    rx, ry, rz = rh[:,0], rh[:,1], rh[:,2]
    vx, vy, vz = vh[:,0], vh[:,1], vh[:,2]

    # Make memory view of the numpy arrays
    cdef cnp.float64_t[::1] mu_v = mu
    cdef cnp.float64_t[::1] rx_v = rx
    cdef cnp.float64_t[::1] ry_v = ry
    cdef cnp.float64_t[::1] rz_v = rz
    cdef cnp.float64_t[::1] vx_v = vx
    cdef cnp.float64_t[::1] vy_v = vy
    cdef cnp.float64_t[::1] vz_v = vz

    # Create arrays for outputs
    _a = np.empty(nbody, dtype=np.float64)
    _e = np.empty(nbody, dtype=np.float64)
    _inc = np.empty(nbody, dtype=np.float64)
    _capom = np.empty(nbody, dtype=np.float64)
    _omega = np.empty(nbody, dtype=np.float64)
    _capm = np.empty(nbody, dtype=np.float64)
    _lam = np.empty(nbody, dtype=np.float64)
    _f = np.empty(nbody, dtype=np.float64)
    _cape = np.empty(nbody, dtype=np.float64)
    _capf = np.empty(nbody, dtype=np.float64)

    cdef cnp.float64_t[::1] a = _a
    cdef cnp.float64_t[::1] e = _e
    cdef cnp.float64_t[::1] inc = _inc
    cdef cnp.float64_t[::1] capom = _capom
    cdef cnp.float64_t[::1] omega = _omega
    cdef cnp.float64_t[::1] capm = _capm
    cdef cnp.float64_t[::1] lam = _lam
    cdef cnp.float64_t[::1] f = _f
    cdef cnp.float64_t[::1] cape = _cape
    cdef cnp.float64_t[::1] capf = _capf

    try:
        bindings_orbel_xv2el(nbody, &mu_v[0], &rx_v[0], &ry_v[0], &rz_v[0], &vx_v[0], &vy_v[0], &vz_v[0], &a[0], &e[0], &inc[0], &capom[0], &omega[0], &capm[0], &lam[0], &f[0], &cape[0], &capf[0])
    except Exception as e:  
        traceback_details = traceback.format_exc()
        warning_message = f"An unexpected error occurred in bindings_orbel_xv2el: {e}\n{traceback_details}"
        warnings.warn(warning_message, RuntimeWarning)

    inc = np.rad2deg(inc)
    capom = np.rad2deg(capom)
    omega = np.rad2deg(omega)
    capm = np.rad2deg(capm)
    lam = np.rad2deg(lam)
    f = np.rad2deg(f)
    cape = np.rad2deg(cape)
    capf = np.rad2deg(capf)

    return a, e, inc, capom, omega, capm, lam, f, cape, capf


