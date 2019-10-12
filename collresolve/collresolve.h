/**
 * The library's public interface in C.
 *
 * Copyright (c) 2016-2017 University of Bern, Switzerland
 * Copyright (c) 2018-2019 Arizona Board of Regents
 *
 * Permission to use, copy, modify, and distribute this software for any
 * purpose with or without fee is hereby granted, provided that the above
 * copyright notice and this permission notice appear in all copies.
 *
 * THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
 * WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
 * MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
 * ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
 * WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
 * ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
 * OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
 *
 * In addition to the above license terms, if you use this library in scientific work
 * that leads to publication, you are kindly requested to cite the following publication:
 * Emsenhuber, A., Cambioni S., Asphaug, E., Gabriel, T. S. J., Schwartz, S. R., and Furfaro, R. (in prep.). Realistic On-the-fly Outcomes of Planetary Collisions II: Bringing Machine Learning to N-body Simulations.
 *
 * If you use the LS2012 model, you should also cite:
 * Leinhardt, Z. M. and Stewart, S. T. (2012). Collisions between Gravity-dominated Bodies. I. Outcome Regimes and Scaling Laws. The Astrophysical Journal, 745(1), 79. doi:10.1088/0004-637X/745/1/79 bib:2012ApJ...745...79L
 *
 * If you use the SL2012 model, you should also cite the same publication as for LS2012, and:
 * Stewart, S. T. and Leinhardt, Z. M. (2012). Collisions between Gravity-dominated Bodies. II. The Diversity of Impact Outcomes during the End Stage of Planet Formation. The Astrophysical Journal, 751(1), 32. doi:10.1088/0004-637X/751/1/32 bib:2012ApJ...751...32S
 * Genda, H., Kokubo, E., and Ida, S. (2012). Merging Criteria for Giant Impacts of Protoplanets. The Astrophysical Journal, 744(2), 137. doi:10.1088/0004-637X/744/2/137 bib:2012ApJ...744..137G
 *
 * If you use the C2019 model, you should also cite:
 * Cambioni, S., Asphaug, E., Emsenhuber, A., Gabriel, T. S. J., Furfaro, R., and Schwartz, S. R. (2019). Realistic On-the-fly Outcomes of Planetary Collisions: Machine Learning Applied to Simulations of Giant Impacts. The Astrophysical Journal, 875(1), 40. doi:10.3847/1538-4357/ab0e8a bib:2019ApJ...875...40C
 *
 * @file
 * @author Alexandre Emsenhuber
 */

#ifndef COLLRESOLVE_H
#define COLLRESOLVE_H

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Model to use for collisions
 */
enum collresolve_model {
	COLLRESOLVE_MODEL_NONE = 0,
	COLLRESOLVE_MODEL_PERFECT_MERGE, /* Always perfect merging */
	COLLRESOLVE_MODEL_LS2012, /* According to Leinhardt & Stewart (2012) */
	COLLRESOLVE_MODEL_SL2012, /* According to Stewart & Leinhardt (2012): as previously, but with the addition of the Graze-and-Merge regime from Kokubo & Ida (2012) */
	COLLRESOLVE_MODEL_C2019 /* According to Cambioni et al. (2019) */
};

/**
 * Type of collision.
 *
 * These are the possible values of the status returned by collresolve_resolve().
 * @{
 */
#define COLLRESOLVE_REGIME_MERGE              1
#define COLLRESOLVE_REGIME_DISRUPTION         2
#define COLLRESOLVE_REGIME_SUPERCATASTROPHIC  3
#define COLLRESOLVE_REGIME_GRAZE_AND_MERGE    4
#define COLLRESOLVE_REGIME_HIT_AND_RUN        5
/** @} */

/**
 * Error codes.
 *
 * @{
 */
#define COLLRESOLVE_ERROR_GENERAL             -1
#define COLLRESOLVE_ERROR_NO_CONF             -2
#define COLLRESOLVE_ERROR_INCORRECT_PARAMETER -3
#define COLLRESOLVE_ERROR_INCORRECT_MODEL     -4
#define COLLRESOLVE_ERROR_INCORRECT_UNIT      -5
#define COLLRESOLVE_ERROR_NON_CROSSING        -6
/** @} */

/**
 * General configuration for the library.
 *
 * It should be treated as an opaque type; the internal state can (and should!)
 * be modified using the different collresolve_conf_* functions.
 */
struct collresolve_conf;

/**
 * Representation of a body involved in a collision.
 *
 * The different fields should be self-explanatory.
 */
struct collresolve_body {
	double mass;
	double radius;
	double pos[3];
	double vel[3];
};

/**
 * Get an human-readable description of the given model code
 */
char* collresolve_model_desc( int );

/**
 * Get an human-readable description of the given regime code
 */
char* collresolve_regime_desc( int );

/**
 * Get an error message for the given error code
 */
char* collresolve_error_message( int );

/**
 * Create a new configuration object.
 *
 * The object is set to an invalid state; it cannot be used directly.
 */
struct collresolve_conf* collresolve_conf_new( void );

/**
 * Frees a configuration object
 */
void collresolve_conf_free( struct collresolve_conf* conf );

/**
 * Set the unit system to SI.
 */
int collresolve_conf_unit_si( struct collresolve_conf* conf );

/**
 * Set the unit system to M_sol, AU, day using the standard IAU values.
 */
int collresolve_conf_unit_msun_au_day( struct collresolve_conf* conf );

/**
 * Set the unit system to M_earth, metre, second.
 */
int collresolve_conf_unit_m_earth( struct collresolve_conf* conf );

/**
 * Set the unit system to M_sol, AU, day using Mercury's values.
 */
int collresolve_conf_unit_merc( struct collresolve_conf* conf );

/**
 * Set the collision model to use.
 *
 * model can take one of the values from the collresolve_model enum.
 */
int collresolve_conf_model( struct collresolve_conf* conf, int model );

/**
 * Set the distance separating the objects after a collision.
 *
 * This is given in terms of the sum of the body radii.
 *
 * This is needed by some algorithms to prevent the re-discovery of the same collision over and over.
 */
int collresolve_conf_sep_after( struct collresolve_conf* conf, double drel );

/**
 * Retrieve the bulk density for a given model and mass
 *
 * For now this always provide the bulk density associated with the mass of the Cambioni et al. (2019) model (i.e. COLLRESOLVE_MODEL_C2019).
 */
double collresolve_bulk_density( struct collresolve_conf* conf, double mass );

/**
 * Set a consistent body radius depending on the model and other body properties.
 *
 * For now this always provide the radius associated with the mass of the Cambioni et al. (2019) model (i.e. COLLRESOLVE_MODEL_C2019).
 */
int collresolve_body_radius( struct collresolve_conf* conf, struct collresolve_body* body );

/**
 * Collision setup function.
 *
 * This can be used when when the directions do not matter. The bodies will
 * be positioned at initial contact (distance is the sum of the radii), and
 * their relative velocity is set according to the two last parameters.
 */
void collresolve_setup( struct collresolve_conf* conf, struct collresolve_body* big, struct collresolve_body* small, double velocity, double angle );

/**
 * Similar, but the relative distance at moment of setup is provided.
 *
 * Useful when collisions aren't detected exactly at initial contact.
 */
void collresolve_setup_dist( struct collresolve_conf* conf, struct collresolve_body* big, struct collresolve_body* small, double distance, double velocity, double angle );

/**
 * Get the specific kinetic energy at initial contact.
 */
double collresolve_specific_energy( struct collresolve_conf* conf, struct collresolve_body big, struct collresolve_body small );

/**
 * Get the relative distance at initial contact.
 */
double collresolve_impact_distance( struct collresolve_conf* conf, struct collresolve_body big, struct collresolve_body small );

/**
 * Get the relative velocity at initial contact.
 */
double collresolve_impact_velocity( struct collresolve_conf* conf, struct collresolve_body big, struct collresolve_body small );

/**
 * Get the escape velocity at initial contact.
 */
double collresolve_escape_velocity( struct collresolve_conf* conf, struct collresolve_body big, struct collresolve_body small );

/**
 * Get the relative velocity at infinity. This will return a negative number if the objects are gravitationally bound.
 */
double collresolve_infinity_velocity( struct collresolve_conf* conf, struct collresolve_body big, struct collresolve_body small );

/**
 * Get the impact parameter (usually called "b").
 */
double collresolve_impact_parameter( struct collresolve_conf* conf, struct collresolve_body big, struct collresolve_body small );

/**
 * Get the impact angle (0 = head-on, pi/2 = grazing).
 *
 * This is simply the arc-sinus of the impact parameter.
 */
double collresolve_impact_angle( struct collresolve_conf* conf, struct collresolve_body big, struct collresolve_body small );

/**
 * Compute the outcome of the collision according to a model.
 *
 * "big" and "small" are the two bodies involved in the collision.
 * The former *must* be the most massive one. The desired number
 * of resulting bodies is given by "n". The last body always represents
 * the other material, so there are actually ( n + 1 ) objects in "ret".
 * "ret" must point to a pre-allocated area of minimal size of
 * ( n + 1 ) * sizeof( struct collresolve_body ).
 */
int collresolve_resolve( struct collresolve_conf* conf, struct collresolve_body big, struct collresolve_body small, int n, struct collresolve_body ret[] );

#ifdef __cplusplus
}
#endif

#endif
