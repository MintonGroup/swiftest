/**
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

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdlib.h>
#include <math.h>

#include "collresolve.h"

/**
 * Cambioni et al. (2019) surrogate model.
 */
#include "cambioni2019/rt_nonfinite.h"
#include "cambioni2019/accretion_efficiency.h"
#include "cambioni2019/collision_classifier.h"
#include "cambioni2019/orbital_hnr.h"

/**
 * The gravity constant in SI units.
 * Value is from CODATA 2014: G = 6.67408(31) * 10^-11  kg^-1  m^3  s^-2
 */
#define GRAV_SI 6.67408e-11

/**
 * An astronomical unit expressed in SI unit; definition of 2012
 */
#define AU_SI 149597870700.

/**
 * A day in SI units
 */
#define DAY_SI 86400.

/**
 * Following values are from IAU 2015 Resolution B3.
 * They are conversion factors only; not true values.
 */
#define R_E_E  6.3781e6
#define R_E_P  6.3568e6
#define R_J_E  7.1492e7
#define R_J_P  6.6854e7
#define R_S    6.957e8
#define GM_E   3.986004e14
#define GM_J   1.2668653e17
#define GM_S   1.3271244e20

/**
 * The configuration object.
 *
 * It is not part of the public interface.
 */

struct collresolve_conf {
	size_t model;
	double G;
	double mEarth;
	double refDens;
	double drel;
};

/**
 * Common quantities for a collision
 */
struct collresolve_quant {
	double total_mass;
	double reduced_mass;
	double total_radius;

	double dpos_sq;
	double dvel_sq;

	double specific_energy;
	double impact_velocity_sq;
	double escape_velocity_sq;
	double impact_parameter;
};

/**
 * Function to get an user-friendly error message
 */
char* collresolve_model_desc( int model ) {
	if ( model == COLLRESOLVE_MODEL_PERFECT_MERGE ) {
		return "Perfect mergering";
	} else if ( model == COLLRESOLVE_MODEL_LS2012 ) {
		return "Leinhardt & Stewart (2012)";
	} else if ( model == COLLRESOLVE_MODEL_SL2012 ) {
		return "Stewart & Leinhardt (2012)";
	} else if ( model == COLLRESOLVE_MODEL_C2019 ) {
		return "Cambioni et al. (2019)";
	} else {
		return NULL;
	}
}

/**
 * Function to get an user-friendly error message
 */
char* collresolve_regime_desc( int regime ) {
	if ( regime == COLLRESOLVE_REGIME_MERGE ) {
		return "Merger";
	} else if ( regime == COLLRESOLVE_REGIME_DISRUPTION ) {
		return "Disruption";
	} else if ( regime == COLLRESOLVE_REGIME_SUPERCATASTROPHIC ) {
		return "Super catastrophic";
	} else if ( regime == COLLRESOLVE_REGIME_GRAZE_AND_MERGE ) {
		return "Graze and merge";
	} else if ( regime == COLLRESOLVE_REGIME_HIT_AND_RUN ) {
		return "Hit and run";
	} else {
		return NULL;
	}
}

/**
 * Function to get an user-friendly error message
 */
char* collresolve_error_message( int error ) {
	if ( error == COLLRESOLVE_ERROR_GENERAL ) {
		return "General error";
	} else if ( error == COLLRESOLVE_ERROR_NO_CONF ) {
		return "No configuration set";
	} else if ( error == COLLRESOLVE_ERROR_INCORRECT_PARAMETER ) {
		return "Incorrect parameter provided";
	} else if ( error == COLLRESOLVE_ERROR_INCORRECT_MODEL ) {
		return "Incorrect model";
	} else if ( error == COLLRESOLVE_ERROR_INCORRECT_UNIT ) {
		return "Incorrect unit system configured";
	} else if ( error == COLLRESOLVE_ERROR_NON_CROSSING ) {
		return "Non-crossing orbit";
	} else {
		return "Unknown error code";
	}
}

/**
 * Functions to set the configuration
 */

struct collresolve_conf* collresolve_conf_new() {
	struct collresolve_conf* ret = calloc( 1, sizeof( struct collresolve_conf ) );
	ret->model = COLLRESOLVE_MODEL_NONE;
	return ret;
}

void collresolve_conf_free( struct collresolve_conf* conf ) {
	free( conf );
}

int collresolve_conf_unit_si( struct collresolve_conf* conf ) {
	if ( conf == NULL ) {
		return COLLRESOLVE_ERROR_NO_CONF;
	}

	conf->G = GRAV_SI;
	conf->mEarth = GM_E / GRAV_SI;
	conf->refDens = 1000.;

	return 1;
}

int collresolve_conf_unit_msun_au_day( struct collresolve_conf* conf ) {
	if ( conf == NULL ) {
		return COLLRESOLVE_ERROR_NO_CONF;
	}

	const double mass = GM_S / GRAV_SI; /* Mass of the sun; the IAU way */
	const double dist = AU_SI;
	const double time = DAY_SI;
	const double vol = dist * dist * dist;

	conf->G = GM_S * time * time / vol;
	conf->mEarth = GM_E / GM_S;
	conf->refDens = 1000. * vol / mass;

	return 1;
}

int collresolve_conf_unit_m_earth( struct collresolve_conf* conf ) {
	if ( conf == NULL ) {
		return COLLRESOLVE_ERROR_NO_CONF;
	}

	const double mass = GM_E / GRAV_SI; /* Mass of the earth; the IAU way */

	conf->G = GM_E;
	conf->mEarth = mass;
	conf->refDens = 1000. / mass;

	return 1;
}

int collresolve_conf_unit_merc( struct collresolve_conf* conf ) {
	if ( conf == NULL ) {
		return COLLRESOLVE_ERROR_NO_CONF;
	}

	const double mass = 1.9891e30;
	const double dist = 1.4959787e11;
	const double vol = dist * dist * dist;

	conf->G = 2.959122082855911e-4;
	conf->mEarth = GM_E / GM_S;
	conf->refDens = 1000. * vol / mass;

	return 1;
}

int collresolve_conf_model( struct collresolve_conf* conf, int model ) {
	if ( conf == NULL ) {
		return COLLRESOLVE_ERROR_NO_CONF;
	}

	if ( model == COLLRESOLVE_MODEL_NONE || model == COLLRESOLVE_MODEL_PERFECT_MERGE || model == COLLRESOLVE_MODEL_LS2012 || model == COLLRESOLVE_MODEL_SL2012 || model == COLLRESOLVE_MODEL_C2019 ) {
		conf->model = model;
		return 1;
	} else {
		return COLLRESOLVE_ERROR_INCORRECT_MODEL;
	}
}

int collresolve_conf_sep_after( struct collresolve_conf* conf, double drel ) {
	if ( conf == NULL ) {
		return COLLRESOLVE_ERROR_NO_CONF;
	}

	if ( drel >= 1. ) {
		conf->drel = drel;
		return 1;
	} else {
		return COLLRESOLVE_ERROR_INCORRECT_PARAMETER;
	}
}

/**
 * Body setup functions
 * @{
 */

/**
 * Retrieve the bulk density for a given model and mass
 *
 * For now this always provide the bulk density associated with the mass of the Cambioni et al. (2019) model (i.e. MODEL_C2019).
 */
double collresolve_bulk_density( struct collresolve_conf* conf, double mass ) {
	if ( conf == NULL ) {
		return (double)COLLRESOLVE_ERROR_NO_CONF;
	}

	if ( conf->G == 0. ) {
		return (double)COLLRESOLVE_ERROR_INCORRECT_UNIT;
	}

	double mscale = mass / conf->mEarth;

	/* Masses in Earth mass */
	double ms[] = {
		1.0e-3,
		2.0e-3,
		3.5e-3,
		5.0e-3,
		7.0e-3,
		9.0e-3,
		1.0e-2,
		2.0e-2,
		3.5e-2,
		5.0e-2,
		7.0e-2,
		9.0e-2,
		1.0e-1,
		2.0e-1,
		3.5e-1,
		5.0e-1,
		7.0e-1,
		9.0e-1,
		1.0e0
	};

	/* Corresponding density (in terms of the reference density) */
	double ds[] = {
		3.1017,
		3.1168,
		3.1348,
		3.1499,
		3.1688,
		3.1850,
		3.1928,
		3.2582,
		3.3344,
		3.3986,
		3.4688,
		3.5303,
		3.5605,
		3.8758,
		4.2942,
		4.5702,
		4.8341,
		5.0417,
		5.1251
	};

	if ( mscale <= ms[ 0 ] ) {
		return ds[ 0 ] * conf->refDens;
	}

	int i;
	for ( i = 0; i < 19; i++ ) {
		if ( mscale >= ms[ i ] && mscale < ms[ i + 1 ] ) {
			double lm = log( mscale );
			double ll = log( ms[ i ] );
			double lh = log( ms[ i + 1 ] );
			double f = ( lm - ll ) / ( lh - ll );
			return ( ds[ i + 1 ] * f + ds[ i ] * ( 1. - f ) ) * conf->refDens;
		}
	}

	/* Extrapolation using the low last data points. */
	double lm = log( mscale );
	double ll = log( ms[ 17 ] );
	double lh = log( ms[ 18 ] );
	double f = ( lm - ll ) / ( lh - ll );
	return ( ds[ 18 ] * f + ds[ 17 ] * ( 1. - f ) ) * conf->refDens;
}

/**
 * Set a consistent body radius depending on the model and other body properties.
 *
 * For now this always provide the radius associated with the mass of the Cambioni et al. (2019) model (i.e. MODEL_C2019).
 */
int collresolve_body_radius( struct collresolve_conf* conf, struct collresolve_body* body ) {
	if ( conf == NULL ) {
		return COLLRESOLVE_ERROR_NO_CONF;
	}

	if ( body == NULL ) {
		return COLLRESOLVE_ERROR_INCORRECT_PARAMETER;
	}

	double dens = collresolve_bulk_density( conf, body->mass );

	if ( dens < 0 ) {
		return (int)( dens - 0.5 );
	}

	body->radius = cbrt( 3. / 4. / M_PI * body->mass / dens );

	return 1;
}

/** @} */

/**
 * Collision setup functions
 * @{
 */

void collresolve_setup( struct collresolve_conf* conf, struct collresolve_body* big, struct collresolve_body* small, double velocity, double angle ) {
	const double total_radius = big->radius + small->radius;

	collresolve_setup_dist( conf, big, small, total_radius, velocity, angle );
}

void collresolve_setup_dist( struct collresolve_conf* conf, struct collresolve_body* big, struct collresolve_body* small, double distance, double velocity, double angle ) {
	const double total_mass = big->mass + small->mass;

	const double big_factor = -small->mass / total_mass;
	const double small_factor = big->mass / total_mass;

	const double dpos_x = distance;
	const double dpos_y = 0.;
	const double dpos_z = 0.;

	const double dvel_x = -cos( angle ) * velocity;
	const double dvel_y = sin( angle ) * velocity;
	const double dvel_z = 0.;

	big->pos[0] = dpos_x * big_factor;
	big->pos[1] = dpos_y * big_factor;
	big->pos[2] = dpos_z * big_factor;
	big->vel[0] = dvel_x * big_factor;
	big->vel[1] = dvel_y * big_factor;
	big->vel[2] = dvel_z * big_factor;

	small->pos[0] = dpos_x * small_factor;
	small->pos[1] = dpos_y * small_factor;
	small->pos[2] = dpos_z * small_factor;
	small->vel[0] = dvel_x * small_factor;
	small->vel[1] = dvel_y * small_factor;
	small->vel[2] = dvel_z * small_factor;
}

/**
 * @}
 */

/**
 * Functions for quantities
 * @{
 */

void collresolve_quant_init( struct collresolve_quant* pquant ) {
	const double nan = strtod( "NAN", NULL );

	pquant->total_mass = nan;
	pquant->reduced_mass = nan;
	pquant->total_radius = nan;

	pquant->dpos_sq = nan;
	pquant->dvel_sq = nan;

	pquant->specific_energy = nan;
	pquant->impact_velocity_sq = nan;
	pquant->escape_velocity_sq = nan;
	pquant->impact_parameter = nan;
}

#define INIT_QUANT struct collresolve_quant quant; \
	struct collresolve_quant* pquant = &quant; \
	collresolve_quant_init( pquant ); \
	struct collresolve_body* pbig = &big; \
	struct collresolve_body* psmall = &small;

double collresolve_quant_total_mass( struct collresolve_conf* conf, struct collresolve_body* pbig, struct collresolve_body* psmall, struct collresolve_quant* pquant ) {
	if ( pquant->total_mass != pquant->total_mass ) {
		pquant->total_mass = pbig->mass + psmall->mass;
	}

	return pquant->total_mass;
}

#define TOTAL_MASS collresolve_quant_total_mass( conf, pbig, psmall, pquant )

double collresolve_quant_reduced_mass( struct collresolve_conf* conf, struct collresolve_body* pbig, struct collresolve_body* psmall, struct collresolve_quant* pquant ) {
	if ( pquant->reduced_mass != pquant->reduced_mass ) {
		pquant->reduced_mass = ( pbig->mass * psmall->mass ) / TOTAL_MASS;
	}

	return pquant->reduced_mass;
}

#define REDUCED_MASS collresolve_quant_reduced_mass( conf, pbig, psmall, pquant )

double collresolve_quant_total_radius( struct collresolve_conf* conf, struct collresolve_body* pbig, struct collresolve_body* psmall, struct collresolve_quant* pquant ) {
	if ( pquant->total_radius != pquant->total_radius ) {
		pquant->total_radius = pbig->radius + psmall->radius;
	}

	return pquant->total_radius;
}

#define TOTAL_RADIUS collresolve_quant_total_radius( conf, pbig, psmall, pquant )

double collresolve_quant_dpos_sq( struct collresolve_conf* conf, struct collresolve_body* pbig, struct collresolve_body* psmall, struct collresolve_quant* pquant ) {
	if ( pquant->dpos_sq != pquant->dpos_sq ) {
		pquant->dpos_sq = ( pbig->pos[0] - psmall->pos[0] ) * ( pbig->pos[0] - psmall->pos[0] ) + ( pbig->pos[1] - psmall->pos[1] ) * ( pbig->pos[1] - psmall->pos[1] ) + ( pbig->pos[2] - psmall->pos[2] ) * ( pbig->pos[2] - psmall->pos[2] );
	}

	return pquant->dpos_sq;
}

#define DPOS_SQ collresolve_quant_dpos_sq( conf, pbig, psmall, pquant )

double collresolve_quant_dvel_sq( struct collresolve_conf* conf, struct collresolve_body* pbig, struct collresolve_body* psmall, struct collresolve_quant* pquant ) {
	if ( pquant->dvel_sq != pquant->dvel_sq ) {
		pquant->dvel_sq = ( pbig->vel[0] - psmall->vel[0] ) * ( pbig->vel[0] - psmall->vel[0] ) + ( pbig->vel[1] - psmall->vel[1] ) * ( pbig->vel[1] - psmall->vel[1] ) + ( pbig->vel[2] - psmall->vel[2] ) * ( pbig->vel[2] - psmall->vel[2] );
	}

	return pquant->dvel_sq;
}

#define DVEL_SQ collresolve_quant_dvel_sq( conf, pbig, psmall, pquant )

double collresolve_quant_specific_energy( struct collresolve_conf* conf, struct collresolve_body* pbig, struct collresolve_body* psmall, struct collresolve_quant* pquant ) {
	if ( pquant->specific_energy != pquant->specific_energy ) {
		pquant->specific_energy = REDUCED_MASS * ( DVEL_SQ / TOTAL_MASS / 2. + conf->G * ( 1. / TOTAL_RADIUS - 1. / sqrt( DPOS_SQ ) ) );
	}

	return pquant->specific_energy;
}

#define SPECIFIC_ENERGY collresolve_quant_specific_energy( conf, pbig, psmall, pquant )

double collresolve_quant_impact_velocity_sq( struct collresolve_conf* conf, struct collresolve_body* pbig, struct collresolve_body* psmall, struct collresolve_quant* pquant ) {
	if ( pquant->impact_velocity_sq != pquant->impact_velocity_sq ) {
		pquant->impact_velocity_sq = 2. * SPECIFIC_ENERGY * TOTAL_MASS / REDUCED_MASS;
	}

	return pquant->impact_velocity_sq;
}

#define IMPACT_VELOCITY_SQ collresolve_quant_impact_velocity_sq( conf, pbig, psmall, pquant )

double collresolve_quant_escape_velocity_sq( struct collresolve_conf* conf, struct collresolve_body* pbig, struct collresolve_body* psmall, struct collresolve_quant* pquant ) {
	if ( pquant->escape_velocity_sq != pquant->escape_velocity_sq ) {
		pquant->escape_velocity_sq = 2. * conf->G * TOTAL_MASS / TOTAL_RADIUS;
	}

	return pquant->escape_velocity_sq;
}

#define ESCAPE_VELOCITY_SQ collresolve_quant_escape_velocity_sq( conf, pbig, psmall, pquant )

double collresolve_quant_impact_parameter( struct collresolve_conf* conf, struct collresolve_body* pbig, struct collresolve_body* psmall, struct collresolve_quant* pquant ) {
	if ( pquant->impact_parameter != pquant->impact_parameter ) {
		const double dpos_x = psmall->pos[0] - pbig->pos[0];
		const double dpos_y = psmall->pos[1] - pbig->pos[1];
		const double dpos_z = psmall->pos[2] - pbig->pos[2];

		const double dvel_x = psmall->vel[0] - pbig->vel[0];
		const double dvel_y = psmall->vel[1] - pbig->vel[1];
		const double dvel_z = psmall->vel[2] - pbig->vel[2];

		const double h_x = dpos_y * dvel_z - dpos_z * dvel_y;
		const double h_y = dpos_z * dvel_x - dpos_x * dvel_z;
		const double h_z = dpos_x * dvel_y - dpos_y * dvel_x;

		const double h_sq = h_x * h_x + h_y * h_y + h_z * h_z;

		pquant->impact_parameter = sqrt( h_sq / ( TOTAL_RADIUS * TOTAL_RADIUS * IMPACT_VELOCITY_SQ ) );
	}

	return pquant->impact_parameter;
}

#define IMPACR_PARAMETER collresolve_quant_impact_parameter( conf, pbig, psmall, pquant )

/**
 * @}
 */

/**
 * Collision functions
 * @{
 */

double collresolve_specific_energy( struct collresolve_conf* conf, struct collresolve_body big, struct collresolve_body small ) {
	if ( conf == NULL ) {
		return (double)COLLRESOLVE_ERROR_NO_CONF;
	}

	if ( conf->G == 0. ) {
		return (double)COLLRESOLVE_ERROR_INCORRECT_UNIT;
	}

	INIT_QUANT;

	return SPECIFIC_ENERGY;
}

double collresolve_impact_distance( struct collresolve_conf* conf, struct collresolve_body big, struct collresolve_body small ) {
	INIT_QUANT;

	return TOTAL_RADIUS;
}

double collresolve_impact_velocity( struct collresolve_conf* conf, struct collresolve_body big, struct collresolve_body small ) {
	if ( conf == NULL ) {
		return (double)COLLRESOLVE_ERROR_NO_CONF;
	}

	if ( conf->G == 0. ) {
		return (double)COLLRESOLVE_ERROR_INCORRECT_UNIT;
	}

	INIT_QUANT;

	return sqrt( IMPACT_VELOCITY_SQ );
}

double collresolve_escape_velocity( struct collresolve_conf* conf, struct collresolve_body big, struct collresolve_body small ) {
	if ( conf == NULL ) {
		return (double)COLLRESOLVE_ERROR_NO_CONF;
	}

	if ( conf->G == 0. ) {
		return (double)COLLRESOLVE_ERROR_INCORRECT_UNIT;
	}

	INIT_QUANT;

	return sqrt( ESCAPE_VELOCITY_SQ );
}

double collresolve_infinity_velocity( struct collresolve_conf* conf, struct collresolve_body big, struct collresolve_body small ) {
	if ( conf == NULL ) {
		return (double)COLLRESOLVE_ERROR_NO_CONF;
	}

	if ( conf->G == 0. ) {
		return (double)COLLRESOLVE_ERROR_INCORRECT_UNIT;
	}

	INIT_QUANT;

	const double dvel_inf_sq = DVEL_SQ - 2. * conf->G * TOTAL_MASS / sqrt( DPOS_SQ );
	return dvel_inf_sq >= 0. ? sqrt( dvel_inf_sq ) : -sqrt( -dvel_inf_sq );
}

double collresolve_impact_parameter( struct collresolve_conf* conf, struct collresolve_body big, struct collresolve_body small ) {
	if ( conf == NULL ) {
		return (double)COLLRESOLVE_ERROR_NO_CONF;
	}

	INIT_QUANT;

	return IMPACR_PARAMETER;
}

double collresolve_impact_angle( struct collresolve_conf* conf, struct collresolve_body big, struct collresolve_body small ) {
	if ( conf == NULL ) {
		return (double)COLLRESOLVE_ERROR_NO_CONF;
	}

	INIT_QUANT;

	const double impact_para = IMPACR_PARAMETER;

	if ( impact_para > 1. ) {
		return (double)COLLRESOLVE_ERROR_NON_CROSSING;
	}

	return asin( impact_para );
}

/**
 * @}
 */

/**
 * The perfect merging case.
 *
 * Returns a single object, containing all the mass and follows the trajectory of the center of mass.
 */
int collresolve_resolve_merge( struct collresolve_conf* conf, struct collresolve_body big, struct collresolve_body small, int n, struct collresolve_body ret[], int regime ) {
	if ( n < 0 ) {
		return COLLRESOLVE_ERROR_INCORRECT_PARAMETER;
	}

	const double total_mass = big.mass + small.mass;
	const double big_factor = big.mass / total_mass;
	const double small_factor = small.mass / total_mass;

	/* Perfect merging */
	/* Radius of the resulting body is determined assuming it has the same density as the most massive body. */
	ret[0].mass = total_mass;
	ret[0].radius = big.radius * cbrt( total_mass / big.mass );
	ret[0].pos[0] = big.pos[0] * big_factor + small.pos[0] * small_factor;
	ret[0].pos[1] = big.pos[1] * big_factor + small.pos[1] * small_factor;
	ret[0].pos[2] = big.pos[2] * big_factor + small.pos[2] * small_factor;
	ret[0].vel[0] = big.vel[0] * big_factor + small.vel[0] * small_factor;
	ret[0].vel[1] = big.vel[1] * big_factor + small.vel[1] * small_factor;
	ret[0].vel[2] = big.vel[2] * big_factor + small.vel[2] * small_factor;

	/* Nothing else remaining */
	int i;
	for ( i = 1; i <= n; i++ ) {
		ret[i].mass = 0.;
		ret[i].radius = 0.;
		ret[i].pos[0] = 0.;
		ret[i].pos[1] = 0.;
		ret[i].pos[2] = 0.;
		ret[i].vel[0] = 0.;
		ret[i].vel[1] = 0.;
		ret[i].vel[2] = 0.;
	}

	return regime;
}

/**
 * Compute the position and velocity of the second largest remnant in the general case.
 *
 * A two-bodies problem is assumed.
 */
void collresolve_resolve_posvel_slr( struct collresolve_conf* conf, struct collresolve_body* pbig, struct collresolve_body* psmall, struct collresolve_body* plr, struct collresolve_body* pslr ) {
	const double drel = conf->drel == 0. ? 1.01 : conf->drel;

	const double mu = conf->G * ( pbig->mass + psmall->mass );

	const double dpos_x = psmall->pos[0] - pbig->pos[0];
	const double dpos_y = psmall->pos[1] - pbig->pos[1];
	const double dpos_z = psmall->pos[2] - pbig->pos[2];

	const double dvel_x = psmall->vel[0] - pbig->vel[0];
	const double dvel_y = psmall->vel[1] - pbig->vel[1];
	const double dvel_z = psmall->vel[2] - pbig->vel[2];

	const double h_x = dpos_y * dvel_z - dpos_z * dvel_y;
	const double h_y = dpos_z * dvel_x - dpos_x * dvel_z;
	const double h_z = dpos_x * dvel_y - dpos_y * dvel_x;

	const double dpos_sq = dpos_x * dpos_x + dpos_y * dpos_y + dpos_z * dpos_z;
	/* const double dvel_sq = dvel_x * dvel_x + dvel_y * dvel_y + dvel_z * dvel_z; */
	const double h_sq = h_x * h_x + h_y * h_y + h_z * h_z;

	const double dpos = sqrt( dpos_sq );
	const double p = h_sq / mu;

	/* Eccentricity vector, point along the semi-major axis */
	const double e_x = ( dvel_y * h_z - dvel_z * h_y ) / mu - dpos_x / dpos;
	const double e_y = ( dvel_z * h_x - dvel_x * h_z ) / mu - dpos_y / dpos;
	const double e_z = ( dvel_x * h_y - dvel_y * h_x ) / mu - dpos_z / dpos;
	const double e = sqrt( e_x * e_x + e_y * e_y + e_z * e_z );
	const double a = p / ( 1. - e * e );

	const double r = drel * ( plr->radius + pslr->radius );
	const double cf = ( ( a / r ) * ( 1 - e * e ) - 1. ) / e; /* cos(true anomaly) */
	const double sf = sqrt( 1. - cf * cf ); /* sin(true anomaly) */

	const double v = sqrt( mu / p );
	const double cv = -sf;
	const double sv = e + cf;

	/* Direction of the minor axis */
	const double q_x = h_y * e_z - h_z * e_y;
	const double q_y = h_z * e_x - h_x * e_z;
	const double q_z = h_x * e_y - h_y * e_x;
	const double q = sqrt( q_x * q_x + q_y * q_y + q_z * q_z );

	pslr->pos[0] = plr->pos[0] + r * ( e_x / e * cf + q_x / q * sf );
	pslr->pos[1] = plr->pos[1] + r * ( e_y / e * cf + q_y / q * sf );
	pslr->pos[2] = plr->pos[2] + r * ( e_z / e * cf + q_z / q * sf );

	pslr->vel[0] = plr->vel[0] + v * ( e_x / e * cv + q_x / q * sv );
	pslr->vel[1] = plr->vel[1] + v * ( e_y / e * cv + q_y / q * sv );
	pslr->vel[2] = plr->vel[2] + v * ( e_z / e * cv + q_z / q * sv );
}

/**
 * Compute the position and velocity of the second largest remnant in case there is no mass loss.
 *
 * We use momentum and angular momentum conservation to determine its orbit.
 */
void collresolve_resolve_posvel_slr_tot( struct collresolve_conf* conf, struct collresolve_body* pbig, struct collresolve_body* psmall, struct collresolve_body* plr, struct collresolve_body* pslr ) {
	const double drel = conf->drel == 0. ? 1.01 : conf->drel;

	const double dpos_x = psmall->pos[0] - pbig->pos[0];
	const double dpos_y = psmall->pos[1] - pbig->pos[1];
	const double dpos_z = psmall->pos[2] - pbig->pos[2];

	const double dvel_x = psmall->vel[0] - pbig->vel[0];
	const double dvel_y = psmall->vel[1] - pbig->vel[1];
	const double dvel_z = psmall->vel[2] - pbig->vel[2];

	const double h_x = dpos_y * dvel_z - dpos_z * dvel_y;
	const double h_y = dpos_z * dvel_x - dpos_x * dvel_z;
	const double h_z = dpos_x * dvel_y - dpos_y * dvel_x;

	/* Momentum conservation */
	pslr->vel[0] = (pbig->vel[0] * pbig->mass + psmall->vel[0] * psmall->mass - plr->vel[0] * plr->mass) / pslr->mass;
	pslr->vel[1] = (pbig->vel[1] * pbig->mass + psmall->vel[1] * psmall->mass - plr->vel[1] * plr->mass) / pslr->mass;
	pslr->vel[2] = (pbig->vel[2] * pbig->mass + psmall->vel[2] * psmall->mass - plr->vel[2] * plr->mass) / pslr->mass;

	/* The new velocity difference */
	const double avel_x = pslr->vel[0] - plr->vel[0];
	const double avel_y = pslr->vel[1] - plr->vel[1];
	const double avel_z = pslr->vel[2] - plr->vel[2];

	/* Direction where to put the new object, perpendicular to both h and v (i.e. the component that gives angular momentum) */
	const double apos_x = avel_y * h_z - avel_z * h_y;
	const double apos_y = avel_z * h_x - avel_x * h_z;
	const double apos_z = avel_x * h_y - avel_y * h_x;

	const double h_sq = h_x * h_x + h_y * h_y + h_z * h_z;
	const double avel_sq = avel_x * avel_x + avel_y * avel_y + avel_z * avel_z;
	const double apos_sq = apos_x * apos_x + apos_y * apos_y + apos_z * apos_z;

	/* The distance (squared) along apos and avel where to place the object */
	const double rpos_sq = h_sq / avel_sq;
	const double rvel_sq = drel * drel * (plr->radius + pslr->radius) * (plr->radius + pslr->radius) - rpos_sq;

	const double rpos = sqrt( rpos_sq / apos_sq );
	const double rvel = rvel_sq > 0. ? sqrt( rvel_sq / avel_sq ) : 0.; // Just to be sure

	/* Angular momentum conservation */
	pslr->pos[0] = plr->pos[0] + apos_x * rpos + avel_x * rvel;
	pslr->pos[1] = plr->pos[1] + apos_y * rpos + avel_y * rvel;
	pslr->pos[2] = plr->pos[2] + apos_z * rpos + avel_z * rvel;
}

/**
 * The Leinhardt & Stewart (2012) and Stewart & Leinhardt (2012) model.
 *
 * References to equations are to the ones in the former of the above mentionned articles.
 * References to the procedure are to the appendix of the same article.
 */
int collresolve_resolve_ls2012( struct collresolve_conf* conf, struct collresolve_body big, struct collresolve_body small, int n, struct collresolve_body ret[], int flags ) {
	if ( conf == NULL ) {
		return COLLRESOLVE_ERROR_NO_CONF;
	}

	if ( conf->G == 0. ) {
		return COLLRESOLVE_ERROR_INCORRECT_UNIT;
	}

	if ( n < 0 ) {
		return COLLRESOLVE_ERROR_INCORRECT_PARAMETER;
	} else if ( n == 0 ) {
		return collresolve_resolve_merge( conf, big, small, n, ret, COLLRESOLVE_REGIME_MERGE );
	}

	INIT_QUANT;

	/* Common stuff */
	const double m_tot = TOTAL_MASS;
	const double gamma = small.mass / big.mass;
	const double mu = REDUCED_MASS;
	const double r_tot = TOTAL_RADIUS;
	const double Q_R = SPECIFIC_ENERGY;

	/* Impact quantities */
	const double b = IMPACR_PARAMETER;
	const double vel_sq = IMPACT_VELOCITY_SQ;

	if ( b > 1. ) {
		return COLLRESOLVE_ERROR_NON_CROSSING;
	}

	/* Pt 1 */
	const double l = r_tot * ( 1. - b );
	const double alpha = l > 2. * small.radius ? 1 : 0.25 * l * l * ( 3 * small.radius - l ) / ( small.radius * small.radius * small.radius );
	const double m_interact = alpha * small.mass;
	const double m_coll = big.mass + m_interact; /* M' in LS2012 */
	const double big_rho = big.mass / ( 4. / 3. * M_PI * big.radius * big.radius * big.radius );
	const double r_coll = cbrt( 0.75 / M_PI * m_coll / big_rho ); /* R' in LS2012 */
	const double r_coll_tot = cbrt( 0.75 / M_PI * m_tot / big_rho ); /* same as R' but for the total mass */

	/* Pt 2 */
	const double v_p_sq = 2. * m_tot * ( DVEL_SQ / m_tot / 2. + conf->G * ( 1. / r_coll_tot - 1. / sqrt( DPOS_SQ ) ) ); /* TODO: CHECK THIS */
	const double v_esc_p_sq = 2. * conf->G * m_coll / r_coll; /* Eq. (55) but squared */

	if ( v_p_sq < v_esc_p_sq ) {
		/* Perfect merging */
		return collresolve_resolve_merge( conf, big, small, n, ret, COLLRESOLVE_REGIME_MERGE );
	}

	/* Pt 3 */
	const double b_crit = big.radius / r_tot; /* Eq. (6) */
	const int grazing = ( b > b_crit ) ? 1 : 0;

	/* Pt 4 */
	/* a */
	const double rc1 = cbrt( 0.75 / M_PI * m_tot / conf->refDens );

	/* b */
	const double c_star = 1.9; // Fixme, set text bottom right page 9
	const double Q_star_RD_1 = c_star * 0.8 * M_PI * conf->refDens * conf->G * rc1 * rc1; /* Eq. (28) */
	/* const double V_star_1 = sqrt( 6.4 * M_PI * c_star * conf->refDens * conf->G ) * rc1; */ /* Eq. (30) */

	/* c */
	const double mu_alpha = ( alpha * big.mass * small.mass ) / ( big.mass + alpha * small.mass ); // Eq. (12)

	/* d */
	const double mu_bar = 0.36; // Same as for c_star
	const double gamma_term = ( ( gamma + 1. ) * ( gamma + 1. ) / ( 4. * gamma ) );
	const double Q_star_RD = Q_star_RD_1 * pow( gamma_term,  2. / ( 3. * mu_bar ) - 1. ); /* Eq. (23) */
	/* const double V_star = V_star_1 * pow( gamma_term, 1. / ( 3. * mu_bar ) ); */ /* Eq. (22) */

	/* e */
	const double Q_prime_star_RD = Q_star_RD * pow( mu / mu_alpha, 2. - 3. * mu_bar / 2. ); /* Eq. (15) */
	/* const double V_prime_star = V_star * sqrt( 2. * Q_prime_star_RD * m_tot / mu ); */ /* Eq. (16) */

	/* Pt 5 */
	const double Q_R_erosion = Q_prime_star_RD * 2. * gamma / ( gamma + 1. ); // From eq. (5) and a sheet of paper
	/* const double V_erosion_sq = sqrt( 2. * Q_R_erosion * m_tot / mu ); */ /* From eq. (1) and a sheet of paper, but squared */

	/* Pt 7 */
	const double Q_R_supercat = 1.8 * Q_prime_star_RD; /* From eq. (5) and a sheet of paper */
	//const double V_supercat_sq = 2. * Q_R_supercat * m_tot / mu; */ /* From eq. (1) and a sheet of paper, but squared */

	const double p_big = b > 0.7 ? 1. : 1. - ( b / 0.7 ) * ( small.mass / m_tot );
	const double p_small = b > 0.7 ? 0. : ( b / 0.7 ) * ( small.mass / m_tot );

	/* Graze-and-Merge regime from Kokubo & Ida (2012); only used in the Stewart & Leinhardt (2012) model */
	const double xi = ( big.mass - small.mass ) / m_tot;
	const double b_para_gnm = pow( 1. - b, 2.5 );
	const double vel_hr = 2.43 * xi * xi * b_para_gnm - 0.0408 * xi * xi + 1.86 * b_para_gnm + 1.08;

	if ( Q_R > Q_R_supercat ) {
		/* Supercatastrophic; pt 9 */
		ret[0].mass = 0.1 * m_tot * pow( Q_R / Q_prime_star_RD / 1.8, -1.5 );
		ret[0].radius = cbrt( 0.75 / M_PI * ret[0].mass / big_rho );
		ret[0].pos[0] = big.pos[0] * p_big + small.pos[0] * p_small;
		ret[0].pos[1] = big.pos[1] * p_big + small.pos[1] * p_small;
		ret[0].pos[2] = big.pos[2] * p_big + small.pos[2] * p_small;
		ret[0].vel[0] = big.vel[0] * p_big + small.vel[0] * p_small;
		ret[0].vel[1] = big.vel[1] * p_big + small.vel[1] * p_small;
		ret[0].vel[2] = big.vel[2] * p_big + small.vel[2] * p_small;
		int i;
		for ( i = 1; i < n; i++ ) {
			ret[i].mass = 0.;
			ret[i].radius = 0.;
			ret[i].pos[0] = 0.;
			ret[i].pos[1] = 0.;
			ret[i].pos[2] = 0.;
			ret[i].vel[0] = 0.;
			ret[i].vel[1] = 0.;
			ret[i].vel[2] = 0.;
		}
		ret[n].mass = m_tot - ret[0].mass;
		ret[n].radius = 0.;
		ret[n].pos[0] = 0.;
		ret[n].pos[1] = 0.;
		ret[n].pos[2] = 0.;
		ret[n].vel[0] = 0.;
		ret[n].vel[1] = 0.;
		ret[n].vel[2] = 0.;

		return COLLRESOLVE_REGIME_SUPERCATASTROPHIC;
	} else if ( Q_R > Q_R_erosion || !grazing ) {
		/* Disruption; pt 8 */
		const double m_lr = m_tot * ( 1. - Q_R / ( 2. * Q_prime_star_RD ) );
		if ( m_lr >= m_tot ) {
			return collresolve_resolve_merge( conf, big, small, n, ret, COLLRESOLVE_REGIME_MERGE );
		} else {
			ret[0].mass = m_lr;
			ret[0].radius = cbrt( 0.75 / M_PI * m_lr / big_rho );
			ret[0].pos[0] = big.pos[0] * p_big + small.pos[0] * p_small;
			ret[0].pos[1] = big.pos[1] * p_big + small.pos[1] * p_small;
			ret[0].pos[2] = big.pos[2] * p_big + small.pos[2] * p_small;
			ret[0].vel[0] = big.vel[0] * p_big + small.vel[0] * p_small;
			ret[0].vel[1] = big.vel[1] * p_big + small.vel[1] * p_small;
			ret[0].vel[2] = big.vel[2] * p_big + small.vel[2] * p_small;
			double m_r = m_lr;
			if ( n > 1 ) {
				double m_slr = ( 0.15 * ( m_lr / m_tot ) / ( 2.85 * 2 ) ) * m_tot; /* Eq. (37) */
				int all_mass = 0;
				if ( m_lr + m_slr > m_tot ) {
					m_slr = m_tot - m_lr;
					all_mass = 1;
				}
				m_r += m_slr;
				ret[1].mass = m_slr;
				ret[1].radius = small.radius * cbrt( m_slr / small.mass ); /* Keeps density of the impactor */
				if ( all_mass ) {
					collresolve_resolve_posvel_slr_tot( conf, pbig, psmall, ret, ret + 1 );
				} else {
					collresolve_resolve_posvel_slr( conf, pbig, psmall, ret, ret + 1 );
				}
				/* Todo: Implement eq. (31) */
				int i;
				for ( i = 2; i < n; i++ ) {
					ret[i].mass = 0.;
					ret[i].radius = 0.;
					ret[i].pos[0] = 0.;
					ret[i].pos[1] = 0.;
					ret[i].pos[2] = 0.;
					ret[i].vel[0] = 0.;
					ret[i].vel[1] = 0.;
					ret[i].vel[2] = 0.;
				}
			}
			ret[n].mass = m_tot - m_r;
			ret[n].radius = 0.;
			ret[n].pos[0] = 0.;
			ret[n].pos[1] = 0.;
			ret[n].pos[2] = 0.;
			ret[n].vel[0] = 0.;
			ret[n].vel[1] = 0.;
			ret[n].vel[2] = 0.;
		}

		return COLLRESOLVE_REGIME_DISRUPTION;
	} else if ( sqrt( vel_sq / ESCAPE_VELOCITY_SQ ) < vel_hr && ( flags & 1 ) ) {
		/* Graze and Merge; same outcome as perfect merging */
		return collresolve_resolve_merge( conf, big, small, n, ret, COLLRESOLVE_REGIME_GRAZE_AND_MERGE );
	} else {
		/* Hit-and-run; pt 6 */
		ret[0].mass = big.mass;
		ret[0].radius = big.radius;
		ret[0].pos[0] = big.pos[0];
		ret[0].pos[1] = big.pos[1];
		ret[0].pos[2] = big.pos[2];
		ret[0].vel[0] = big.vel[0];
		ret[0].vel[1] = big.vel[1];
		ret[0].vel[2] = big.vel[2];
		double m_slr = 0.;
		if ( n > 1 ) {
			const double phi_rev = 2. * acos( ( l - small.radius ) / small.radius );
			const double A_interact = small.radius * small.radius * ( M_PI - ( phi_rev - sin( phi_rev ) ) / 2. ); /* Eq. (46) */
			const double L_interact = 2. * sqrt( big.radius * big.radius - ( big.radius - l / 2. ) * ( big.radius - l / 2. ) ); /* Eq. (47) */
			const double M_interact = A_interact * L_interact; /* Eq. (48) */

			const double mass1_rev = small.mass;
			const double mass2_rev = M_interact;
			const double rho1_rev = small.mass / ( 4. / 3. * M_PI * small.radius * small.radius * small.radius );

			const double m_tot_rev = mass1_rev + mass2_rev;
			const double mu_rev = mass1_rev * mass2_rev / m_tot_rev; /* Eq. (49) */
			const double gamma_rev = mass2_rev / mass1_rev; /* Eq. (50) */
			const double Q_R_rev = mu_rev * vel_sq / ( 2. * m_tot_rev );

			/* a */
			const double rc1_rev = cbrt( 0.75 / M_PI * m_tot_rev / conf->refDens );

			/* b */
			const double Q_star_RD_1_rev = c_star * 0.8 * M_PI * rho1_rev * conf->G * rc1_rev * rc1_rev; /* Eq. (28) */
			/* const double V_star_1_rev = sqrt( 6.4 * M_PI * c_star * rho1_rev * conf->G ) * rc1_rev; */ /* Eq. (30) */

			/* c */
			const double mu_alpha_rev = ( alpha * big.mass * small.mass ) / ( big.mass + alpha * small.mass ); /* Eq. (12) */

			/* d */
			const double gamma_term_rev = ( ( gamma_rev + 1. ) * ( gamma_rev + 1. ) / ( 4. * gamma_rev ) );
			const double Q_star_RD_rev = Q_star_RD_1_rev * pow( gamma_term_rev,  2. / ( 3. * mu_bar ) - 1. ); /* Eq. (23) */
			/* const double V_star_rev = V_star_1_rev * pow( gamma_term_rev, 1. / ( 3. * mu_bar ) ); */ /* Eq. (22) */

			/* e */
			const double Q_prime_star_RD_rev = Q_star_RD_rev * pow( mu_rev / mu_alpha_rev, 2. - 3. * mu_bar / 2. ); /* Eq. (15) */
			/* const double V_prime_star_rev = V_star_rev * sqrt( 2. * Q_prime_star_RD_rev * m_tot_rev / mu_rev ); */ /* Eq. (16) */

			const double Q_R_supercat_rev = Q_prime_star_RD_rev * 2. * ( 0.1 * gamma_rev - 0.9 ) / ( gamma_rev + 1. ); /* From eq. (5) and a sheet of paper */
			/* const double V_supercat_rev_sq = 2. * Q_R_supercat_rev * m_tot_rev / mu_rev; */ /* From eq. (1) and a sheet of paper, but squared */

			if ( Q_R_rev > Q_R_supercat_rev ) {
				m_slr = 0.1 * m_tot_rev * pow( Q_R_rev / Q_prime_star_RD_rev / 1.8, -1.5 );
			} else {
				m_slr = m_tot_rev * ( 1. - Q_R_rev / ( 2. * Q_prime_star_RD_rev ) );
			}

			/* Mass conservation */
			if ( m_slr > small.mass ) {
				m_slr = small.mass;
			}

			ret[1].mass = m_slr;
			ret[1].radius = small.radius * cbrt( m_slr / small.mass ); /* Keeps density of the impactor */
			collresolve_resolve_posvel_slr( conf, pbig, psmall, ret, ret + 1 );

			int i;
			for ( i = 2; i < n; i++ ) {
				ret[i].mass = 0.;
				ret[i].radius = 0.;
				ret[i].pos[0] = 0.;
				ret[i].pos[1] = 0.;
				ret[i].pos[2] = 0.;
				ret[i].vel[0] = 0.;
				ret[i].vel[1] = 0.;
				ret[i].vel[2] = 0.;
			}
		}
		ret[n].mass = small.mass - m_slr;
		ret[n].radius = 0.;
		ret[n].pos[0] = 0.;
		ret[n].pos[1] = 0.;
		ret[n].pos[2] = 0.;
		ret[n].vel[0] = 0.;
		ret[n].vel[1] = 0.;
		ret[n].vel[2] = 0.;

		return COLLRESOLVE_REGIME_HIT_AND_RUN;
	}
}

int collresolve_resolve_c2019( struct collresolve_conf* conf, struct collresolve_body big, struct collresolve_body small, int n, struct collresolve_body ret[] ) {
	/* Inputs are: Target mass [Earth mass]; Projectile-to-target mass ratio; Impact angle [degree]; Impact-to-escape velocity ratio. */
	/* Regimes are: -1: Erosion, 0: Accretion, 1: Hit-and-run. */

	if ( conf == NULL ) {
		return COLLRESOLVE_ERROR_NO_CONF;
	}

	if ( conf->G == 0. ) {
		return COLLRESOLVE_ERROR_INCORRECT_UNIT;
	}

	/* Check for specific values of n */
	if ( n < 0 ) {
		return COLLRESOLVE_ERROR_INCORRECT_PARAMETER;
	} else if ( n == 0 ) {
		collresolve_resolve_merge( conf, big, small, n, ret, COLLRESOLVE_REGIME_MERGE );
		collresolve_body_radius( conf, &( ret[0] ) );
		return COLLRESOLVE_REGIME_MERGE;
	}

	INIT_QUANT;

	double impact_para = IMPACR_PARAMETER;
	if ( impact_para > 1. ) {
		return COLLRESOLVE_ERROR_NON_CROSSING;
	}

	double params[4];

	params[0] = log10( big.mass / conf->mEarth );
	params[1] = small.mass / big.mass;
	params[2] = asin( impact_para ) * 180. / M_PI;
	params[3] = sqrt( IMPACT_VELOCITY_SQ / ESCAPE_VELOCITY_SQ );

	int regime;
	double score[3];
	collision_classifier( params, &regime, score );

	double accs[2];
	accretion_efficiency( params, accs );

	double orb[3];
	orbital_hnr( params, orb );

	int grazing = 0;
	if ( regime == 1 ) {
		/* For now, only hit and run is grazing, i.e. has a second remnant. */
		if ( orb[0] > 0. ) {
			grazing = 1;
		} else {
			/* This is a consistency check to be sure that the orbit is unbound. */
			grazing = 0;
		}
	} else {
		grazing = 0;
	}

	const double total_mass = big.mass + small.mass;
	const double big_factor = big.mass / total_mass;
	const double small_factor = small.mass / total_mass;

	if ( grazing && n > 1 ) {
		/* Two remnants */

		/* First let's compute the resulting mass and radius */
		double acclr = accs[0];
		double accsr = accs[1];
		double acctr = 0. - acclr - accsr;

		/* Mass conservation */
		if ( acclr > 1. ) {
			acclr = 1.;
			accsr = -1.;
			acctr = 0.;
		} else if ( acclr + accsr > 0. ) {
			accsr = -acclr;
			acctr = 0.;
		}

		ret[0].mass = big.mass + small.mass * acclr;
		collresolve_body_radius( conf, &( ret[0] ) );

		ret[1].mass = small.mass * ( 1. + accsr );
		collresolve_body_radius( conf, &( ret[1] ) );

		/* Now, we do quite a bit of orbital mechanics. */
		const double drel = conf->drel == 0. ? 1.01 : conf->drel;

		const double mu = conf->G * total_mass;

		const double dpos_x = small.pos[0] - big.pos[0];
		const double dpos_y = small.pos[1] - big.pos[1];
		const double dpos_z = small.pos[2] - big.pos[2];

		const double dvel_x = small.vel[0] - big.vel[0];
		const double dvel_y = small.vel[1] - big.vel[1];
		const double dvel_z = small.vel[2] - big.vel[2];

		const double h_x = dpos_y * dvel_z - dpos_z * dvel_y;
		const double h_y = dpos_z * dvel_x - dpos_x * dvel_z;
		const double h_z = dpos_x * dvel_y - dpos_y * dvel_x;

		const double dpos_sq = dpos_x * dpos_x + dpos_y * dpos_y + dpos_z * dpos_z;
		const double dpos = sqrt( dpos_sq );
		/* const double dvel_sq = dvel_x * dvel_x + dvel_y * dvel_y + dvel_z * dvel_z; */
		const double h_sq = h_x * h_x + h_y * h_y + h_z * h_z;
		const double h = sqrt( h_sq );

		/* Eccentricity vector, point along the semi-major axis */
		const double e_x = ( dvel_y * h_z - dvel_z * h_y ) / mu - dpos_x / dpos;
		const double e_y = ( dvel_z * h_x - dvel_x * h_z ) / mu - dpos_y / dpos;
		const double e_z = ( dvel_x * h_y - dvel_y * h_x ) / mu - dpos_z / dpos;
		const double e = sqrt( e_x * e_x + e_y * e_y + e_z * e_z );

		/* Rotate the eccentricity vector, by the shift of the argument of pericenter, about the direction given by h. */
		const double peri_rot = orb[2];
		const double srot = sin( peri_rot );
		const double crot = cos( peri_rot );
		const double omcrot = 1. - crot;

		const double axis_x = h_x / h;
		const double axis_y = h_y / h;
		const double axis_z = h_z / h;

		const double maj_x = ( e_x / e ) * ( crot + axis_x * axis_x * omcrot )          + ( e_y / e ) * ( axis_x * axis_y * omcrot - axis_z * srot ) + ( e_z / e ) * ( axis_x * axis_z * omcrot + axis_y * srot );
		const double maj_y = ( e_x / e ) * ( axis_y * axis_x * omcrot + axis_z * srot ) + ( e_y / e ) * ( crot + axis_y * axis_y * omcrot )          + ( e_z / e ) * ( axis_y * axis_z * omcrot - axis_x * srot );
		const double maj_z = ( e_x / e ) * ( axis_z * axis_x * omcrot - axis_y * srot ) + ( e_y / e ) * ( axis_z * axis_y * omcrot + axis_x * srot ) + ( e_z / e ) * ( crot + axis_z * axis_z * omcrot );
		const double m = sqrt( maj_x * maj_x + maj_y * maj_y + maj_z * maj_z );

		/* Direction of the minor axis */
		const double q_x = axis_y * maj_z - axis_z * maj_y;
		const double q_y = axis_z * maj_x - axis_x * maj_z;
		const double q_z = axis_x * maj_y - axis_y * maj_x;
		const double q = sqrt( q_x * q_x + q_y * q_y + q_z * q_z );

		/* Orbital characteristics of the remnants */
		const double eorb = orb[0];
		const double b = orb[1];

		const double p = 2. * b * b * ( 1. + eorb );
		const double e_sq = 1. + 4. * b * b * eorb * ( 1. + eorb );
		const double e_af = e_sq < 0. ? 0. : sqrt( e_sq );

		const double cf = ( p / drel - 1. ) / e_af; /* cos(true anomaly) */
		const double sf = sqrt( 1. - cf * cf ); /* sin(true anomaly) */

		const double cv = -sf;
		const double sv = e_af + cf;

		/* Characteristics of the bodies */
		const double mu_af = conf->G * ( ret[0].mass + ret[1].mass );
		const double r_ref = ret[0].radius + ret[1].radius;

		/* Relative position */
		const double r = drel * r_ref;
		const double rpos_x = r * ( maj_x / m * cf + q_x / q * sf );
		const double rpos_y = r * ( maj_y / m * cf + q_y / q * sf );
		const double rpos_z = r * ( maj_z / m * cf + q_z / q * sf );

		/* Relative velocity */
		const double v = sqrt( mu_af / ( r_ref * p ) );
		const double rvel_x = v * ( maj_x / m * cv + q_x / q * sv );
		const double rvel_y = v * ( maj_y / m * cv + q_y / q * sv );
		const double rvel_z = v * ( maj_z / m * cv + q_z / q * sv );

		/* Set the positions and velocities */
		const double lr_factor = -ret[1].mass / ( ret[0].mass + ret[1].mass );
		const double sr_factor =  ret[0].mass / ( ret[0].mass + ret[1].mass );

		ret[0].pos[0] = big.pos[0] * big_factor + small.pos[0] * small_factor + lr_factor * rpos_x;
		ret[0].pos[1] = big.pos[1] * big_factor + small.pos[1] * small_factor + lr_factor * rpos_y;
		ret[0].pos[2] = big.pos[2] * big_factor + small.pos[2] * small_factor + lr_factor * rpos_z;
		ret[0].vel[0] = big.vel[0] * big_factor + small.vel[0] * small_factor + lr_factor * rvel_x;
		ret[0].vel[1] = big.vel[1] * big_factor + small.vel[1] * small_factor + lr_factor * rvel_y;
		ret[0].vel[2] = big.vel[2] * big_factor + small.vel[2] * small_factor + lr_factor * rvel_z;

		ret[1].pos[0] = big.pos[0] * big_factor + small.pos[0] * small_factor + sr_factor * rpos_x;
		ret[1].pos[1] = big.pos[1] * big_factor + small.pos[1] * small_factor + sr_factor * rpos_y;
		ret[1].pos[2] = big.pos[2] * big_factor + small.pos[2] * small_factor + sr_factor * rpos_z;
		ret[1].vel[0] = big.vel[0] * big_factor + small.vel[0] * small_factor + sr_factor * rvel_x;
		ret[1].vel[1] = big.vel[1] * big_factor + small.vel[1] * small_factor + sr_factor * rvel_y;
		ret[1].vel[2] = big.vel[2] * big_factor + small.vel[2] * small_factor + sr_factor * rvel_z;

		int i;
		for ( i = 2; i < n; i++ ) {
			ret[i].mass = 0.;
			ret[i].radius = 0.;
			ret[i].pos[0] = 0.;
			ret[i].pos[1] = 0.;
			ret[i].pos[2] = 0.;
			ret[i].vel[0] = 0.;
			ret[i].vel[1] = 0.;
			ret[i].vel[2] = 0.;
		}

		ret[n].mass = small.mass * acctr;
		ret[n].radius = 0.;
		ret[n].pos[0] = 0.;
		ret[n].pos[1] = 0.;
		ret[n].pos[2] = 0.;
		ret[n].vel[0] = 0.;
		ret[n].vel[1] = 0.;
		ret[n].vel[2] = 0.;
	} else {
		/* Non-grazing collision; only a single remnant */
		double acc = accs[0];

		double rem_mass = big.mass + small.mass * acc;
		if ( rem_mass > total_mass ) {
			rem_mass = total_mass;
		}

		ret[0].mass = rem_mass;
		collresolve_body_radius( conf, &( ret[0] ) );
		ret[0].pos[0] = big.pos[0] * big_factor + small.pos[0] * small_factor;
		ret[0].pos[1] = big.pos[1] * big_factor + small.pos[1] * small_factor;
		ret[0].pos[2] = big.pos[2] * big_factor + small.pos[2] * small_factor;
		ret[0].vel[0] = big.vel[0] * big_factor + small.vel[0] * small_factor;
		ret[0].vel[1] = big.vel[1] * big_factor + small.vel[1] * small_factor;
		ret[0].vel[2] = big.vel[2] * big_factor + small.vel[2] * small_factor;

		int i;
		for ( i = 1; i < n; i++ ) {
			ret[i].mass = 0.;
			ret[i].radius = 0.;
			ret[i].pos[0] = 0.;
			ret[i].pos[1] = 0.;
			ret[i].pos[2] = 0.;
			ret[i].vel[0] = 0.;
			ret[i].vel[1] = 0.;
			ret[i].vel[2] = 0.;
		}

		ret[n].mass = total_mass - rem_mass;
		ret[n].radius = 0.;
		ret[n].pos[0] = 0.;
		ret[n].pos[1] = 0.;
		ret[n].pos[2] = 0.;
		ret[n].vel[0] = 0.;
		ret[n].vel[1] = 0.;
		ret[n].vel[2] = 0.;
	}

	/* Type of collision */
	if ( regime == 1 ) {
		if ( orb[0] > 0. ) {
			return COLLRESOLVE_REGIME_HIT_AND_RUN;
		} else {
			/* This is a special value to indicate the "inconsistent" outcome. */
			return COLLRESOLVE_REGIME_GRAZE_AND_MERGE;
		}
	} else if ( regime == 0 ) {
		return COLLRESOLVE_REGIME_MERGE;
	} else if ( regime == -1 ) {
		return COLLRESOLVE_REGIME_DISRUPTION;
	} else {
		return 0;
	}
}

int collresolve_resolve( struct collresolve_conf* conf, struct collresolve_body big, struct collresolve_body small, int n, struct collresolve_body ret[] ) {
	if ( conf == NULL ) {
		return COLLRESOLVE_ERROR_NO_CONF;
	} else if ( conf->model == COLLRESOLVE_MODEL_PERFECT_MERGE ) {
		return collresolve_resolve_merge( conf, big, small, n, ret, COLLRESOLVE_REGIME_MERGE );
	} else if ( conf->model == COLLRESOLVE_MODEL_LS2012 ) {
		return collresolve_resolve_ls2012( conf, big, small, n, ret, 0 );
	} else if ( conf->model == COLLRESOLVE_MODEL_SL2012 ) {
		return collresolve_resolve_ls2012( conf, big, small, n, ret, 1 );
	} else if ( conf->model == COLLRESOLVE_MODEL_C2019 ) {
		return collresolve_resolve_c2019( conf, big, small, n, ret );
	} else {
		return COLLRESOLVE_ERROR_INCORRECT_MODEL;
	}
}
