/**
 * FORTRAN bindings for the collresolve library.
 *
 * Copyright (c) 2016-2017 University of Bern, Switzerland
 * Copyright (c) 2018-2019 Arizona Board of Regents
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 * @file
 * @author Alexandre Emsenhuber
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdlib.h> /* For malloc()/free() */
#include <math.h> /* For M_PI/cbrt() */

#include "collresolve.h"

/**
 * FORTRAN bindings for the the collresolve_resolve() function.
 *
 * Configuration is the following:
 * - Collision model is provided as first parameter
 * - Unit system is fixed to Mercury's one.
 * - Relative distance after collision factor is fixed to 1.25.
 *
 * Interface:
 * INTEGER, INTENT(IN) :: model                            ! collision model to apply
 * DOUBLE PRECISION, INTENT(IN) :: m1                      ! mass of the target
 * DOUBLE PRECISION, INTENT(IN) :: m2                      ! mass of the impactor
 * DOUBLE PRECISION, INTENT(IN) :: r1                      ! radius of the target
 * DOUBLE PRECISION, INTENT(IN) :: r2                      ! radius of the impactor
 * DOUBLE PRECISION, DIMENSION(3), INTENT(IN) :: p1        ! position of the target
 * DOUBLE PRECISION, DIMENSION(3), INTENT(IN) :: p2        ! position of the impactor
 * DOUBLE PRECISION, DIMENSION(3), INTENT(IN) :: v1        ! velocity of the target
 * DOUBLE PRECISION, DIMENSION(3), INTENT(IN) :: v2        ! velocity of the impactor
 * INTEGER, INTENT(IN) :: n                                ! number of bodies to return
 * DOUBLE PRECISION, DIMENSION(n+1), INTENT(OUT) :: mres   ! mass of the resulting bodies
 * DOUBLE PRECISION, DIMENSION(n+1), INTENT(OUT) :: rres   ! radius of the resulting bodies
 * DOUBLE PRECISION, DIMENSION(3,n+1), INTENT(OUT) :: pres ! position of the resulting bodies
 * DOUBLE PRECISION, DIMENSION(3,n+1), INTENT(OUT) :: vres ! velocity of the resulting bodies
 */
int collresolve_resolve_( int* model, double* m1, double* m2, double* r1, double* r2, double* p1, double* p2, double* v1, double* v2, int* n,
	double* mres, double* rres, double* pres, double* vres
) {
	struct collresolve_conf* conf = collresolve_conf_new();
	collresolve_conf_unit_merc( conf );
	collresolve_conf_model( conf, *model );
	collresolve_conf_sep_after( conf, 1.25 ); /* There's a factor 1.2 hardcoded in Mercury. */

	struct collresolve_body big, small;
	big.mass = *m1;
	small.mass = *m2;
	big.radius = *r1;
	small.radius = *r2;
	for ( int i = 0; i < 3; i++ ) {
		big.pos[i] = p1[i];
		small.pos[i] = p2[i];
		big.vel[i] = v1[i];
		small.vel[i] = v2[i];
	}

	struct collresolve_body* ret = malloc( sizeof( struct collresolve_body ) * ( *n + 1 ) );

	int res = collresolve_resolve( conf, big, small, *n, ret );

	for ( int i = 0; i <= *n; i++ ) {
		mres[i] = ret[i].mass;
		rres[i] = ret[i].radius;
		for ( int j = 0; j < 3; j++ ) {
			pres[3*i+j] = ret[i].pos[j];
			vres[3*i+j] = ret[i].vel[j];
		}
	}

	collresolve_conf_free( conf );

	free( ret );

	return res;
}


/**
 * FORTRAN binding to retrieve the radius of a body of a given mass.
 *
 * Configuration is the following:
 * - Collision model is provided as first parameter
 * - Unit system is fixed to Mercury's one.
 *
 * Interface:
 * INTEGER, INTENT(IN) :: model          ! collision model to apply
 * DOUBLE PRECISION, INTENT(IN) :: mass  ! Mass of the body
 */
double collresolve_radius_( int* model, double* mass ) {
	struct collresolve_conf* conf = collresolve_conf_new();
	collresolve_conf_unit_merc( conf );
	collresolve_conf_model( conf, *model );

	double dens = collresolve_bulk_density( conf, *mass );

	collresolve_conf_free( conf );

	return cbrt( 3. / 4. / M_PI * (*mass) / dens );
}
