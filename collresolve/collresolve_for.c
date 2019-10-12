#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

/* #include <stdio.h> */
#include <stdlib.h> /* For malloc()/free() */

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
	/*
	printf( "collresolve_resolve_:\n" );
	printf( "model=%i\n", *model );
	printf( "m1=%g m2=%g\n", *m1, *m2 );
	printf( "r1=%g r2=%g\n", *r1, *r2 );
	printf( "p1=(%g, %g, %g)\n", p1[0], p1[1], p1[2] );
	printf( "p2=(%g, %g, %g)\n", p2[0], p2[1], p2[2] );
	printf( "v1=(%g, %g, %g)\n", v1[0], v1[1], v1[2] );
	printf( "v2=(%g, %g, %g)\n", v2[0], v2[1], v2[2] );
	printf( "n=%d\n", *n );
	*/

	struct collresolve_conf* conf = collresolve_conf_new();
	collresolve_conf_unit_merc( conf );
	collresolve_conf_model( conf, *model );
	collresolve_conf_sep_after( conf, 1.25 ); /* There's a factor 1.2 hardcoded in Mercury. */

	struct collresolve_body big, small;
	big.mass = *m1;
	small.mass = *m2;
	big.radius = *r1;
	small.radius = *r2;
	int i;
	for ( i = 0; i < 3; i++ ) {
		big.pos[i] = p1[i];
		small.pos[i] = p2[i];
		big.vel[i] = v1[i];
		small.vel[i] = v2[i];
	}

	struct collresolve_body* ret = malloc( sizeof( struct collresolve_body ) * ( *n + 1 ) );

	int res = collresolve_resolve( conf, big, small, *n, ret );
	for ( i = 0; i <= *n; i++ ) {
		mres[i] = ret[i].mass;
		rres[i] = ret[i].radius;
		int j;
		for ( j = 0; j < 3; j++ ) {
			pres[3*i+j] = ret[i].pos[j];
			vres[3*i+j] = ret[i].vel[j];
		}
	}

	collresolve_conf_free( conf );

	free( ret );

	return res;
}
