/**
 * The library's interface as a Python module.
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
 * If you use this library in a scientific work that lead to publication,
 * we would like you to acknowledge the following article:
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

#include "Python.h"
#include "structmember.h"

#include <stddef.h>

#include "collresolve.h"

/**
 * Common stuff
 */

static PyObject* collresolve_Error;

/**
 * Body object
 */

struct collresolve_BodyObject {
	PyObject_HEAD
	struct collresolve_body data;
};

static PyObject* collresolve_BodyObject_repr( struct collresolve_BodyObject* obj ) {
	PyObject* pymass = PyFloat_FromDouble( obj->data.mass );
	PyObject* pyradius = PyFloat_FromDouble( obj->data.radius );
	PyObject* pyposx = PyFloat_FromDouble( obj->data.pos[0] );
	PyObject* pyposy = PyFloat_FromDouble( obj->data.pos[1] );
	PyObject* pyposz = PyFloat_FromDouble( obj->data.pos[2] );
	PyObject* pyvelx = PyFloat_FromDouble( obj->data.vel[0] );
	PyObject* pyvely = PyFloat_FromDouble( obj->data.vel[1] );
	PyObject* pyvelz = PyFloat_FromDouble( obj->data.vel[2] );

	PyObject* ret = PyUnicode_FromFormat( "collresolve.Body(mass=%R, radius=%R, pos_x=%R, pos_y=%R, pos_z=%R, vel_x=%R, vel_y=%R, vel_z=%R)", pymass, pyradius, pyposx, pyposy, pyposz, pyvelx, pyvely, pyvelz );

	Py_DECREF( pymass );
	Py_DECREF( pyradius );
	Py_DECREF( pyposx );
	Py_DECREF( pyposy );
	Py_DECREF( pyposz );
	Py_DECREF( pyvelx );
	Py_DECREF( pyvely );
	Py_DECREF( pyvelz );

	return ret;
}

static int collresolve_BodyObject_init( PyObject* self, PyObject* args, PyObject* kwds ) {
	static char *kwlist[] = { "mass", "radius", "pos_x", "pos_y", "pos_z", "vel_x", "vel_y", "vel_z", NULL };

	struct collresolve_BodyObject* pybody = (struct collresolve_BodyObject*) self;

	if ( !PyArg_ParseTupleAndKeywords( args, kwds, "|dddddddd", kwlist, &pybody->data.mass, &pybody->data.radius, &pybody->data.pos[0], &pybody->data.pos[1], &pybody->data.pos[2], &pybody->data.vel[0], &pybody->data.vel[1], &pybody->data.vel[2] ) ) {
		return -1;
	}

	return 0;
}

static PyMemberDef collresolve_BodyObject_members[] = {
	{ "mass", T_DOUBLE, offsetof(struct collresolve_BodyObject, data.mass), 0, "body's mass" },
	{ "radius", T_DOUBLE, offsetof(struct collresolve_BodyObject, data.radius), 0, "body's radius" },
	{ "pos_x", T_DOUBLE, offsetof(struct collresolve_BodyObject, data.pos[0]), 0, "body's position" },
	{ "pos_y", T_DOUBLE, offsetof(struct collresolve_BodyObject, data.pos[1]), 0, "body's position" },
	{ "pos_z", T_DOUBLE, offsetof(struct collresolve_BodyObject, data.pos[2]), 0, "body's position" },
	{ "vel_x", T_DOUBLE, offsetof(struct collresolve_BodyObject, data.vel[0]), 0, "body's velocity" },
	{ "vel_y", T_DOUBLE, offsetof(struct collresolve_BodyObject, data.vel[1]), 0, "body's velocity" },
	{ "vel_z", T_DOUBLE, offsetof(struct collresolve_BodyObject, data.vel[2]), 0, "body's velocity" },
	{ NULL }  /* Sentinel */
};

/**
 * Object type
 */
static PyTypeObject collresolve_BodyType = {
	PyVarObject_HEAD_INIT( NULL, 0 )
	"collresolve.Body",     /* tp_name */
	sizeof( struct collresolve_BodyObject ), /* tp_basicsize */
	0,                      /* tp_itemsize */
	0,                      /* tp_dealloc */
	0,                      /* tp_print */
	0,                      /* tp_getattr */
	0,                      /* tp_setattr */
	0,                      /* tp_reserved */
	(reprfunc) collresolve_BodyObject_repr, /* tp_repr */
	0,                      /* tp_as_number */
	0,                      /* tp_as_sequence */
	0,                      /* tp_as_mapping */
	0,                      /* tp_hash  */
	0,                      /* tp_call */
	0,                      /* tp_str */
	0,                      /* tp_getattro */
	0,                      /* tp_setattro */
	0,                      /* tp_as_buffer */
	Py_TPFLAGS_DEFAULT,     /* tp_flags */
	"Collision body object", /* tp_doc */
	0,                      /* tp_traverse */
	0,                      /* tp_clear */
	0,                      /* tp_richcompare */
	0,                      /* tp_weaklistoffset */
	0,                      /* tp_iter */
	0,                      /* tp_iternext */
	0,                      /* tp_methods */
	collresolve_BodyObject_members, /* tp_members */
	0,                      /* tp_getset */
	0,                      /* tp_base */
	0,                      /* tp_dict */
	0,                      /* tp_descr_get */
	0,                      /* tp_descr_set */
	0,                      /* tp_dictoffset */
	(initproc) collresolve_BodyObject_init, /* tp_init */
	0,                      /* tp_alloc */
	0,                      /* tp_new */
};

/**
 * Configuration object
 */

struct collresolve_ConfObject {
	PyObject_HEAD
	struct collresolve_conf* data;
};

/**
 * Object methods
 */

static void python_collresolve_conf_dealloc( struct collresolve_ConfObject* obj ) {
	collresolve_conf_free( obj->data );
	Py_TYPE( obj )->tp_free( obj );
}

/**
 * Object type
 */
static PyTypeObject collresolve_ConfType = {
	PyVarObject_HEAD_INIT( NULL, 0 )
	"collresolve.Conf",     /* tp_name */
	sizeof( struct collresolve_ConfObject ), /* tp_basicsize */
	0,                      /* tp_itemsize */
	(destructor)python_collresolve_conf_dealloc, /* tp_dealloc */
	0,                      /* tp_print */
	0,                      /* tp_getattr */
	0,                      /* tp_setattr */
	0,                      /* tp_reserved */
	0,                      /* tp_repr */
	0,                      /* tp_as_number */
	0,                      /* tp_as_sequence */
	0,                      /* tp_as_mapping */
	0,                      /* tp_hash  */
	0,                      /* tp_call */
	0,                      /* tp_str */
	0,                      /* tp_getattro */
	0,                      /* tp_setattro */
	0,                      /* tp_as_buffer */
	Py_TPFLAGS_DEFAULT,     /* tp_flags */
	"Collision configuration object", /* tp_doc */
};

static PyObject* python_collresolve_conf_new_internal( void ) {
	PyObject* self = collresolve_ConfType.tp_alloc( &collresolve_ConfType, 0 );

	if ( self != NULL ) {
		struct collresolve_ConfObject* pyconf = (struct collresolve_ConfObject*) self;
		pyconf->data = collresolve_conf_new();
	}

	return self;
}

static PyObject* python_collresolve_conf_new( PyTypeObject* type, PyObject* args, PyObject* kwds ) {
	return python_collresolve_conf_new_internal();
}

/**
 * Module methods
 */

static PyObject* python_collresolve_model_desc( PyObject* self, PyObject* args ) {
	int model;

	if ( !PyArg_ParseTuple( args, "i", &model ) ) {
		return NULL;
	}

	const char* desc = collresolve_model_desc( model );

	if ( desc == NULL ) {
		PyErr_SetString( collresolve_Error, "argument is an invalid model code" );
		return NULL;
	} else {
		return Py_BuildValue( "s", desc );
	}
}

static PyObject* python_collresolve_regime_desc( PyObject* self, PyObject* args ) {
	int regime;

	if ( !PyArg_ParseTuple( args, "i", &regime ) ) {
		return NULL;
	}

	const char* desc = collresolve_regime_desc( regime );

	if ( desc == NULL ) {
		PyErr_SetString( collresolve_Error, "argument is an invalid regime code" );
		return NULL;
	} else {
		return Py_BuildValue( "s", desc );
	}
}

static PyObject* python_collresolve_error_desc( PyObject* self, PyObject* args ) {
	int code;

	if ( !PyArg_ParseTuple( args, "i", &code ) ) {
		return NULL;
	}

	const char* desc = collresolve_error_message( code );

	if ( desc == NULL ) {
		PyErr_SetString( collresolve_Error, "argument is an invalid error code" );
		return NULL;
	} else {
		return Py_BuildValue( "s", desc );
	}
}

static PyObject* python_collresolve_common_conf( PyObject* self, PyObject* args, int ( *func )( struct collresolve_conf* ) ) {
	PyObject* pyconfobj = NULL;

	if ( self == NULL || PyModule_CheckExact( self ) ) {
		if ( !PyArg_ParseTuple( args, "|O", &pyconfobj ) ) {
			return NULL;
		}
	} else {
		if ( !PyArg_ParseTuple( args, "" ) ) {
			return NULL;
		}
		pyconfobj = self;
	}

	if ( pyconfobj == NULL ) {
		pyconfobj = python_collresolve_conf_new_internal();
	} else {
		if ( !PyObject_TypeCheck( pyconfobj, &collresolve_ConfType ) ) {
			PyErr_SetString( PyExc_TypeError, "arg #1 not a collresolve.Conf object" );
			return NULL;
		}
		// We will also return it, so reference count must be incremented
		Py_INCREF( pyconfobj );
	}

	if ( pyconfobj != NULL ) {
		struct collresolve_ConfObject* pyconf = (struct collresolve_ConfObject*) pyconfobj;
		int res = func( pyconf->data );
		if ( res < 0 ) {
			Py_DECREF( pyconfobj );
			PyErr_SetString( collresolve_Error, collresolve_error_message( res ) );
			return NULL;
		}
	}

	return pyconfobj;
}

static PyObject* python_collresolve_conf_unit_si( PyObject* self, PyObject* args ) {
	return python_collresolve_common_conf( self, args, collresolve_conf_unit_si );
}

static PyObject* python_collresolve_conf_unit_msun_au_day( PyObject* self, PyObject* args ) {
	return python_collresolve_common_conf( self, args, collresolve_conf_unit_msun_au_day );
}

static PyObject* python_collresolve_conf_unit_m_earth( PyObject* self, PyObject* args ) {
	return python_collresolve_common_conf( self, args, collresolve_conf_unit_m_earth );
}

static PyObject* python_collresolve_conf_unit_merc( PyObject* self, PyObject* args ) {
	return python_collresolve_common_conf( self, args, collresolve_conf_unit_merc );
}

static PyObject* python_collresolve_conf_model( PyObject* self, PyObject* args ) {
	PyObject* pyconfobj = NULL;
	int model;

	if ( self == NULL || PyModule_CheckExact( self ) ) {
		if ( !PyArg_ParseTuple( args, "Oi", &pyconfobj, &model ) ) {
			return NULL;
		}
	} else {
		if ( !PyArg_ParseTuple( args, "i", &model ) ) {
			return NULL;
		}
		pyconfobj = self;
	}

	if ( !PyObject_TypeCheck( pyconfobj, &collresolve_ConfType ) ) {
		PyErr_SetString( PyExc_TypeError, "arg #1 not a collresolve.Conf object" );
		return NULL;
	}

	struct collresolve_ConfObject* pyconf = (struct collresolve_ConfObject*) pyconfobj;
	int res = collresolve_conf_model( pyconf->data, model );
	if ( res < 0 ) {
		PyErr_SetString( collresolve_Error, collresolve_error_message( res ) );
		return NULL;
	}

	return Py_BuildValue( "" );
}

static PyObject* python_collresolve_conf_sep_after( PyObject* self, PyObject* args ) {
	PyObject* pyconfobj = NULL;
	double drel;

	if ( self == NULL || PyModule_CheckExact( self ) ) {
		if ( !PyArg_ParseTuple( args, "Od", &pyconfobj, &drel ) ) {
			return NULL;
		}
	} else {
		if ( !PyArg_ParseTuple( args, "d", &drel ) ) {
			return NULL;
		}
		pyconfobj = self;
	}

	if ( !PyObject_TypeCheck( pyconfobj, &collresolve_ConfType ) ) {
		PyErr_SetString( PyExc_TypeError, "arg #1 not a collresolve.Conf object" );
		return NULL;
	}

	struct collresolve_ConfObject* pyconf = (struct collresolve_ConfObject*) pyconfobj;
	int res = collresolve_conf_sep_after( pyconf->data, drel );
	if ( res < 0 ) {
		PyErr_SetString( collresolve_Error, collresolve_error_message( res ) );
		return NULL;
	}

	return Py_BuildValue( "" );
}

static PyObject* python_collresolve_bulk_density( PyObject* self, PyObject* args ) {
	PyObject* pyconfobj = NULL;
	double mass;

	if ( !PyArg_ParseTuple( args, "Od", &pyconfobj, &mass ) ) {
		return NULL;
	}

	if ( !PyObject_TypeCheck( pyconfobj, &collresolve_ConfType ) ) {
		PyErr_SetString( PyExc_TypeError, "arg #1 not a collresolve.Conf object" );
		return NULL;
	}

	struct collresolve_ConfObject* pyconf = (struct collresolve_ConfObject*) pyconfobj;

	double res = collresolve_bulk_density( pyconf->data, mass );
	if ( res < 0. ) {
		int code = (int)( res - 0.5 );
		PyErr_SetString( collresolve_Error, collresolve_error_message( code ) );
		return NULL;
	}

	return Py_BuildValue( "d", res );
}

static PyObject* python_collresolve_body_radius( PyObject* self, PyObject* args ) {
	PyObject* pyconfobj = NULL;
	PyObject* pybodyobj = NULL;

	if ( !PyArg_ParseTuple( args, "OO", &pyconfobj, &pybodyobj ) ) {
		return NULL;
	}

	if ( !PyObject_TypeCheck( pyconfobj, &collresolve_ConfType ) ) {
		PyErr_SetString( PyExc_TypeError, "arg #1 not a collresolve.Conf object" );
		return NULL;
	}

	if ( !PyObject_TypeCheck( pybodyobj, &collresolve_BodyType ) ) {
		PyErr_SetString( PyExc_TypeError, "arg #2 not a collresolve.Body object" );
		return NULL;
	}

	struct collresolve_ConfObject* pyconf = (struct collresolve_ConfObject*) pyconfobj;
	struct collresolve_BodyObject* pybody = (struct collresolve_BodyObject*) pybodyobj;

	int res = collresolve_body_radius( pyconf->data, &( pybody->data ) );

	if ( res < 0 ) {
		PyErr_SetString( collresolve_Error, collresolve_error_message( res ) );
		return NULL;
	}

	return Py_BuildValue( "" );
}

static PyObject* python_collresolve_setup( PyObject* self, PyObject* args, PyObject* keywds ) {
	PyObject* pyconfobj = NULL;
	PyObject* pybigobj = NULL;
	PyObject* pysmallobj = NULL;
	double velocity;
	double angle;

	static char* kwlist[] = { "conf", "big", "small", "velocity", "angle", NULL };

	if ( !PyArg_ParseTupleAndKeywords( args, keywds, "OOOdd", kwlist, &pyconfobj, &pybigobj, &pysmallobj, &velocity, &angle ) ) {
		return NULL;
	}

	if ( !PyObject_TypeCheck( pyconfobj, &collresolve_ConfType ) ) {
		PyErr_SetString( PyExc_TypeError, "arg #1 not a collresolve.Conf object" );
		return NULL;
	}

	if ( !PyObject_TypeCheck( pybigobj, &collresolve_BodyType ) ) {
		PyErr_SetString( PyExc_TypeError, "arg #2 not a collresolve.Body object" );
		return NULL;
	}

	if ( !PyObject_TypeCheck( pysmallobj, &collresolve_BodyType ) ) {
		PyErr_SetString( PyExc_TypeError, "arg #3 not a collresolve.Body object" );
		return NULL;
	}

	struct collresolve_conf* conf = ( ( struct collresolve_ConfObject* ) pyconfobj )->data;
	struct collresolve_body* big = &( ( struct collresolve_BodyObject* ) pybigobj )->data;
	struct collresolve_body* small = &( ( struct collresolve_BodyObject* ) pysmallobj )->data;

	collresolve_setup( conf, big, small, velocity, angle );

	return Py_BuildValue( "" );
}

static PyObject* python_collresolve_setup_dist( PyObject* self, PyObject* args, PyObject* keywds ) {
	PyObject* pyconfobj = NULL;
	PyObject* pybigobj = NULL;
	PyObject* pysmallobj = NULL;
	double distance;
	double velocity;
	double angle;

	static char* kwlist[] = { "conf", "big", "small", "distance", "velocity", "angle", NULL };

	if ( !PyArg_ParseTupleAndKeywords( args, keywds, "OOOddd", kwlist, &pyconfobj, &pybigobj, &pysmallobj, &distance, &velocity, &angle ) ) {
		return NULL;
	}

	if ( !PyObject_TypeCheck( pyconfobj, &collresolve_ConfType ) ) {
		PyErr_SetString( PyExc_TypeError, "arg #1 not a collresolve.Conf object" );
		return NULL;
	}

	if ( !PyObject_TypeCheck( pybigobj, &collresolve_BodyType ) ) {
		PyErr_SetString( PyExc_TypeError, "arg #2 not a collresolve.Body object" );
		return NULL;
	}

	if ( !PyObject_TypeCheck( pysmallobj, &collresolve_BodyType ) ) {
		PyErr_SetString( PyExc_TypeError, "arg #3 not a collresolve.Body object" );
		return NULL;
	}

	struct collresolve_conf* conf = ( ( struct collresolve_ConfObject* ) pyconfobj )->data;
	struct collresolve_body* big = &( ( struct collresolve_BodyObject* ) pybigobj )->data;
	struct collresolve_body* small = &( ( struct collresolve_BodyObject* ) pysmallobj )->data;

	collresolve_setup_dist( conf, big, small, distance, velocity, angle );

	return Py_BuildValue( "" );
}

static PyObject* python_collresolve_common_quant( PyObject* self, PyObject* args, PyObject* keywds, double ( *func )( struct collresolve_conf*, struct collresolve_body, struct collresolve_body ), int neg_err ) {
	PyObject* pyconfobj = NULL;
	PyObject* pybigobj = NULL;
	PyObject* pysmallobj = NULL;

	static char* kwlist[] = { "conf", "big", "small", NULL };

	if ( !PyArg_ParseTupleAndKeywords( args, keywds, "OOO", kwlist, &pyconfobj, &pybigobj, &pysmallobj ) ) {
		return NULL;
	}

	if ( !PyObject_TypeCheck( pyconfobj, &collresolve_ConfType ) ) {
		PyErr_SetString( PyExc_TypeError, "arg #1 not a collresolve.Conf object" );
		return NULL;
	}

	if ( !PyObject_TypeCheck( pybigobj, &collresolve_BodyType ) ) {
		PyErr_SetString( PyExc_TypeError, "arg #2 not a collresolve.Body object" );
		return NULL;
	}

	if ( !PyObject_TypeCheck( pysmallobj, &collresolve_BodyType ) ) {
		PyErr_SetString( PyExc_TypeError, "arg #3 not a collresolve.Body object" );
		return NULL;
	}

	struct collresolve_conf* conf = ( ( struct collresolve_ConfObject* ) pyconfobj )->data;
	struct collresolve_body big = ( ( struct collresolve_BodyObject* ) pybigobj )->data;
	struct collresolve_body small = ( ( struct collresolve_BodyObject* ) pysmallobj )->data;

	double res = func( conf, big, small );

	if ( neg_err && res < 0. ) {
		int code = (int)( res - 0.5 );
		PyErr_SetString( collresolve_Error, collresolve_error_message( code ) );
	}

	return Py_BuildValue( "d", res );
}

static PyObject* python_collresolve_specific_energy( PyObject* self, PyObject* args, PyObject* keywds ) {
	return python_collresolve_common_quant( self, args, keywds, collresolve_specific_energy, 1 );
}

static PyObject* python_collresolve_impact_distance( PyObject* self, PyObject* args, PyObject* keywds ) {
	return python_collresolve_common_quant( self, args, keywds, collresolve_impact_distance, 1 );
}

static PyObject* python_collresolve_impact_velocity( PyObject* self, PyObject* args, PyObject* keywds ) {
	return python_collresolve_common_quant( self, args, keywds, collresolve_impact_velocity, 1 );
}

static PyObject* python_collresolve_escape_velocity( PyObject* self, PyObject* args, PyObject* keywds ) {
	return python_collresolve_common_quant( self, args, keywds, collresolve_escape_velocity, 1 );
}

static PyObject* python_collresolve_infinity_velocity( PyObject* self, PyObject* args, PyObject* keywds ) {
	return python_collresolve_common_quant( self, args, keywds, collresolve_infinity_velocity, 0 );
}

static PyObject* python_collresolve_impact_parameter( PyObject* self, PyObject* args, PyObject* keywds ) {
	return python_collresolve_common_quant( self, args, keywds, collresolve_impact_parameter, 1 );
}

static PyObject* python_collresolve_impact_angle( PyObject* self, PyObject* args, PyObject* keywds ) {
	return python_collresolve_common_quant( self, args, keywds, collresolve_impact_angle, 1 );
}

/**
 * God himself.
 */

static PyObject* python_collresolve_resolve( PyObject* self, PyObject* args, PyObject* keywds ) {
	PyObject* pyconfobj = NULL;
	PyObject* pybigobj = NULL;
	PyObject* pysmallobj = NULL;
	int n;
	int r = 0;

	static char* kwlist[] = { "conf", "big", "small", "n", "r", NULL };

	if ( !PyArg_ParseTupleAndKeywords( args, keywds, "OOOi|i", kwlist, &pyconfobj, &pybigobj, &pysmallobj, &n, &r ) ) {
		return NULL;
	}

	if ( !PyObject_TypeCheck( pyconfobj, &collresolve_ConfType ) ) {
		PyErr_SetString( PyExc_TypeError, "arg #1 not a collresolve.Conf object" );
		return NULL;
	}

	if ( !PyObject_TypeCheck( pybigobj, &collresolve_BodyType ) ) {
		PyErr_SetString( PyExc_TypeError, "arg #2 not a collresolve.Body object" );
		return NULL;
	}

	if ( !PyObject_TypeCheck( pysmallobj, &collresolve_BodyType ) ) {
		PyErr_SetString( PyExc_TypeError, "arg #3 not a collresolve.Body object" );
		return NULL;
	}

	struct collresolve_conf* conf = ( ( struct collresolve_ConfObject* ) pyconfobj )->data;
	struct collresolve_body big = ( ( struct collresolve_BodyObject* ) pybigobj )->data;
	struct collresolve_body small = ( ( struct collresolve_BodyObject* ) pysmallobj )->data;

	struct collresolve_body* res = malloc( sizeof( struct collresolve_body ) * ( n + 1 ) );

	int status = collresolve_resolve( conf, big, small, n, res );
	if ( status < 0 ) {
		free( res );
		PyErr_SetString( collresolve_Error, collresolve_error_message( status ) );
		return NULL;
	}

	PyObject* ret = PyList_New( n + 1 );

	for ( int i = 0; i <= n; i++ ) {
		PyObject* entry = collresolve_BodyType.tp_alloc( &collresolve_BodyType, 0 );
		memcpy( (void*) entry + offsetof( struct collresolve_BodyObject, data ), ( void* )( res + i ), sizeof( struct collresolve_body ) );
		PyList_SET_ITEM( ret, i, entry );
	}

	free( res );

	if ( r == 0 ) {
		return ret;
	} else {
		return Py_BuildValue( "(Oi)", ret, status );
	}
}

static PyMethodDef collresolve_methods[] = {
	{ "model_desc", python_collresolve_model_desc, METH_VARARGS, "Get a human-readable description of the model code." },
	{ "regime_desc", python_collresolve_regime_desc, METH_VARARGS, "Get a human-readable description of the model code." },
	{ "error_desc", python_collresolve_error_desc, METH_VARARGS, "Get a human-readable description of the error code." },
	{ "conf_unit_si", python_collresolve_conf_unit_si, METH_VARARGS, "Get a configuration object for SI units." },
	{ "conf_unit_msun_au_day", python_collresolve_conf_unit_msun_au_day, METH_VARARGS, "Get a configuration object for Solar mass, AU and day using their IAU definitions." },
	{ "conf_unit_m_earth", python_collresolve_conf_unit_m_earth, METH_VARARGS, "Get a configuration object for Earth mass and SI units." },
	{ "conf_unit_merc", python_collresolve_conf_unit_merc, METH_VARARGS, "Get a configuration object for Mercury units." },
	{ "conf_model", python_collresolve_conf_model, METH_VARARGS, "Set the collision model to use." },
	{ "conf_sep_after", python_collresolve_conf_sep_after, METH_VARARGS, "Set the relative distance factor for spacing after a collision." },
	{ "bulk_density", python_collresolve_bulk_density, METH_VARARGS, "Retrieve the bulk density from a mass." },
	{ "body_radius", python_collresolve_body_radius, METH_VARARGS, "Set a consistent radius from the body properties." },
	{ "setup", (PyCFunction)python_collresolve_setup, METH_VARARGS | METH_KEYWORDS, "Setup collision geometry." },
	{ "setup_dist", (PyCFunction)python_collresolve_setup_dist, METH_VARARGS | METH_KEYWORDS, "Setup collision geometry (with custom initial distance)." },
	{ "specific_energy", (PyCFunction)python_collresolve_specific_energy, METH_VARARGS | METH_KEYWORDS, "Compute collision specific energy." },
	{ "impact_distance", (PyCFunction)python_collresolve_impact_distance, METH_VARARGS | METH_KEYWORDS, "Compute relative distance at initial contact." },
	{ "impact_velocity", (PyCFunction)python_collresolve_impact_velocity, METH_VARARGS | METH_KEYWORDS, "Compute relative velocity at initial contact." },
	{ "escape_velocity", (PyCFunction)python_collresolve_escape_velocity, METH_VARARGS | METH_KEYWORDS, "Compute the mutual escape velocity." },
	{ "infinity_velocity", (PyCFunction)python_collresolve_infinity_velocity, METH_VARARGS | METH_KEYWORDS, "Compute relative velocity at infinity." },
	{ "impact_parameter", (PyCFunction)python_collresolve_impact_parameter, METH_VARARGS | METH_KEYWORDS, "Compute collision specific energy." },
	{ "impact_angle", (PyCFunction)python_collresolve_impact_angle, METH_VARARGS | METH_KEYWORDS, "Compute collision specific energy." },
	{ "resolve", (PyCFunction)python_collresolve_resolve, METH_VARARGS | METH_KEYWORDS, "Compute collisional outcome." },
	{ NULL, NULL, 0, NULL }
};

#if PY_MAJOR_VERSION >= 3
/**
 * Module definition
 */

static struct PyModuleDef collresolvemodule = {
	PyModuleDef_HEAD_INIT,
	"collresolve",  /* name of module */
	NULL, /* module documentation, may be NULL */
	-1,   /* size of per-interpreter state of the module, or -1 if the module keeps state in global variables. */
	collresolve_methods
};
#endif

/**
 * Module initialisation
 */

#if PY_MAJOR_VERSION >= 3
PyMODINIT_FUNC PyInit_collresolve( void )
#else
PyMODINIT_FUNC initcollresolve( void )
#endif
{
	collresolve_BodyType.tp_new = PyType_GenericNew;
	collresolve_ConfType.tp_new = python_collresolve_conf_new;

	if ( PyType_Ready( &collresolve_BodyType ) < 0 ) {
#if PY_MAJOR_VERSION >= 3
		return NULL;
#else
		return;
#endif
	}

	if ( PyType_Ready( &collresolve_ConfType ) < 0 ) {
#if PY_MAJOR_VERSION >= 3
		return NULL;
#else
		return;
#endif
	}

#if PY_MAJOR_VERSION >= 3
	PyObject* module = PyModule_Create( &collresolvemodule );
#else
	PyObject* module = Py_InitModule( "collresolve", collresolve_methods );
#endif

	if ( module == NULL ) {
#if PY_MAJOR_VERSION >= 3
		return NULL;
#else
		return;
#endif
	}

	collresolve_Error = PyErr_NewException( "collresolve.Error", NULL, NULL );

	Py_INCREF( &collresolve_BodyType );
	Py_INCREF( &collresolve_ConfType );
	Py_INCREF( collresolve_Error );

	PyModule_AddObject( module, "Body", (PyObject*) &collresolve_BodyType );
	PyModule_AddObject( module, "Conf", (PyObject*) &collresolve_ConfType );
	PyModule_AddObject( module, "Error", collresolve_Error );

	PyModule_AddIntConstant( module, "MODEL_NONE", COLLRESOLVE_MODEL_NONE );
	PyModule_AddIntConstant( module, "MODEL_PERFECT_MERGE", COLLRESOLVE_MODEL_PERFECT_MERGE );
	PyModule_AddIntConstant( module, "MODEL_LS2012", COLLRESOLVE_MODEL_LS2012 );
	PyModule_AddIntConstant( module, "MODEL_SL2012", COLLRESOLVE_MODEL_SL2012 );
	PyModule_AddIntConstant( module, "MODEL_C2019", COLLRESOLVE_MODEL_C2019 );

	PyModule_AddIntConstant( module, "REGIME_MERGE", COLLRESOLVE_REGIME_MERGE );
	PyModule_AddIntConstant( module, "REGIME_DISRUPTION", COLLRESOLVE_REGIME_DISRUPTION );
	PyModule_AddIntConstant( module, "REGIME_SUPERCATASTROPHIC", COLLRESOLVE_REGIME_SUPERCATASTROPHIC );
	PyModule_AddIntConstant( module, "REGIME_GRAZE_AND_MERGE", COLLRESOLVE_REGIME_GRAZE_AND_MERGE );
	PyModule_AddIntConstant( module, "REGIME_HIT_AND_RUN", COLLRESOLVE_REGIME_HIT_AND_RUN );

	PyModule_AddIntConstant( module, "ERROR_GENERAL", COLLRESOLVE_ERROR_GENERAL );
	PyModule_AddIntConstant( module, "ERROR_NO_CONF", COLLRESOLVE_ERROR_NO_CONF );
	PyModule_AddIntConstant( module, "ERROR_INCORRECT_PARAMETER", COLLRESOLVE_ERROR_INCORRECT_PARAMETER );
	PyModule_AddIntConstant( module, "ERROR_INCORRECT_MODEL", COLLRESOLVE_ERROR_INCORRECT_MODEL );
	PyModule_AddIntConstant( module, "ERROR_NON_CROSSING", COLLRESOLVE_ERROR_NON_CROSSING );

#if PY_MAJOR_VERSION >= 3
	return module;
#endif
}
