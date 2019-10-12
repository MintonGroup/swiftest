# Example script that shows the Python interface of the collresolve library

from __future__ import print_function

import math
import collresolve

vel = 20.e3
angle = math.radians( 45. )

conf = collresolve.Conf()

collresolve.conf_unit_m_earth( conf )
collresolve.conf_model( conf, 2 )

big = collresolve.Body( mass = 1., radius = 6.371e6 )
small = collresolve.Body( mass = 0.1, radius = 3.4e6 )

collresolve.setup( conf, big, small, vel, angle )

models = [ collresolve.MODEL_PERFECT_MERGE, collresolve.MODEL_LS2012, collresolve.MODEL_SL2012, collresolve.MODEL_C2019 ]

print( "impact velocity = {0:.1f} m/s".format( collresolve.impact_velocity( conf, big, small ) ) )
print( "velocity ratio = {0:.2f}".format( collresolve.impact_velocity( conf, big, small ) / collresolve.escape_velocity( conf, big, small ) ) )
print( "impact angle = {0:.1f} deg".format( math.degrees( collresolve.impact_angle( conf, big, small ) ) ) )

for model in models:
	print()
	print( "resulting bodies using the {0:} model:".format( collresolve.model_desc( model ) ) )

	collresolve.conf_model( conf, model )
	res, regime = collresolve.resolve( conf, big, small, 2, 1 )

	print( "  regime =", collresolve.regime_desc( regime ) )

	tot_m = - big.mass - small.mass
	for body in res:
		tot_m += body.mass
		print( "   ", body )

	print( "  mass change =", tot_m )
