import distutils.core

module = distutils.core.Extension(
	"collresolve",
	sources = [ "collresolve.c", "collresolve_python.c",
			"cambioni2019/accretion_efficiency.c", "cambioni2019/orbital_hnr.c", "cambioni2019/collision_classifier.c" ]
)

distutils.core.setup(
	name = "collresolve",
	version = "1.1",
	description = "Analyse and predict outcomes of collision in N-body",
	ext_modules = [ module ]
)
