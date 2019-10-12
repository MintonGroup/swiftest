import distutils.core

module = distutils.core.Extension(
	"collresolve",
	sources = [ "collresolve.c", "collresolve_python.c",
			"cambioni2019/accretion_efficiency.c", "cambioni2019/orbital_hnr.c",
			"cambioni2019/bsxfun.c", "cambioni2019/CompactClassificationECOC.c", "cambioni2019/CompactSVM.c", "cambioni2019/collision_classifier.c", "cambioni2019/minOrMax.c", "cambioni2019/Poly.c", "cambioni2019/rtGetInf.c", "cambioni2019/rtGetNaN.c", "cambioni2019/rt_nonfinite.c" ]
)

distutils.core.setup(
	name = "collresolve",
	version = "1.0",
	description = "Analyse and predict outcomes of collision in N-body",
	ext_modules = [ module ]
)
