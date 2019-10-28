from distutils.core import setup, Extension

plumed = "/Users/smarti/Devel/plumed2/dist_gnu"

setup( 
	name = "Plumed hooks", 
	ext_modules = [ 
		Extension( "_plumed",
			sources = [ "qm3/engines/plumed.c" ],
			include_dirs = [ plumed + "/include/plumed/wrapper", plumed + "/include" ],
			library_dirs = [ plumed + "/lib" ],
			libraries = [ "plumed" ]
		)
	]
)
