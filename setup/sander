from distutils.core import setup, Extension

amber = "/Users/smarti/Devel/amber/amber20_src" 

setup( 
	name = "Sander hooks", 
	ext_modules = [ 
		Extension( "_sander",
			sources = [ "qm3/engines/sander.c" ],
			include_dirs = [ amber + "/include" ],
			library_dirs = [ amber + "/lib" ],
			libraries = [ "sander" ]
		)
	]
)
