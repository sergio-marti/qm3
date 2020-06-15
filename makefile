include config


all: matrix minimize cions volume colvar_v mmint mol_mech ode fitpack grids dynamo shm conn


cions:
	$(PYX) setup/cions build_ext --build-lib qm3/utils


mpi:
	$(MPI) $(PYX) setup/mpi build_ext --build-lib qm3/utils


shm:
	$(PYX) setup/shm build_ext --build-lib qm3/utils


volume:
	$(PYX) setup/volume build_ext --build-lib qm3/utils


conn:
	$(PYX) setup/conn build_ext --build-lib qm3/utils


grids:
	$(PYX) setup/grids build_ext --build-lib qm3/maths


plumed:
	$(PYX) setup/plumed build_ext --build-lib qm3/engines


sander:
	$(PYX) setup/sander build_ext --build-lib qm3/engines


dynamo:
	$(PYX) setup/dynamo build_ext --build-lib qm3/engines


colvar_v:
	$(PYX) setup/colvar_v build_ext --build-lib qm3/engines


long_elec:
	$(PYX) setup/long_elec build_ext --build-lib qm3/engines


mmint:
	$(PYX) setup/mmint build_ext --build-lib qm3/engines


mol_mech:
	$(PYX) setup/mol_mech build_ext --build-lib qm3/engines


ode:
	$(PYX) setup/ode build_ext --build-lib qm3/maths


fitpack: fitpack.o
	$(PYX) setup/fitpack build
	gcc -lpthread $(SHD) -o qm3/maths/_fitpack.so \
		fitpack.o \
		`find build -name fitpack.o` \
		-lm $(PYL)


matrix: ilaenv_fix.o lapack_deps.o
	CFLAGS='-DMCORE=\"$(MTH)\"' $(PYX) setup/matrix build
	gcc $(SHD) -o qm3/maths/_matrix.so \
		ilaenv_fix.o \
		`find build -name matrix.o` \
		$(MLB) -lm $(PYL)


minimize: cgplus.o
	$(PYX) setup/minimize build
	gcc -lpthread $(SHD) -o qm3/actions/_minimize.so \
		cgplus.o \
		`find build -name minimize.o` \
		-lm $(PYL)


clean:
	rm -rf *.o build
	
	
cgplus.o:
	gfortran -c -w -O1 -fPIC qm3/actions/cgplus.f


fitpack.o:
	gfortran -c -w -O1 -fPIC qm3/maths/fitpack.f


ilaenv_fix.o:
	gfortran -c -w -O1 -fPIC qm3/maths/ilaenv_fix.f

lapack_deps.o:
	gfortran -c -w -O2 -fPIC qm3/maths/lapack_deps.f


dftb.so:
	gfortran $(SHD) -o dftb.so qm3/engines/dftb.f90 \
		-I/Users/smarti/Devel/dftb+/dftbplus/_build/prog/dftb+/include \
		/Users/smarti/Devel/dftb+/dftbplus/_build/prog/dftb+/libdftbplus.a \
		/Users/smarti/Devel/dftb+/dftbplus/_build/external/dftd3/origin/lib/libdftd3.a \
		-framework Accelerate -fopenmp


xtb.so:
	gfortran $(SHD) -o xtb.so qm3/engines/xtb.f90 \
		-I/Users/smarti/Devel/xtb/dist/xtb-mod \
		/Users/smarti/Devel/xtb/dist/libxtb.a \
		-framework Accelerate


dftd3.so:
	gfortran $(SHD) -o dftd3.so qm3/engines/dftd3.f90 \
		-I/Users/smarti/Devel/dftd3/3.2r0/dftd3-lib-0.9/lib \
		/Users/smarti/Devel/dftd3/3.2r0/dftd3-lib-0.9/lib/libdftd3.a


sqm.so:
	gfortran $(SHD) -o sqm.so qm3/engines/sqm_18.f90 \
		-I/Users/smarti/Devel/amber/18/amber18/AmberTools/src/sqm \
		/Users/smarti/Devel/amber/18/amber18/AmberTools/src/sqm/*.o \
		-L/Users/smarti/Devel/amber/18/amber18/lib -larpack -lxblas-amb \
		/Users/smarti/Devel/amber/18/amber18/AmberTools/src/lib/sys.a \
		-framework Accelerate
