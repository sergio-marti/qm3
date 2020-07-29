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
	gfortran $(SHD) -o dftb.so qm3/engines/dftb.f90 -Jbuild -Ibuild \
		-I$(DFTB)/prog/dftb+/include \
		$(DFTB)/prog/dftb+/libdftbplus.a \
		$(DFTB)/external/dftd3/origin/lib/libdftd3.a \
		$(DFTB)external/mudpack/libmudpack.a $(MLB) -fopenmp


xtb.so:
	gfortran $(SHD) -o xtb.so qm3/engines/xtb.f90 -Jbuild -Ibuild \
		-I$(XTB)/xtb-mod $(XTB)/libxtb.a $(MLB)


dftd3.so:
	gfortran $(SHD) -o dftd3.so qm3/engines/dftd3.f90 -Jbuild -Ibuild \
		-I$(DFTD3) $(DFTD3)/libdftd3.a


dftd4.so:
	gfortran $(SHD) -o dftd4.so qm3/engines/dftd4.f90 -Jbuild -Ibuild \
		-I$(DFTD4_INC) $(DFTD4)/libdftd4.a -fopenmp $(MLB)


sqm.so:
	@gfortran $(SHD) -o sqm.so qm3/engines/sqm.f90 -Jbuild -Ibuild \
		-I$(AMBER)/AmberTools/src/sqm \
		$(AMBER)/AmberTools/src/sqm/file_io_dat.o \
		$(AMBER)/AmberTools/src/sqm/constants.o \
		$(AMBER)/AmberTools/src/sqm/findmask.o \
		$(AMBER)/AmberTools/src/sqm/sqm.SQM.o \
		$(AMBER)/AmberTools/src/sqm/xmin.o \
		$(AMBER)/AmberTools/src/sqm/xminC.o \
		$(AMBER)/AmberTools/src/sqm/amopen.o \
		$(AMBER)/AmberTools/src/sqm/mexit.o \
		$(AMBER)/AmberTools/src/sqm/assert.o \
		$(AMBER)/AmberTools/src/sqm/timer_dummy.o \
		$(AMBER)/AmberTools/src/sqm/nmlsrc.o \
		$(AMBER)/AmberTools/src/sqm/qm_print_info.SQM.o \
		$(AMBER)/AmberTools/src/sqm/qm2_energy.SQM.o \
		$(AMBER)/AmberTools/src/sqm/qm2_read_nm_and_alloc.SQM.o \
		$(AMBER)/AmberTools/src/sqm/qm2_scf.SQM.o \
		$(AMBER)/AmberTools/src/sqm/qm2_dftb_scf.SQM.o \
		$(AMBER)/AmberTools/src/sqm/qm2_allocate_e_repul.o \
		$(AMBER)/AmberTools/src/sqm/qm2_calc_charges.o \
		$(AMBER)/AmberTools/src/sqm/qm2_calc_dipole.o \
		$(AMBER)/AmberTools/src/sqm/qm2_calc_rij_and_eqns.o \
		$(AMBER)/AmberTools/src/sqm/qm2_dihed.o \
		$(AMBER)/AmberTools/src/sqm/qm2_fock.o \
		$(AMBER)/AmberTools/src/sqm/qm2_get_qm_forces.o \
		$(AMBER)/AmberTools/src/sqm/qm2_get_qmmm_forces.o \
		$(AMBER)/AmberTools/src/sqm/qm2_h1elec.o \
		$(AMBER)/AmberTools/src/sqm/qm2_hcore_qmqm.o \
		$(AMBER)/AmberTools/src/sqm/qm2_hcore_qmmm.o \
		$(AMBER)/AmberTools/src/sqm/qm2_identify_peptide_links.o \
		$(AMBER)/AmberTools/src/sqm/qm2_load_params_and_allocate.o \
		$(AMBER)/AmberTools/src/sqm/qm2_repp.o \
		$(AMBER)/AmberTools/src/sqm/qm2_rotate_qmqm.o \
		$(AMBER)/AmberTools/src/sqm/qm2_core_core_repulsion.o \
		$(AMBER)/AmberTools/src/sqm/qm2_core_core_repulsion_dxyz.o \
		$(AMBER)/AmberTools/src/sqm/qm2_iterator_mod.o \
		$(AMBER)/AmberTools/src/sqm/qm2_diagonalizer_module.o \
		$(AMBER)/AmberTools/src/sqm/qm2_setup_orb_exp.o \
		$(AMBER)/AmberTools/src/sqm/qm2_smallest_number.o \
		$(AMBER)/AmberTools/src/sqm/qm2_fock_predict.o \
		$(AMBER)/AmberTools/src/sqm/qm_gb.o \
		$(AMBER)/AmberTools/src/sqm/qm_zero_charges.o \
		$(AMBER)/AmberTools/src/sqm/qm_assign_atom_types.o \
		$(AMBER)/AmberTools/src/sqm/qm_link_atoms.o \
		$(AMBER)/AmberTools/src/sqm/qm2_print_charges.o \
		$(AMBER)/AmberTools/src/sqm/qmmm_qmtheorymodule.o \
		$(AMBER)/AmberTools/src/sqm/qm2_print_bondorders.o \
		$(AMBER)/AmberTools/src/sqm/qm2_pm6_hof_module.o \
		$(AMBER)/AmberTools/src/sqm/qm2_params_module.o \
		$(AMBER)/AmberTools/src/sqm/qm2_dftb_module.o \
		$(AMBER)/AmberTools/src/sqm/qm2_dftb_broyden.o \
		$(AMBER)/AmberTools/src/sqm/qm2_dftb_dispersion_egr.o \
		$(AMBER)/AmberTools/src/sqm/qm2_dftb_dispersion_params.o \
		$(AMBER)/AmberTools/src/sqm/qm2_dftb_dispersionread.o \
		$(AMBER)/AmberTools/src/sqm/qm2_dftb_energy.o \
		$(AMBER)/AmberTools/src/sqm/qm2_dftb_ewevge.o \
		$(AMBER)/AmberTools/src/sqm/qm2_dftb_externalshift.o \
		$(AMBER)/AmberTools/src/sqm/qm2_dftb_fermi.o \
		$(AMBER)/AmberTools/src/sqm/qm2_dftb_get_qm_forces.o \
		$(AMBER)/AmberTools/src/sqm/qm2_dftb_gamma.o \
		$(AMBER)/AmberTools/src/sqm/qm2_dftb_gammamat.o \
		$(AMBER)/AmberTools/src/sqm/qm2_dftb_get_qmmm_forces.o \
		$(AMBER)/AmberTools/src/sqm/qm2_dftb_gettab.o \
		$(AMBER)/AmberTools/src/sqm/qm2_dftb_load_params.o \
		$(AMBER)/AmberTools/src/sqm/qm2_dftb_mulliken.o \
		$(AMBER)/AmberTools/src/sqm/qm2_dftb_repulsiv.o \
		$(AMBER)/AmberTools/src/sqm/qm2_dftb_self.o \
		$(AMBER)/AmberTools/src/sqm/qm2_dftb_shift.o \
		$(AMBER)/AmberTools/src/sqm/qm2_dftb_skpar.o \
		$(AMBER)/AmberTools/src/sqm/qm2_dftb_slkode.o \
		$(AMBER)/AmberTools/src/sqm/qm2_dftb_slktrafo.o \
		$(AMBER)/AmberTools/src/sqm/qm2_dftb_read_cm3.o \
		$(AMBER)/AmberTools/src/sqm/qm2_dftb_cm3.o \
		$(AMBER)/AmberTools/src/sqm/qm2_dftb_gb.o \
		$(AMBER)/AmberTools/src/sqm/qm2_dftb_3rd_order.o \
		$(AMBER)/AmberTools/src/sqm/qmmm_module.o \
		$(AMBER)/AmberTools/src/sqm/qmmm_nml_module.o \
		$(AMBER)/AmberTools/src/sqm/qmmm_struct_module.o \
		$(AMBER)/AmberTools/src/sqm/qmmm_vsolv_module.o \
		$(AMBER)/AmberTools/src/sqm/qm2_parameters.o \
		$(AMBER)/AmberTools/src/sqm/dh_correction_module.o \
		$(AMBER)/AmberTools/src/sqm/utilitiesModule.o \
		$(AMBER)/AmberTools/src/sqm/elementOrbitalIndex.o \
		$(AMBER)/AmberTools/src/sqm/parameterReader.o \
		$(AMBER)/AmberTools/src/sqm/slater_overlap.o \
		$(AMBER)/AmberTools/src/sqm/qm2_h1elec_d.o \
		$(AMBER)/AmberTools/src/sqm/rotation.o \
		$(AMBER)/AmberTools/src/sqm/qm2_repp_d.o \
		$(AMBER)/AmberTools/src/sqm/opnq_Edisp.o \
		$(AMBER)/AmberTools/src/sqm/opnq_Erep.o \
		$(AMBER)/AmberTools/src/sqm/opnq_Evdw.o \
		$(AMBER)/AmberTools/src/sqm/opnq.o \
		$(AMBER)/AmberTools/src/sqm/opnq_SwitchMod.o \
		$(AMBER)/AmberTools/src/sqm/qm2_fock_d.o \
		$(AMBER)/AmberTools/src/sqm/MNDOChargeSeparation.o \
		$(AMBER)/AmberTools/src/sqm/qm2_print_energy.o \
		$(AMBER)/AmberTools/src/sqm/qm2_dftb_get_dftb3_parameters.o \
		$(AMBER)/AmberTools/src/sqm/qm2_dftb_gamma_dftb3.o \
		-L$(AMBER)/lib -larpack -lxblas-amb $(AMBER)/AmberTools/src/lib/sys.a $(MLB)
