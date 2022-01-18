open:
	gcc prefine.c libopenblas.a -lpthread -lm
	echo "export OMP_NUM_THREADS=1"


macos:
	gcc prefine.c -framework Accelerate


mkl:
	gcc prefine.c -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_rt -lpthread -lm -ldl
	echo "export MKL_NUM_THREADS=1"
