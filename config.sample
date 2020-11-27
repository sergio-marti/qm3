# ------------------------------------------------------------------
# Python

#PYX = python
#PYL = -lpython

PYX = python3
PYL = -L/usr/local/Cellar/python@3.9/3.9.0_2/Frameworks/Python.framework/Versions/3.9/lib -lpython3.9

#PYL = -L/opt/local/Cellar/python/3.7.1/Frameworks/Python.framework/Versions/3.7/lib -lpython3.7
#PYL = -lpython3.5m

# ------------------------------------------------------------------
# Linker

# Darwin/macOS
SHD = -dynamiclib

# Linux
#SHD = -shared

# ------------------------------------------------------------------
# BLAS/LAPACK

MTH = macos_Accelerate
MLB = -framework Accelerate

#MTH = Lapack
#MLB = lapack_deps.o -L/opt/local/Cellar/gcc/8.2.0/lib/gcc/8 -lgfortran

#MTH = Intel_MKL
#MLB = -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_intel_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -liomp5 -lpthread -lm -ldl
#MLB = -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_rt -lpthread -lm -ldl

#MTH = OpenBLAS
#MLB = libopenblas.a -lpthread

# ------------------------------------------------------------------
# MPI

# Darwin/macOS
MPI = CC=mpicc

# Linux
#MPI = CC=mpicc LDSHARED="mpicc -shared"


# ------------------------------------------------------------------
# External programs
DFTB = /Users/smarti/Devel/dftbplus/_build/
XTB = /Users/smarti/Devel/xtb/dist
DFTD3 = /Users/smarti/Devel/dftd3/3.2r0/dftd3-lib-0.9/lib
DFTD4 = /Users/smarti/Devel/dftd4/build
DFTD4_INC = /Users/smarti/Devel/dftd4/build/dftd4@sta
AMBER = /Users/smarti/Devel/amber/amber20_src
LIO = /home/gpuser/smarti/slio
