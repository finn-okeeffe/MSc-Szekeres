# Compiler
FC=gfortran

# Healpix location
HP=/usr/local/Healpix_3.30
# cfitsio location
CF=/usr/local

# flags to include healpix, cfitsio, and openmp
HPFLAGS=-I${HP}/include-gfortran \
		-L${HP}/lib-gfortran \
		-I${CF}/include \
		-L${CF}/lib \
		-lhealpix -lcfitsio #-lhpxgif -lsharp

# compiler flags (standards, warnings, defaults, optimizations, etc)
FFLAGS= -O3 -Wall -Wextra -std=f2008 -fdefault-real-8 -finit-real=NaN -fcheck=all \
		-fopenmp -fPIC -fbacktrace

# source files
LIBS=solveode.f90 szmodel.f90 initialize_k.f90 \
		ray_propagator.f90 CMB_rays.f90 trace_composite.f90 model_statistics.f90
SRC=${LIBS} main.f90
OBJ=${SRC:.f90=.o}

# fortran file for python library
PY_LIB=py_szekeres.f90

PY_IL_FLAGS = -I${HP}/include-gfortran \
	      -L${HP}/lib-gfortran \
		-I${CF}/include\
		-L${CF}/lib \
		-lhealpix -lcfitsio -lgomp #-lhpxgif -lsharp


%.o: %.f90
	${FC} -c ${FFLAGS} -o $@ $< ${HPFLAGS}

main: ${OBJ}
	${FC} ${FFLAGS} -o $@ ${OBJ} ${HPFLAGS}

clean:
	@rm *.mod *.o

main.o: ${LIBS:.f90=.o}

# build python library
pylib: main
	f2py --f90flags="${FFLAGS}" ${PY_IL_FLAGS}\
	 -c ${LIBS:.f90=.o} ${PY_LIB} -m szekeres
