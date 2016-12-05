## use ifort and f2py at NERSC
ifeq (${NERSC_HOST}, edison)
        FF = gfortran
        FPY = f2py
        F90FLAGS = "-fopenmp"
        OPT = --opt=-O3 -lgfortran -lifcore
else ifeq (${NERSC_HOST}, cori)
        FF = gfortran
        FPY = f2py
        F90FLAGS = "-fopenmp"
        OPT = --opt=-O3 -lgfortran -lifcore
else
        FF = gfortran
        FPY = f2py-2.7
        F90FLAGS = "-fopenmp"
        OPT = --opt=-ffixed-line-length-none --opt=-O3
endif

all: biases

biases:
		${FPY} --fcompiler=${FF} --f90flags=${F90FLAGS} -lgomp -c LensingBiases.f90 -m LensingBiases_f ${OPT}

clean:
		-rm *.so
