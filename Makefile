## use ifort and f2py at NERSC
ifeq (${NERSC_HOST}, edison)
        N_ONE_MAT = ifort -openmp Biases_n1mat.f90 -o getn1_mat -lgfortran -lifcore
        FF = ifort
        #FF = gfortran
        FPY = f2py
		F90FLAGS = "-fopenmp"
        OPT = --opt=-O3 -lgfortran -lifcore
else ifeq (${NERSC_HOST}, cori)
        N_ONE_MAT = ftn -openmp Biases_n1mat.f90 -o getn1_mat
        FF = ifort
        FF = gfortran
        FPY = f2py
		F90FLAGS = "-fopenmp"
        OPT = --opt=-O3
else
        N_ONE_MAT =
        FF = gfortran
        FPY = f2py-2.7
		F90FLAGS = "-fopenmp"
        OPT = --opt=-ffixed-line-length-none --opt=-O3
endif

all: biases

biases:
		${FPY} --fcompiler=${FF} --f90flags=${F90FLAGS} -lgomp -c LensingBiases.f90 -m LensingBiases ${OPT}

clean:
		-rm *.so
