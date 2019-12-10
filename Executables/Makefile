include Makefile_ObjectFiles

default:
	f2py -c --f90exec=mpif90 $(F90_MODULES) \
	--f90flags="-fopenmp -g -fbacktrace -ffpe-trap=invalid,zero -ffpe-summary=invalid,zero -fbounds-check" -lgomp \
	-m Thor2SQA ThornadoSQAInterfaceModule.f90 ProgramStartEndModule.f90 
clean:
	rm -f *.so