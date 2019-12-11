include Makefile_ObjectFiles

default:
	f2py -c --f90exec=mpif90 $(F90_MODULES) \
	--f90flags="-fopenmp -g -O2" -lgomp \
	-m Thor2SQA ThornadoSQAInterfaceModule.f90 ProgramStartEndModule.f90 
clean:
	rm -f *.so