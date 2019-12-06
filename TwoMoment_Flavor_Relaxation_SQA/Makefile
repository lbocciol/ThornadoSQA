MICROPHYSICS = WEAKLIB
MOMENT_CLOSURE = MAXIMUM_ENTROPY_CB

THORNADO_DIR ?= ../../../
include $(THORNADO_DIR)/Build/Makefile_Build
include $(THORNADO_DIR)/Build/Makefile_Thornado_Dependencies

WEAKLIB_DIR ?= $(HOME)/weaklib
include $(WEAKLIB_DIR)/Distributions/Build/Makefile_Path
include $(WEAKLIB_DIR)/Distributions/Build/Makefile_WeakLib_ObjectFiles
include $(WEAKLIB_DIR)/Distributions/Build/Makefile_WeakLib_Dependencies


.DEFAULT_GOAL := all

all: Relaxation_SQA

Relaxation_SQA: \
	$(weaklib) \
	$(thornado) \
	InterpolationModule.o \
	ReadProfileModule.o \
	InitializationModule.o \
	InputOutputRelaxationModule.o \
	ImplicitSolverModule.o \
	IntegrationModule.o \
	OscillationsUtilsModule.o \
    OscillationsModule.o \
    ThornadoSQAInterfaceModule.o \
	Relaxation_SQA.o
	$(FLINKER) -fPIC $(FLAGS) -o Relaxation_SQA_$(MACHINE) \
	$(weaklib) \
	$(thornado) \
	InterpolationModule.o \
	ReadProfileModule.o \
	InitializationModule.o \
	InputOutputRelaxationModule.o \
	ImplicitSolverModule.o \
	IntegrationModule.o \
	OscillationsUtilsModule.o \
    OscillationsModule.o \
    ThornadoSQAInterfaceModule.o \
	Relaxation_SQA.o \
	$(LIBRARIES)

clean:
	rm -f *.o *.mod *.ld

clobber: clean
	rm -f Relaxation_SQA_$(MACHINE)

Relaxation_SQA.o: \
  Relaxation_SQA.f90
