# ThornadoSQA

- ImplicitSolverModule: it's where the implicit step is carried out. 
It's a simple backward Euler that only includes Emission and Absorption

- InitializationModule: it's where Fields, Radiation and Oscillations are initialized.

- InputOutputRelaxationModule: dumps fMatrixOscillations and also attempts to create restart files, 
although I couldn't quite succeed, so restart doesn't work yet. 
However, since I was trying to set up restarts, the Restart Subroutine compiles 
only if FileNumber in 
/Modules/InputOutput/InputOutputModuleHDF.f90 is made public. 

- IntegrationModule: this is basically TimeSteppingCASTRO.f90 with my implicit solver instead of the original one 
and without any matter feedback.

- InterpolationModule: there's a 1D linear and quadratic interpolation (for the matter profile) and a 2D linear interpolation 
(in space and energy, for the Radiation). The 2D interpolation routine does a simple linear extrapolation outside the domain 
which is really bad, try not to use and remain within the domain.

- OscillationsModule: the main routine here is EvolveOscillations (i.e. the "evolve_oscillations" function in evolve.h from IsotropicSQA) 
which calls the RK_Step (i.e. the "K" function in misc.h from IsotropicSQA). 
The only difference is that I integrate in time instead that in space, but it's only a "factor of c" difference. 
Also, keep in mind that SMatrixOsc is the "local S", meaning that SMatrixOsc = S(t_n, t_n+1), 
therefore to get the total S one has to do S(t_0,t_n) = S(t_0,t_1) * S(t_1,t_2) * ... * S(t_n-1,t_n)

- OscillationsUtilsModule: just some routines that calculate eigenvalues, eigenvectors and things like that. 
Few things worth mentioning are that in general S = B(Y) * W(Y),
 and that JInverse is used to solve the self-interaction part, but I don't know exactly how.

- ProgramStartEndModule: Used ONLY for the python wrapping. It does the Initialization and the cleanups doen in Relaxation_SQA, 
plus telling python some of the units used in Thornado

- ReadProfileModule: the routines used to read the matter and radiation 1D slices from CHIMERA

- Relaxation_SQA: this is the driver where all the general initializations are performed, and where the main integration loop is.

- ThornadoSQAInterfaceModule: this is used both in the "normal" fortran code and in the python wrapping, it's where the SQA driver is. 
It calls the SQA in the appropriate zones, and also takes care of calculating the total S every time, 
which is then used to initialize the radiation field for the next zone as described in Stapleford (2019). 
At the end, it calculates the opacities as described in that same paper. 
You'll notice that there some subroutines whose only purpose is to interface with python: 
GetShockRadius, Get_Profile, OscillationsInterface, InitializeOscInterface, FinalizeOscInterface. 
